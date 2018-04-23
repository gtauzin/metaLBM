#pragma once

#include <chrono>

#include "Algorithm.h"
#include "AnalysisList.h"
#include "Commons.h"
#include "Communication.h"
#include "Computation.h"
#include "Distribution.h"
#include "FieldList.h"
#include "FourierDomain.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Options.h"
#include "Transformer.h"
#include "Writer.h"

namespace lbm {

/// Initializes all fields and starts the LBM method
template <class T, Architecture architecture, Implementation implementation>
class Routine {
 protected:
  using Clock = std::chrono::high_resolution_clock;
  using Seconds = std::chrono::duration<double>;

  Communication_ communication;

  FieldWriter_ fieldWriter;
  DistributionWriter_ distributionWriter;
  FieldList<T, architecture> fieldList;
  Distribution<T, architecture> distribution;
  const unsigned int numberElements;

  Curl<double, Architecture::CPU, PartitionningType::OneD, L::dimD, L::dimD>
      curlVelocity;
  ScalarAnalysisList<T, architecture> scalarAnalysisList;
  SpectralAnalysisList<T, architecture> spectralAnalysisList;

  Algorithm<dataT, algorithmT, architecture, implementation> algorithm;
  Stream<architecture> stream;

  double initialMass;
  double finalMass;
  double differenceMass;

  double computationTime;
  double communicationTime;
  double writeFieldTime;
  double writeAnalysisTime;
  double totalTime;

 public:
  Routine(const MathVector<int, 3>& rankMPI_in,
          const MathVector<int, 3>& sizeMPI_in,
          const std::string& processorName_in,
          const unsigned int numberElements_in)
      : communication(rankMPI_in, sizeMPI_in, processorName_in),
        fieldWriter(prefix, rankMPI_in),
        distributionWriter(prefix, rankMPI_in),
        fieldList(rankMPI_in, numberElements_in, fieldWriter),
        curlVelocity(
            fieldList.velocity.getLocalData(),
            fieldList.vorticity.getLocalData(),
            numberElements,
            Cast<unsigned int, ptrdiff_t, 3>::Do(gSD::sLength()).data(),
            gFD::offset(rankMPI_in)),
        distribution(initLocalDistribution<T, architecture>(fieldList.density,
                                                            fieldList.velocity,
                                                            rankMPI_in)),
        numberElements(numberElements_in),
        scalarAnalysisList(fieldList, numberElements, communication),
        spectralAnalysisList(fieldList, numberElements, communication),
        algorithm(fieldList, distribution, numberElements, communication),
        stream(true),
        initialMass(0.0),
        finalMass(0.0),
        differenceMass(0.0),
        computationTime(0.0),
        communicationTime(0.0),
        writeFieldTime(0.0),
        writeAnalysisTime(0.0),
        totalTime(0.0) {}

  void compute() {
    {LBM_INSTRUMENT_ON("Routine<T>::compute", 1)} algorithm.unpack(
        stream);

    Clock::time_point t0;
    Clock::time_point t1;

    if (writeFieldInit) {
      t0 = Clock::now();
      writeFields(startIteration);
      t1 = Clock::now();
      writeFieldTime += Seconds(t1 - t0).count();
    }

    if (writeAnalysisInit) {
      t0 = Clock::now();
      writeAnalyses(startIteration);
      t1 = Clock::now();
      writeAnalysisTime += Seconds(t1 - t0).count();
    }

    printInputs();

    initialMass = communication.reduce(fieldList.density.getLocalData());

    // Execute LBM algorithm
    for (int iteration = startIteration + 1; iteration <= endIteration;
         ++iteration) {
      algorithm.isStored = fieldWriter.getIsWritten(iteration) ||
                           scalarAnalysisList.getIsAnalyzed(iteration) ||
                           spectralAnalysisList.getIsAnalyzed(iteration);

      algorithm.iterate(iteration, stream);

      if (algorithm.isStored) {
        if (writeVorticity)
          curlVelocity.executeSpace();
      }

      t0 = Clock::now();
      writeFields(iteration);
      t1 = Clock::now();
      writeFieldTime += Seconds(t1 - t0).count();

      t0 = Clock::now();
      writeAnalyses(iteration);
      t1 = Clock::now();
      writeAnalysisTime += Seconds(t1 - t0).count();

      communicationTime += algorithm.getCommunicationTime();
      computationTime += algorithm.getComputationTime();
    }

    finalMass = communication.reduce(fieldList.density.getLocalData());

    differenceMass = fabs(initialMass - finalMass) / initialMass;
    totalTime = computationTime + communicationTime + writeFieldTime +
                writeAnalysisTime;
    printOutputs();
  }

 protected:
  void printInputs() {
    if (communication.rankMPI == MathVector<int, 3>({0, 0, 0})) {
      std::cout.precision(15);
      // clang-format off
      std::cout << "-------------------OPTIONS-------------------\n"
                << "Lattice                  : D" << L::dimD << "Q" << L::dimQ << "\n"
                << "Global lengths           : " << gSD::sLength() << "\n"
                << "Global memory            : " << gSD::sVolume() * sizeof(dataT) << "B\n"
                << "----------------------------------------------\n"
                << "NPROCS                   : " << numProcs << "\n"
                << "NTHREADS                 : " << numThreads << "\n"
                << "-------------------PARAMETERS-----------------\n"
                << "Relaxation time          : " << relaxationTime << "\n"
                << "Viscosity                : " << L::cs2 * (relaxationTime - 0.5) << "\n"
                << "Start iteration          : " << startIteration << "\n"
                << "End iteration            : " << endIteration << "\n"
                << "----------------------------------------------\n";
      // clang-format on
    }
  }

  void printOutputs() {
    if (communication.rankMPI == MathVector<int, 3>({0, 0, 0})) {
      std::cout << "-------------------OUTPUTS--------------------\n"
                << "Total time               : " << totalTime << " s\n"
                << "Computatation time       : " << computationTime << " s\n"
                << "Communication time       : " << communicationTime << " s\n"
                << "Analysis time            : " << writeAnalysisTime << " s\n"
                << "Write time               : " << writeFieldTime << " s\n";

      const double mlups = (gSD::sVolume() * 1e-6) /
                           (totalTime / (endIteration - startIteration + 1));

      std::cout << "MLUPS                   : " << mlups << "\n"
                << "Initial mass            : " << initialMass << "\n"
                << "Final mass              : " << finalMass << "\n"
                << "% mass diff.            : " << differenceMass << "\n"
                << "----------------------------------------------\n";
    }
  }

  void writeFields(const unsigned int iteration) {
    { LBM_INSTRUMENT_ON("Routine<T>::writeFields", 2) }

    if (fieldWriter.getIsWritten(iteration)) {
      fieldWriter.openFile(iteration);

      fieldList.writeFields();
      fieldWriter.closeFile();
    }

    if (distributionWriter.getIsBackedUp(iteration)) {
      algorithm.pack(stream);
      distributionWriter.openFile(iteration);
      distributionWriter.writeDistribution(distribution);
      distributionWriter.closeFile();
    }
  }

  void writeAnalyses(const unsigned int iteration) {
    { LBM_INSTRUMENT_ON("Routine<T>::writeFields", 2) }

    if (scalarAnalysisList.getIsAnalyzed(iteration)) {
      scalarAnalysisList.writeAnalyses(iteration);
    }

    if (spectralAnalysisList.getIsAnalyzed(iteration)) {
      spectralAnalysisList.writeAnalyses(iteration);
    }
  }
};

}  // namespace lbm
