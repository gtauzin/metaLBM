#pragma once

#include <chrono>

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "FieldList.h"
#include "Distribution.h"
#include "AnalysisList.h"
#include "Algorithm.h"
#include "Writer.h"
#include "Communication.h"
#include "Computation.h"

namespace lbm {

  template<class T, Architecture architecture, Implementation implementation>
  class Routine {
  protected:
    Communication_ communication;

    FieldWriter_ fieldWriter;
    DistributionWriter_ distributionWriter;
    FieldList<T, architecture> fieldList;
    Distribution<T, architecture> distribution;
    const unsigned int numberElements;

    Curl<double, Architecture::CPU, PartitionningType::OneD,
         L::dimD, L::dimD> curlVelocity;
    ScalarAnalysisList<T, architecture> scalarAnalysisList;
    SpectralAnalysisList<T, architecture> spectralAnalysisList;

    Algorithm<dataT, algorithmT, architecture, implementation> algorithm;
    Computation<Architecture::CPU, L::dimD> computationLocal;

    double initialMass;
    double finalMass;
    double differenceMass;
    double computationTime;
    double communicationTime;
    double writeFieldTime;
    double writeAnalysisTime;
    double totalTime;
    std::chrono::duration<double> dtWriteField;
    std::chrono::duration<double> dtWriteAnalysis;

  public:
    Routine(const MathVector<int, 3>& rankMPI_in,
            const MathVector<int, 3>& sizeMPI_in,
            const std::string& processorName_in,
            const unsigned int numberElements_in)
      : communication(rankMPI_in, sizeMPI_in, processorName_in)
      , fieldWriter(prefix, rankMPI_in)
      , distributionWriter(prefix, rankMPI_in)
      , fieldList(rankMPI_in, numberElements_in, fieldWriter)
      , curlVelocity(fieldList.velocity.getLocalData(), fieldList.vorticity.getLocalData(),
                     numberElements, Cast<unsigned int,
                     ptrdiff_t, 3>::Do(gSD::sLength()).data(),
                     gFD::offset(rankMPI_in))
      , distribution(initLocalDistribution<T, architecture>(fieldList.density,
                                                            fieldList.velocity,
                                                            rankMPI_in))
      , numberElements(numberElements_in)
      , scalarAnalysisList(fieldList, numberElements, communication)
      , spectralAnalysisList(fieldList, numberElements, communication)
      , algorithm(fieldList, distribution, numberElements, communication)
      , computationLocal(lSD::sStart(), lSD::sEnd())
      , initialMass(0.0)
      , finalMass(0.0)
      , differenceMass(0.0)
      , computationTime(0.0)
      , communicationTime(0.0)
      , writeFieldTime(0.0)
      , writeAnalysisTime(0.0)
      , totalTime(0.0)
      , dtWriteField()
      , dtWriteAnalysis()
    {}

    void compute() {
      { INSTRUMENT_ON("Routine<T>::compute",1) }
      algorithm.unpack();

      auto t0 = std::chrono::high_resolution_clock::now();
      auto t1 = std::chrono::high_resolution_clock::now();
      auto t2 = std::chrono::high_resolution_clock::now();

      if(writeFieldInit) {
        t0 = std::chrono::high_resolution_clock::now();
        writeFields(startIteration);
        t1 = std::chrono::high_resolution_clock::now();
        dtWriteField = (t1 - t0);
        writeFieldTime += dtWriteField.count();
      }

      if(writeAnalysisInit) {
        t0 = std::chrono::high_resolution_clock::now();
        writeAnalyses(startIteration);
        t1 = std::chrono::high_resolution_clock::now();
        dtWriteAnalysis = (t1 - t0);
        writeAnalysisTime += dtWriteAnalysis.count();
      }

      printInputs();

      initialMass = communication.reduce(fieldList.density.getLocalData());

      for(int iteration = startIteration+1; iteration <= endIteration; ++iteration) {
        algorithm.isStored = fieldWriter.getIsWritten(iteration)
        || scalarAnalysisList.getIsAnalyzed(iteration)
        || spectralAnalysisList.getIsAnalyzed(iteration);

        algorithm.iterate(iteration);

        if(algorithm.isStored) {
          if(writeVorticity) curlVelocity.executeSpace();
        }

        t0 = std::chrono::high_resolution_clock::now();
        writeFields(iteration);
        t1 = std::chrono::high_resolution_clock::now();
        writeAnalyses(iteration);
        t2 = std::chrono::high_resolution_clock::now();
        dtWriteField = (t1 - t0);
        writeFieldTime += dtWriteField.count();
        dtWriteAnalysis = (t2 - t1);
        writeAnalysisTime += dtWriteAnalysis.count();

        communicationTime += algorithm.getCommunicationTime();
        computationTime += algorithm.getComputationTime();
      }

      finalMass = communication.reduce(fieldList.density.getLocalData());

      differenceMass = fabs(initialMass-finalMass)/initialMass;
      totalTime = computationTime + communicationTime
        + writeFieldTime + writeAnalysisTime;
      printOutputs();
    }


  protected:
    void printInputs() {
      if (communication.rankMPI == MathVector<int, 3>({0, 0, 0})) {
        std::cout.precision(15);
        std::cout << "-------------------OPTIONS-------------------" << std::endl
                  << "Lattice                  : D" << L::dimD << "Q" << L::dimQ << std::endl
                  << "Global lengths           : " << gSD::sLength() << std::endl
                  << "Global memory            : " << gSD::sVolume()*sizeof(dataT) << "B"
                  << std::endl
                  << "----------------------------------------------" << std::endl
                  << "NPROCS                   : " << numProcs << "" << std::endl
                  << "NTHREADS                 : " << numThreads << "" << std::endl
                  << "-------------------PARAMETERS-----------------" << std::endl
                  << "Relaxation time          : " << relaxationTime << std::endl
                  << "Viscosity                : " << L::cs2 * (relaxationTime - 0.5)
                  << std::endl
                  << "Start iteration          : " << startIteration << std::endl
                  << "End iteration            : " << endIteration << std::endl
                  << "----------------------------------------------" << std::endl;
      }
    }

    void printOutputs() {
      if (communication.rankMPI == MathVector<int, 3>({0, 0, 0})) {
        std::cout << "-------------------OUTPUTS-_------------------" << std::endl
                  << "Total time               : " << totalTime << " s" << std::endl
                  << "Computatation time       : " << computationTime << " s" << std::endl
                  << "Communication time       : " << communicationTime << " s" << std::endl
                  << "Analysis time            : " << writeAnalysisTime << " s" << std::endl
                  << "Write time               : " << writeFieldTime << " s" << std::endl;

        const double mlups = (gSD::sVolume() * 1e-6)/(totalTime / (endIteration-startIteration+1));

        std::cout << "MLUPS                   : " << mlups << std::endl
                  << "Initial mass            : " << initialMass << std::endl
                  << "Final mass              : " << finalMass << std::endl
                  << "% mass diff.            : " << differenceMass << std::endl
                  << "----------------------------------------------" << std::endl;
      }
    }

    void writeFields(const unsigned int iteration) {
      { INSTRUMENT_ON("Routine<T>::writeFields",2) }

      if(fieldWriter.getIsWritten(iteration)) {
        fieldWriter.openFile(iteration);

        fieldList.writeFields();
        fieldWriter.closeFile();
      }

      if(distributionWriter.getIsBackedUp(iteration)) {
        algorithm.pack();
        distributionWriter.openFile(iteration);
        distributionWriter.writeDistribution(distribution);
        distributionWriter.closeFile();
      }

    }

    void writeAnalyses(const unsigned int iteration) {
      { INSTRUMENT_ON("Routine<T>::writeFields",2) }

      if(scalarAnalysisList.getIsAnalyzed(iteration)) {
        scalarAnalysisList.writeAnalyses(iteration);
      }

      if(spectralAnalysisList.getIsAnalyzed(iteration)) {
        spectralAnalysisList.writeAnalyses(iteration);
      }
    }


  };

} // namespace lbm
