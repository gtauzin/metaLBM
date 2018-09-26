#pragma once

#include <chrono>

#include "Algorithm.h"
#include "AnalysisList.h"
#include "Commons.h"
#include "Communication.h"
#include "Computation.h"
#include "Event.h"
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
  template <class T, AlgorithmType algorithmType, Architecture architecture,
            MemoryLayout memoryLayout, PartitionningType partitionningType,
            CommunicationType communicationType, Overlapping overlapping>
  class Routine {
  protected:
    using Clock = std::chrono::high_resolution_clock;
    using Seconds = std::chrono::duration<double>;

    using Algorithm_
    = Algorithm<T, algorithmType, architecture, memoryLayout, partitionningType,
                communicationType, overlapping>;

    Communication_ communication;
    Stream<architecture> defaultStream;
    Stream<architecture> bulkStream;
    Stream<architecture> leftStream;
    Stream<architecture> rightStream;

    Event<architecture> leftEvent;
    Event<architecture> rightEvent;

    FieldWriter_ fieldWriter;
    DistributionWriter_ distributionWriter;
    FieldList<T, architecture> fieldList;
    Distribution<T, architecture> distribution;

    Curl<double, Architecture::CPU, PartitionningType::OneD, L::dimD, L::dimD>
    curlVelocity;
    ScalarAnalysisList<T, architecture> scalarAnalysisList;
    SpectralAnalysisList<T, architecture> spectralAnalysisList;

    Algorithm_ algorithm;
    PerformanceAnalysisList performanceAnalysisList;

  public:
    Routine()
      : communication()
      , defaultStream(true)
      , bulkStream(false)
      , leftStream(false)
      , rightStream(false)
      , leftEvent()
      , rightEvent()
      , fieldWriter(prefix)
      , distributionWriter(prefix)
      , fieldList(fieldWriter, defaultStream)
      , curlVelocity(fieldList.velocity.getData(FFTWInit::numberElements),
                     fieldList.vorticity.getData(FFTWInit::numberElements),
                     Cast<unsigned int, ptrdiff_t, 3>::Do(gSD::sLength()).data(),
                     gFD::offset(MPIInit::rank))
      , distribution(initDistribution<T, architecture>(fieldList.density,
                                                       fieldList.velocity,
                                                       defaultStream))
      , scalarAnalysisList(fieldList, communication, scalarAnalysisStep,
                           startIteration)
      , spectralAnalysisList(fieldList, communication, spectralAnalysisStep,
                             startIteration)
      , algorithm(fieldList, distribution, communication)
      , performanceAnalysisList(performanceAnalysisStep, startIteration)
    {
      printInputs();
    }

    ~Routine() {
      printOutputs();
    }

    void compute() {
      LBM_INSTRUMENT_ON("Routine<T>::compute", 1)

        algorithm.unpack(defaultStream);

      Clock::time_point t0;
      Clock::time_point t1;

      if (writeFieldInit || writeAnalysisInit) {
        curlVelocity.executeSpace();
        curlVelocity.normalize();
      }

      if (writeFieldInit) {
        t0 = Clock::now();
        writeFields(startIteration);
        t1 = Clock::now();
        performanceAnalysisList.updateWriteFieldTime(Seconds(t1 - t0).count());
      }

      if (writeAnalysisInit) {
        t0 = Clock::now();
        writeAnalyses(startIteration);
        t1 = Clock::now();
        performanceAnalysisList.updateWriteAnalysisTime(Seconds(t1 - t0).count());
      }

      performanceAnalysisList.setInitialMass(
                                             communication.reduce(fieldList.density.getData(FFTWInit::numberElements)));

      // Execute LBM algorithm
      for (int iteration = startIteration + 1; iteration <= endIteration; ++iteration) {
        algorithm.isStored = (fieldWriter.getIsWritten(iteration)
                              || scalarAnalysisList.getIsAnalyzed(iteration)
                              || spectralAnalysisList.getIsAnalyzed(iteration));

        algorithm.iterate(iteration, defaultStream, bulkStream, leftStream, rightStream,
                          leftEvent, rightEvent);

        if (algorithm.isStored) {
          curlVelocity.executeSpace();
          curlVelocity.normalize();
        }

        t0 = Clock::now();
        writeFields(iteration);
        t1 = Clock::now();
        performanceAnalysisList.updateWriteFieldTime(Seconds(t1 - t0).count());

        t0 = Clock::now();
        writeAnalyses(iteration);
        t1 = Clock::now();
        performanceAnalysisList.updateWriteAnalysisTime(Seconds(t1 - t0).count());

        performanceAnalysisList.updateCommunicationTime(algorithm.getCommunicationTime());
        performanceAnalysisList.updateComputationTime(algorithm.getComputationTime());

      }

      performanceAnalysisList.updateMass(communication.reduce(
                                                              fieldList.density.getData(FFTWInit::numberElements)));

      performanceAnalysisList.updateMLUPS(endIteration - startIteration);

    }

  protected:
    void printInputs() {
      if (MPIInit::rank[d::X] == 0) {
        std::cout.precision(15);
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
      }
    }

    void printOutputs() {
      if (MPIInit::rank[d::X] == 0) {
        std::cout << "-------------------OUTPUTS--------------------\n"
                  << "Total time               : "
                  <<  performanceAnalysisList.getTotalTime() << " s\n"
                  << "Computatation time       : "
                  << performanceAnalysisList.getComputationTime() << " s\n"
                  << "Communication time       : "
                  << performanceAnalysisList.getCommunicationTime() << " s\n"
                  << "Analysis time            : "
                  << performanceAnalysisList.getWriteAnalysisTime() << " s\n"
                  << "Write time               : "
                  << performanceAnalysisList.getWriteFieldTime() << " s\n";

        std::cout << "MLUPS                   : "
                  << performanceAnalysisList.getMLUPS() << "\n"
                  << "Initial mass            : "
                  << performanceAnalysisList.getInitialMass() << "\n"
                  << "Final mass              : "
                  << performanceAnalysisList.getFinalMass() << "\n"
                  << "% mass diff.            : "
                  << performanceAnalysisList.getDifferenceMass() << "\n"
                  << "----------------------------------------------\n";
      }
    }

    void writeFields(const unsigned int iteration) {
      LBM_INSTRUMENT_ON("Routine<T>::writeFields", 2)

        if (fieldWriter.getIsWritten(iteration)) {
          fieldWriter.openFile(iteration);

          fieldList.writeFields();
          fieldWriter.closeFile();
        }

      if (distributionWriter.getIsBackedUp(iteration)) {
        algorithm.pack(defaultStream);
        distributionWriter.openFile(iteration);
        distributionWriter.writeDistribution(distribution);
        distributionWriter.closeFile();
      }
    }

    void writeAnalyses(const unsigned int iteration) {
      LBM_INSTRUMENT_ON("Routine<T>::writeFields", 2)

        if (scalarAnalysisList.getIsAnalyzed(iteration)) {
          scalarAnalysisList.writeAnalyses(iteration);
        }

      if (spectralAnalysisList.getIsAnalyzed(iteration)) {
        spectralAnalysisList.writeAnalyses(iteration);
      }

      if (performanceAnalysisList.getIsAnalyzed(iteration)) {
        performanceAnalysisList.updateMass(communication.reduce(
                                                                fieldList.density.getData(FFTWInit::numberElements)));

        performanceAnalysisList.updateMLUPS(iteration - startIteration);

        performanceAnalysisList.writeAnalyses(iteration);
      }
    }

  };

}  // namespace lbm
