#ifndef ROUTINE_H
#define ROUTINE_H

#include <chrono>

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "Initialize.h"
#include "Domain.h"
#include "Field.h"
#include "Distribution.h"
#include "Algorithm.h"
#include "Writer.h"
#include "Communication.h"
#include "Computation.h"

namespace lbm {

  template<class T, Architecture architecture, Implementation implementation,
           InputOutputType inputOutputType>
  class Routine {};

  template<class T, Architecture architecture, Implementation implementation>
  class Routine<T, architecture, implementation, InputOutputType::Generic> {
  protected:
    Communication<T, latticeT, algorithmT,
      memoryL, partitionningT, implementation, L::dimD> communication;

    Writer_ writer;
    Writer<T, InputOutput::DAT, InputOutputType::Serial,
      InputOutputDataFormat::ascii> analysisWriter;

    double initialMass;
    double finalMass;
    double differenceMass;
    double computationTime;
    double communicationTime;
    double writeTime;
    double totalTime;

  public:
    Routine(const MathVector<int, 3>& rankMPI_in,
            const MathVector<int, 3>& sizeMPI_in,
            const std::string& processorName_in)
      : communication(rankMPI_in, sizeMPI_in, processorName_in)
      , writer(prefix, rankMPI_in, sizeMPI_in)
      , analysisWriter(prefix, rankMPI_in, sizeMPI_in)
      , initialMass(0.0)
      , finalMass(0.0)
      , differenceMass(0.0)
      , computationTime(0.0)
      , communicationTime(0.0)
      , writeTime(0.0)
      , totalTime(0.0)
    {}


  protected:
    void printInputs() {

      communication.printInputs();

      if (communication.getRankMPI() == MathVector<int, 3>({0, 0, 0})) {
        std::cout.precision(15);
        std::cout << "-------------------OPTIONS-------------------" << std::endl
                  << "Lattice         : D" << L::dimD << "Q" << L::dimQ << std::endl
                  << "Global lengths  : " << gSD::sLength() << std::endl
                  << "Global memory   : " << gSD::sVolume()*sizeof(dataT) << "B" << std::endl
                  << "----------------------------------------------" << std::endl
                  << "Local lengths   : " << lSD::pLength() << std::endl
                  << "Local memory    : " << lSD::pVolume()*sizeof(dataT) << "B" << std::endl
                  << "----------------------------------------------" << std::endl
                  << "NPROCS          : " << NPROCS << "" << std::endl
                  << "NTHREADS        : " << NTHREADS << "" << std::endl
                  << "DATA STRUCTURE  : SoA " << std::endl
                  << "-------------------PARAMETERS-----------------" << std::endl
                  << "Relaxation time : " << relaxationTime << std::endl
                  << "Viscosity       : " << L::cs2 * (relaxationTime - 0.5) << std::endl
                  << "Start iteration : " << startIteration << std::endl
                  << "End iteration   : " << endIteration << std::endl
                  << "----------------------------------------------" << std::endl;
      }
    }

    void printOutputs() {
      if (communication.getRankMPI() == MathVector<int, 3>({0, 0, 0})) {
        std::cout << "--------------------OUTPUTS-------------------" << std::endl
                  << "Total time      : " << totalTime << " s" << std::endl
                  << "Comp time       : " << computationTime << " s" << std::endl
                  << "Comm time       : " << communicationTime << " s" << std::endl;

        const double mlups = (gSD::sVolume() * 1e-6)/(totalTime / (endIteration-startIteration+1));

        std::cout << "MLUPS           : " << mlups << std::endl
                  << "Initial mass    : " << initialMass << std::endl
                  << "Final mass      : " << finalMass << std::endl
                  << "% mass diff.    : " << differenceMass << std::endl
                  << "----------------------------------------------" << std::endl;
      }
    }

  };


  template<class T, Architecture architecture, Implementation implementation>
  class Routine<T, architecture, implementation, InputOutputType::Parallel>
    : public Routine<T, architecture, implementation, InputOutputType::Generic> {
  private:
    using Base = Routine<T, architecture, implementation, InputOutputType::Generic>;

    Field<T, 1, DomainType::LocalSpace, architecture, writeDensity> densityField;
    Field<T, L::dimD, DomainType::LocalSpace, architecture, writeVelocity> velocityField;
    Field<T, L::dimD, DomainType::LocalSpace, architecture, writeDensity> forceField;
    Field<T, 1, DomainType::LocalSpace, architecture, writeDensity> alphaField;

    Distribution<T, DomainType::LocalSpace, architecture> distribution;

    Algorithm<dataT, algorithmT, DomainType::LocalSpace, architecture,
              implementation> algorithm;

    using Base::communication;
    using Base:: initialMass;
    using Base:: finalMass;
    using Base:: differenceMass;
    using Base:: computationTime;
    using Base:: communicationTime;
    using Base:: writeTime;
    using Base:: totalTime;

    using Base::writer;

  public:
    Routine(const MathVector<int, 3>& rankMPI_in,
            const MathVector<int, 3>& sizeMPI_in,
            const std::string& processorName_in,
            const unsigned int numberElements_in)
      :  Routine<T, architecture, implementation,
                 InputOutputType::Generic>(rankMPI_in, sizeMPI_in,
                                           processorName_in)
      , densityField(initLocalDensity<T, architecture>(numberElements_in, rankMPI_in))
      , velocityField(initLocalVelocity<T, architecture>(numberElements_in))
      , forceField(initLocalForce<T, architecture>(numberElements_in, rankMPI_in))
      , alphaField(initLocalAlpha<T, architecture>(numberElements_in))
      , distribution("distribution",
                     initLocalDistribution<T, architecture>(densityField,
                                                            velocityField))
      , algorithm(densityField, velocityField, forceField, alphaField,
                  distribution, communication)
    {}

    void compute() {
      INSTRUMENT_ON("Routine<T>::compute",1)

      printInputs();

      initializeLocalDistribution();

      writeFields(startIteration);

      initialMass = communication.reduce(densityField.getLocalData());

      for(int iteration = startIteration+1; iteration <= endIteration; ++iteration) {
        algorithm.setIsWritten(writer.getIsWritten(iteration));
        algorithm.iterate(iteration);

        //writeAnalyses(iteration);
        writeFields(iteration);

        communicationTime += algorithm.getCommunicationTime();
        computationTime += algorithm.getComputationTime();
        totalTime += algorithm.getTotalTime();
      }

      finalMass = communication.reduce(densityField.getLocalData());

      differenceMass = fabs(initialMass-finalMass)/initialMass;
      printOutputs();
    }

  protected:
    using Base::printInputs;
    using Base::printOutputs;

    void initializeLocalDistribution() {
      INSTRUMENT_ON("Routine<T>::initializeLocalDistribution",2)
      algorithm.unpack();
    }

    void writeFields(const int iteration) {
      INSTRUMENT_ON("Routine<T>::writeFields",2)

      if(writer.getIsWritten(iteration)) {
        writer.openFile(iteration);

        writer.writeField(densityField);
        writer.writeField(velocityField);
        writer.writeField(alphaField);
        writer.writeField(forceField);

      if(writer.getIsBackedUp(iteration)) {
        algorithm.pack();
        writer.writeField(distribution);
      }

      writer.closeFile();
      }
    }

    //using Base::writeAnalyses;


  };

} // namespace lbm




#endif // ROUTINE_H
