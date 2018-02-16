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

  template<class T, Architecture architecture, InputOutputType inputOutputType>
  class Routine {};

  template<class T, Architecture architecture>
  class Routine<T, architecture, InputOutputType::Generic> {
  protected:
    Communication<T, latticeT, algorithmT,
                  memoryL, partitionningT, L::dimD> communication;

    Writer_ writer;

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
      , writer(prefix, rankMPI_in)
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
                  << "Global lengths  : " << gD::length() << std::endl
                  << "Global memory   : " << gD::volume()*sizeof(dataT) << "B" << std::endl
                  << "----------------------------------------------" << std::endl
                  << "Local lengths   : " << lD::length() << std::endl
                  << "Local memory    : " << lD::volume()*sizeof(dataT) << "B" << std::endl
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

        const double mlups = (gD::volume() * 1e-6)/(totalTime / (endIteration-startIteration+1));

        std::cout << "MLUPS           : " << mlups << std::endl
                  << "Initial mass    : " << initialMass << std::endl
                  << "Final mass      : " << finalMass << std::endl
                  << "% mass diff.    : " << differenceMass << std::endl
                  << "----------------------------------------------" << std::endl;
      }
    }

  };


  template<class T, Architecture architecture>
  class Routine<T, architecture, InputOutputType::Serial>
    : public Routine<T, architecture, InputOutputType::Generic> {
  private:
    Field<T, 1, DomainType::Global, architecture, writeDensity> densityField;
    Field<T, L::dimD, DomainType::Global, architecture, writeVelocity> velocityField;
    Field<T, L::dimD, DomainType::Global, architecture, writeDensity> forceField;
    Field<T, 1, DomainType::Global, architecture, writeDensity> alphaField;

    Distribution<T, DomainType::Global, architecture> f_Previous;
    Distribution<T, DomainType::Global, architecture> f_Next;
    Algorithm<dataT, algorithmT, DomainType::Global, arch> algorithm;

    using Routine<T, architecture, InputOutputType::Generic>::communication;
    using Routine<T, architecture, InputOutputType::Generic>:: initialMass;
    using Routine<T, architecture, InputOutputType::Generic>:: finalMass;
    using Routine<T, architecture, InputOutputType::Generic>:: differenceMass;
    using Routine<T, architecture, InputOutputType::Generic>:: computationTime;
    using Routine<T, architecture, InputOutputType::Generic>:: communicationTime;
    using Routine<T, architecture, InputOutputType::Generic>:: writeTime;
    using Routine<T, architecture, InputOutputType::Generic>:: totalTime;

    using Routine<T, architecture, InputOutputType::Generic>::writer;

  public:
    Routine(const MathVector<int, 3>& rankMPI_in,
            const MathVector<int, 3>& sizeMPI_in,
            const std::string& processorName_in)
      :  Routine<T, architecture, InputOutputType::Generic>(rankMPI_in, sizeMPI_in,
                                                            processorName_in)
      , densityField(initGlobalDensity<T>())
      , velocityField(initGlobalVelocity<T>())
      , forceField("force")
      , alphaField(initGlobalAlpha<T>())
      , f_Previous("previousDistribution", initGlobalDistribution<T>(densityField,
                                                                     velocityField))
      , f_Next("nextDistribution", initGlobalDistribution<T>(densityField,
                                                            velocityField))
      , algorithm(densityField, velocityField, forceField, alphaField,
                  f_Previous, f_Next, communication)
    {}

    void compute() {
      INSTRUMENT_ON("Routine<T>::compute",1)

      printInputs();

      initializeLocalFields();

      writeFields(startIteration);

      initialMass = communication.reduce(densityField.localDeviceArray(),
                                         densityField.localHostArray());

      for(int iteration = startIteration+1; iteration <= endIteration; ++iteration) {
        algorithm.setIsWritten(writer.getIsWritten(iteration));
        // what is this set(get(itration))????
        algorithm.iterate();

        writeFields(iteration);

        communicationTime += algorithm.getCommunicationTime();
        computationTime += algorithm.getComputationTime();
        totalTime += algorithm.getTotalTime();
      }


      finalMass = communication.reduce(densityField.localDeviceArray(),
                                       densityField.localHostArray());

      differenceMass = fabs(initialMass-finalMass)/initialMass;
      printOutputs();
    }

  protected:
    using Routine<T, architecture, InputOutputType::Generic>::printInputs;
    using Routine<T, architecture, InputOutputType::Generic>::printOutputs;

    void initializeLocalFields() {
      INSTRUMENT_ON("Routine<T>::initializeLocalFields",2)

      communication.sendGlobalToLocal(densityField.globalArray(),
                                      densityField.localHostArray(),
                                      densityField.localDeviceArray(),
                                      densityField.numberComponents);

      communication.sendGlobalToLocal(velocityField.globalArray(),
                                      velocityField.localHostArray(),
                                      velocityField.localDeviceArray(),
                                      velocityField.numberComponents);

      communication.sendGlobalToLocal(alphaField.globalArray(),
                                      alphaField.localHostArray(),
                                      alphaField.localDeviceArray(),
                                      alphaField.numberComponents);

      communication.sendGlobalToLocal(f_Previous.globalArray(),
                                      f_Previous.localHostArray(),
                                      f_Previous.localDeviceArray(),
                                      f_Previous.numberComponents);
      f_Previous.unpackLocal();
      f_Previous.haloDeviceArray().copyFrom(f_Previous.haloHostArray());

      communication.sendGlobalToLocal(f_Next.globalArray(),
                                      f_Next.localHostArray(),
                                      f_Next.localDeviceArray(),
                                      f_Next.numberComponents);
      f_Next.unpackLocal();
      f_Next.haloDeviceArray().copyFrom(f_Next.haloHostArray());

    }

    void writeFields(const int iteration) {
      INSTRUMENT_ON("Routine<T>::writeFields",2)

      if(writer.getIsWritten(iteration)) {
          communication.sendLocalToGlobal(densityField.localDeviceArray(),
                                          densityField.localHostArray(),
                                          densityField.globalArray(),
                                          densityField.numberComponents);
          communication.sendLocalToGlobal(velocityField.localDeviceArray(),
                                          velocityField.localHostArray(),
                                          velocityField.globalArray(),
                                          velocityField.numberComponents);
          communication.sendLocalToGlobal(alphaField.localDeviceArray(),
                                          alphaField.localHostArray(),
                                          alphaField.globalArray(),
                                          alphaField.numberComponents);
          communication.sendLocalToGlobal(forceField.localDeviceArray(),
                                          forceField.localHostArray(),
                                          forceField.globalArray(),
                                          forceField.numberComponents);

        writer.openFile(iteration);

        writer.writeField(densityField);
        writer.writeField(velocityField);
        writer.writeField(alphaField);
        writer.writeField(forceField);

      if(writer.getIsBackedUp(iteration)) {

        f_Previous.packLocal();
        communication.sendLocalToGlobal(f_Previous.localDeviceArray(),
                                        f_Previous.localHostArray(),
                                        f_Previous.globalArray(),
                                        f_Previous.numberComponents);

        writer.writeField(f_Previous);
      }

      writer.closeFile();
      }
    }

  };


  template<class T, Architecture architecture>
  class Routine<T, architecture, InputOutputType::Parallel>
    : public Routine<T, architecture, InputOutputType::Generic> {
  private:
    Field<T, 1, DomainType::Local, architecture, writeDensity> densityField;
    Field<T, L::dimD, DomainType::Local, architecture, writeVelocity> velocityField;
    Field<T, L::dimD, DomainType::Local, architecture, writeDensity> forceField;
    Field<T, 1, DomainType::Local, architecture, writeDensity> alphaField;

    Distribution<T, DomainType::Local, architecture> f_Previous;
    Distribution<T, DomainType::Local, architecture> f_Next;
    Algorithm<dataT, algorithmT, DomainType::Local, arch> algorithm;

    using Routine<T, architecture, InputOutputType::Generic>::communication;
    using Routine<T, architecture, InputOutputType::Generic>:: initialMass;
    using Routine<T, architecture, InputOutputType::Generic>:: finalMass;
    using Routine<T, architecture, InputOutputType::Generic>:: differenceMass;
    using Routine<T, architecture, InputOutputType::Generic>:: computationTime;
    using Routine<T, architecture, InputOutputType::Generic>:: communicationTime;
    using Routine<T, architecture, InputOutputType::Generic>:: writeTime;
    using Routine<T, architecture, InputOutputType::Generic>:: totalTime;

    using Routine<T, architecture, InputOutputType::Generic>::writer;

  public:
    Routine(const MathVector<int, 3>& rankMPI_in,
            const MathVector<int, 3>& sizeMPI_in,
            const std::string& processorName_in)
      :  Routine<T, architecture, InputOutputType::Generic>(rankMPI_in, sizeMPI_in,
                                                            processorName_in)
      , densityField(initLocalDensity<T>())
      , velocityField(initLocalVelocity<T>())
      , forceField("force")
      , alphaField(initLocalAlpha<T>())
      , f_Previous("previousDistribution", initLocalDistribution<T>(densityField,
                                                                     velocityField))
      , f_Next("nextDistribution", initLocalDistribution<T>(densityField,
                                                             velocityField))
      , algorithm(densityField, velocityField, forceField, alphaField,
                  f_Previous, f_Next, communication)
    {}

    void compute() {
      INSTRUMENT_ON("Routine<T>::compute",1)

      printInputs();

      initializeLocalFields();

      writeFields(startIteration);

      initialMass = communication.reduce(densityField.localDeviceArray(),
                                         densityField.localHostArray());

      for(int iteration = startIteration+1; iteration <= endIteration; ++iteration) {
        algorithm.setIsWritten(writer.getIsWritten(iteration));
        // what is this set(get(itration))????
        algorithm.iterate();

        writeFields(iteration);

        communicationTime += algorithm.getCommunicationTime();
        computationTime += algorithm.getComputationTime();
        totalTime += algorithm.getTotalTime();
      }


      finalMass = communication.reduce(densityField.localDeviceArray(),
                                       densityField.localHostArray());

      differenceMass = fabs(initialMass-finalMass)/initialMass;
      printOutputs();
    }

  protected:
    using Routine<T, architecture, InputOutputType::Generic>::printInputs;
    using Routine<T, architecture, InputOutputType::Generic>::printOutputs;

    void initializeLocalFields() {
      INSTRUMENT_ON("Routine<T>::initializeLocalFields",2)

      f_Previous.unpackLocal();
      f_Previous.haloDeviceArray().copyFrom(f_Previous.haloHostArray());

      f_Next.unpackLocal();
      f_Next.haloDeviceArray().copyFrom(f_Next.haloHostArray());

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

          f_Previous.packLocal();

          writer.writeField(f_Previous);
        }

        writer.closeFile();
      }
    }

  };

}




#endif // ROUTINE_H
