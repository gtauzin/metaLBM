#ifndef ROUTINE_H
#define ROUTINE_H

#include <chrono>
#include <omp.h>

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

  template<class T>
  class Routine {
  private:
    Communication_ communication;

    Field<T, 1, writeDensity> densityField;
    Field<T, L::dimD, writeVelocity> velocityField;
    Field<T, L::dimD, writeDensity> forceField;
    Field<T, 1, writeDensity> alphaField;

    Distribution<T> f_Previous;
    Distribution<T> f_Next;
    Algorithm_ algorithm;

    Writer_ writer;

    double mass;
    double computationTime;
    double communicationTime;
    double writeTime;
    double totalTime;

  public:
    Routine(const MathVector<int, 3>& rankMPI_in,
            const MathVector<int, 3>& sizeMPI_in,
            const std::string& processorName_in)
      : communication(rankMPI_in, sizeMPI_in, processorName_in)
      , densityField("density", initGlobalDensity<T>())
      , velocityField("velocity", initGlobalVelocity<T>())
      , forceField("force")
      , alphaField("alpha", initGlobalAlpha<T>())
      , f_Previous("previousDistribution", initGlobalDistribution<T>(densityField,
                                                                     velocityField))
      , f_Next("nextDistribution", initGlobalDistribution<T>(densityField,
                                                             velocityField))
      , algorithm(communication, densityField, velocityField, forceField, alphaField,
                  f_Previous, f_Next)
      , writer(prefix, rankMPI_in)
      , mass(0.0)
      , computationTime(0.0)
      , communicationTime(0.0)
      , writeTime(0.0)
      , totalTime(0.0)
    {}

    void compute() {
      SCOREP_INSTRUMENT("Routine<T>::compute")

      printInputs();

      initializeLocalFields();

      writeFields(startIteration);

      mass = communication.reduce(densityField.localData());

      for(int iteration = startIteration+1; iteration <= endIteration; ++iteration) {

        algorithm.iterate(iteration);

        writeFields(iteration);

        communicationTime += algorithm.getCommunicationTime();
        computationTime += algorithm.getComputationTime();
        totalTime += algorithm.getTotalTime();
      }

      mass = fabs(mass-communication.reduce(densityField.localData()))/mass;
      printOutputs();
    }

  private:
    void printInputs() {

      communication.printInputs();

      if (communication.getRankMPI() == MathVector<int, 3>{0, 0, 0}) {
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
      if (communication.getRankMPI() == MathVector<int, 3>{0, 0, 0}) {
        std::cout << "--------------------OUTPUTS-------------------" << std::endl
                  << "Total time      : " << totalTime << " s" << std::endl
                  << "Comp time       : " << computationTime << " s" << std::endl
                  << "Comm time       : " << communicationTime << " s" << std::endl;

        const double mlups = (gD::volume() * 1e-6)/(totalTime / (endIteration-startIteration+1));

        std::cout << "MLUPS           : " << mlups << std::endl
                  << "% mass diff.    : " << mass << std::endl
                  << "----------------------------------------------" << std::endl;
      }
    }

    void initializeLocalFields() {
      SCOREP_INSTRUMENT("Routine<T>::initializeLocalFields")

      communication.sendGlobalToLocal(densityField.globalData(),
                                      densityField.localData(),
                                      densityField.numberComponents);

      communication.sendGlobalToLocal(velocityField.globalData(),
                                      velocityField.localData(),
                                      velocityField.numberComponents);

      communication.sendGlobalToLocal(alphaField.globalData(),
                                      alphaField.localData(),
                                      alphaField.numberComponents);

      communication.sendGlobalToLocal(f_Previous.globalData(),
                                      f_Previous.localData(),
                                      f_Previous.numberComponents);
      f_Previous.unpackLocal();

      communication.sendGlobalToLocal(f_Next.globalData(),
                                      f_Next.localData(),
                                      f_Next.numberComponents);
      f_Next.unpackLocal();

    }

    void writeFields(const int iteration) {
      if(iteration%writeStep == 0) {
        writer.openFile(iteration);

        if(writer.isSerial) {
          communication.sendLocalToGlobal(densityField.localData(),
                                          densityField.globalData(),
                                          densityField.numberComponents);
          communication.sendLocalToGlobal(velocityField.localData(),
                                          velocityField.globalData(),
                                          velocityField.numberComponents);
          communication.sendLocalToGlobal(alphaField.localData(),
                                          alphaField.globalData(),
                                          alphaField.numberComponents);
          communication.sendLocalToGlobal(forceField.localData(),
                                          forceField.globalData(),
                                          forceField.numberComponents);
        }

        writer.writeField(densityField.fieldName, densityField.numberComponents,
                          densityField.globalData(), gD());
        writer.writeField(velocityField.fieldName, velocityField.numberComponents,
                          velocityField.globalData(), gDD());
        writer.writeField(alphaField.fieldName, alphaField.numberComponents,
                          alphaField.globalData(), gD());
        writer.writeField(forceField.fieldName, forceField.numberComponents,
                          forceField.globalData(), gDD());

        if(iteration%backupStep == 0) {
          f_Previous.packLocal();

          if(writer.isSerial) {
            communication.sendLocalToGlobal(f_Previous.localData(),
                                            f_Previous.globalData(),
                                            f_Previous.numberComponents);
          }

          writer.writeField(f_Previous.fieldName, f_Previous.numberComponents,
                            f_Previous.globalData(), gQD());
        }

        writer.closeFile();
      }
    }

  };


  typedef Routine<dataT> Routine_;

}




#endif // ROUTINE_H
