#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <chrono>
#include <omp.h>

#include "Input.h"
#include "Options.h"
#include "Initialize.h"
#include "Domain.h"
#include "Field.h"
#include "Distribution.h"
#include "Moment.h"
#include "Collision.h"
#include "Force.h"
#include "Boundary.h"
#include "Reader.h"
#include "Writer.h"
#include "Communication.h"
#include "Computation.h"

namespace lbm {

  template<class T>
  class Algorithm {
  private:
    const MathVector<unsigned int, 3> rankMPI;
    const MathVector<unsigned int, 3> sizeMPI;
    const std::string processorName;

    Communication_ communication;

    Field<T, 1, writeDensity> densityField;
    Field<T, L::dimD, writeVelocity> velocityField;
    Field<T, L::dimD, writeDensity> forceField;
    Field<T, 1, writeDensity> alphaField;
    // Field<T, 1, writeDensity> entropyField;

    Distribution<T> f_Previous;
    Distribution<T> f_Next;

    Moment<T> moment;
    Collision_ collision;
    Boundary_ boundary;
    Writer_ writer;

    double mass;

    std::chrono::duration<double> dtComputation;
    std::chrono::duration<double> dtCommunication;
    std::chrono::duration<double> dtTotal;

    double computationTime;
    double communicationTime;
    double writeTime;
    double totalTime;

  public:
    Algorithm(const MathVector<unsigned int, 3>& rankMPI_in,
              const MathVector<unsigned int, 3>& sizeMPI_in,
              const std::string& processorName_in)
      : rankMPI(rankMPI_in)
      , sizeMPI(sizeMPI_in)
      , processorName(processorName_in)
      , communication(rankMPI_in, sizeMPI_in, processorName_in)
      , densityField("density")
      , velocityField("velocity")
      , forceField("force")
      , alphaField("alpha")
      , f_Previous("previousDistribution")
      , f_Next("nextDistribution")
      , collision(relaxationTime, forceAmplitude, forceWaveLength)
      , boundary()
      , writer(prefix, rankMPI_in)
      , mass(0.0)
      , computationTime(0.0)
      , communicationTime(0.0)
      , writeTime(0.0)
      , totalTime(0.0)
    {}

    Algorithm(const Algorithm<T>& algorithm_in)
      : rankMPI(algorithm_in.rankMPI)
      , sizeMPI(algorithm_in.sizeMPI)
      , processorName(algorithm_in.processorName)
      , communication(algorithm_in.communication)
      , densityField(algorithm_in.densityField)
      , velocityField(algorithm_in.velocityField)
      , forceField(algorithm_in.forceField)
      , alphaField(algorithm_in.alphaField)
      , f_Previous(algorithm_in.f_Previous)
      , f_Next(algorithm_in.f_Next)
      , collision(algorithm_in.collision)
      , boundary(algorithm_in.boundary)
      , writer(algorithm_in.writer)
      , mass(algorithm_in.mass)
      , computationTime(algorithm_in.computationTime)
      , communicationTime(algorithm_in.communicationTime)
      , writeTime(algorithm_in.writeTime)
      , totalTime(algorithm_in.totalTime)
    {}

    void computeLBM() {
      printInputs();

      initializeFields();

      writeFields(startIteration);

      mass = communication.reduceLocal(densityField.localData());

      for(int iteration = startIteration+1; iteration <= endIteration; ++iteration) {

        iterateLBM(iteration);

        writeFields(iteration);

        communicationTime += dtCommunication.count();
        computationTime += dtComputation.count();
        totalTime += dtTotal.count();
      }

      mass = fabs(mass-communication.reduceLocal(densityField.localData()))/mass;
      printOutputs();
    }

  private:
    void printInputs() {

      communication.printInputs();

      if (rankMPI == MathVector<unsigned int, 3>{{ 0 }}) {
        std::cout.precision(15);
        std::cout << "-------------------OPTIONS-------------------" << std::endl;
        std::cout << "Lattice         : D" << L::dimD << "Q" << L::dimQ << std::endl;
        std::cout << "Global lengths  : " << gD::length() << std::endl;
        std::cout << "Global memory   : " << (int)(gD::volume()
                                                 *sizeof(dataT))/(1<<30) << std::endl;
        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "Rank MPI        : " << rankMPI << std::endl;
        std::cout << "Local lengths   : " << lD::length() << std::endl;
        std::cout << "Local memory    : " << (int)(lD::volume()
                                                 *sizeof(dataT))/(1<<30) << std::endl;
        std::cout << "----------------------------------------------" << std::endl;
        std::cout << "NPROCS          : " << NPROCS << "" << std::endl;
        std::cout << "NTHREADS        : " << NTHREADS << "" << std::endl;
        std::cout << "STRUCTURE       : SOA " << std::endl;
        std::cout << "-------------------PARAMETERS-----------------" << std::endl;
        std::cout << "Relaxation time : " << relaxationTime << std::endl;
        std::cout << "Viscosity       : " << L::cs2 * (relaxationTime - 0.5) << std::endl;
        std::cout << "Start iteration : " << startIteration << std::endl;
        std::cout << "End iteration   : " << endIteration << std::endl;

      }
    }

    void printOutputs() {
      if (rankMPI == MathVector<unsigned int, 3> {{ 0 }}) {
        std::cout << "--------------------OUTPUTS------------------" << std::endl;
        std::cout << "Total time      : " << totalTime << " s" << std::endl;
        std::cout << "Comp time       : " << computationTime << " s" << std::endl;
        std::cout << "Comm time       : " << communicationTime << " s" << std::endl;

        const double mlups = (gD::volume() * 1e-6)/(totalTime / (endIteration-startIteration+1));

        std::cout << "MLUPS           : " << mlups << std::endl;
        std::cout << "% mass diff.    : " << mass << std::endl;
        std::cout << "---------------------------------------------\n" << std::endl;
      }
    }

    void initializeFields() {
      densityField.setGlobalField(initGlobalDensity<T>());
      communication.sendGlobalToLocal(densityField.globalData(),
                                      densityField.localData(),
                                      densityField.numberComponents);

      velocityField.setGlobalField(initGlobalVelocity<T>());
      communication.sendGlobalToLocal(velocityField.globalData(),
                                      velocityField.localData(),
                                      velocityField.numberComponents);

      alphaField.setGlobalField(initGlobalAlpha<T>());
      communication.sendGlobalToLocal(alphaField.globalData(),
                                      alphaField.localData(),
                                      alphaField.numberComponents);

      f_Previous.setGlobalField(initGlobalDistribution<T>(densityField.getGlobalField(),
                                                          velocityField.getGlobalField()));
      communication.sendGlobalToLocal(f_Previous.globalData(),
                                      f_Previous.localData(),
                                      f_Previous.numberComponents);
      f_Previous.unpackLocal();

      f_Next.setGlobalField(initGlobalDistribution<T>(densityField.getGlobalField(),
                                                      velocityField.getGlobalField()));
      communication.sendGlobalToLocal(f_Next.globalData(),
                                      f_Next.localData(),
                                      f_Next.numberComponents);
      f_Next.unpackLocal();

    }

    void storeLocalFields(const MathVector<unsigned int, 3>& iP,
                          const unsigned int iteration) {
      if(iteration%writeStep == 0) {
        unsigned int indexLocal = hD::getIndexLocal(iP);

        densityField.setLocalField(indexLocal, moment.getDensity());
        velocityField.setLocalField(indexLocal, collision.getHydrodynamicVelocity());
        alphaField.setLocalField(indexLocal, collision.getAlpha());
        forceField.setLocalField(indexLocal, collision.getForce());
      }
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
                          densityField.globalData());
        writer.writeField(velocityField.fieldName, velocityField.numberComponents,
                          velocityField.globalData());
        writer.writeField(alphaField.fieldName, alphaField.numberComponents,
                          alphaField.globalData());
        writer.writeField(forceField.fieldName, forceField.numberComponents,
                          forceField.globalData());

        if(iteration%backupStep == 0) {
          f_Previous.packLocal();

          if(writer.isSerial) {
            communication.sendLocalToGlobal(f_Previous.localData(),
                                            f_Previous.globalData(),
                                            f_Previous.numberComponents);
          }

          writer.writeField(f_Previous.fieldName, f_Previous.numberComponents,
                            f_Previous.globalData());
        }

        writer.closeFile();
      }
    }

    void iterateLBM(const unsigned int iteration) {
      f_Previous.swapHalo(f_Next);

      //force.update(iteration);

      auto t0 = std::chrono::high_resolution_clock::now();
      communication.periodic(f_Previous.haloData());

      //boundary.apply(f_Previous.haloData());

      auto t1 = std::chrono::high_resolution_clock::now();


      Computation::Do([&] (MathVector<unsigned int, 3>& iP) {
          moment.computeMoments(f_Previous.haloData(), iP);

          collision.setForce(iP+gD::offset(rankMPI));
          collision.setVariables(f_Previous.haloData(), iP,
                                 moment.getDensity(), moment.getVelocity());

          UnrolledFor<0, L::dimQ>::Do([&] (unsigned int iQ) {
              f_Next.setHaloField(hD::getIndex(iP, iQ),
                                  collision.postDistribution(f_Previous.haloData(),
                                                             iP-projectionI(L::celerity()[iQ]), iQ));
          });

          storeLocalFields(iP, iteration);
        });

      auto t2 = std::chrono::high_resolution_clock::now();

      dtCommunication = t1 - t0;
      dtComputation = t2 - t1;
      dtTotal = t1 - t0;
    }

  };


  typedef Algorithm<dataT> Algorithm_;

}




#endif // ALGORITHM_H
