#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <chrono>
#include <utility>
#include <cstdio>

#include "Commons.h"
#include "Options.h"
#include "Domain.h"
#include "Lattice.h"
#include "Field.h"
#include "Distribution.h"
#include "Moment.h"
#include "Collision.h"
#include "Boundary.h"
#include "Communication.h"
#include "Computation.h"


namespace lbm {

  template<class T, AlgorithmType algorithmType, Architecture architecture>
  class Algorithm {
  public:
    HOST DEVICE
      void operator()(const MathVector<unsigned int, 3> iP);
  };


  template<class T, Architecture architecture>
    class Algorithm<T, AlgorithmType::Generic, architecture> {
  protected:

    T * RESTRICT localDensity_Ptr;
    T * RESTRICT localVelocity_Ptr;
    T * RESTRICT localForce_Ptr;
    T * RESTRICT localAlpha_Ptr;

    T * RESTRICT haloDistribution_Previous_Ptr;
    T * RESTRICT haloDistribution_Next_Ptr;

    Communication<T, latticeT, algorithmT, memoryL,
                  partitionningT, L::dimD> communication;
    Collision_ collision;
    Moment<T> moment;

    Computation<architecture, L::dimD> computationLocal;
    Computation<architecture, L::dimD> computationHalo;

    std::chrono::duration<double> dtComputation;
    std::chrono::duration<double> dtCommunication;
    std::chrono::duration<double> dtTotal;

    bool isWritten;

    Algorithm(Field<T, 1, architecture, writeDensity>& densityField_in,
              Field<T, L::dimD, architecture, writeVelocity>& velocityField_in,
              Field<T, L::dimD, architecture, writeDensity>& forceField_in,
              Field<T, 1, architecture, writeDensity>& alphaField_in,
              Distribution<T, architecture>& f_Previous_in,
              Distribution<T, architecture>& f_Next_in,
              Communication<T, latticeT, algorithmT, memoryL,
              partitionningT, L::dimD>& communication_in)
      : localDensity_Ptr(densityField_in.localComputedData())
      , localVelocity_Ptr(velocityField_in.localComputedData())
      , localForce_Ptr(forceField_in.localComputedData())
      , localAlpha_Ptr(alphaField_in.localComputedData())
      , haloDistribution_Previous_Ptr(f_Previous_in.haloComputedData())
      , haloDistribution_Next_Ptr(f_Next_in.haloComputedData())
      , communication(communication_in)
      , collision(relaxationTime, forceAmplitude, forceWaveLength)
      , moment()
      , computationLocal(lD::start()+L::halo(),
                    lD::end()+L::halo())
      , computationHalo(hD::start(),
                        hD::end())
      , dtComputation()
      , dtCommunication()
      , dtTotal()
      , isWritten()
    {}

    DEVICE HOST
    void storeLocalFields(const MathVector<unsigned int, 3>& iP) {
      INSTRUMENT_OFF("Algorithm<T, AlgorithmType::Pull>::storeLocalFields",4)

      const unsigned int indexLocal = hD::getIndexLocal(iP);

      localDensity_Ptr[indexLocal] = moment.getDensity();
      localAlpha_Ptr[indexLocal] = collision.getAlpha();

      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        localVelocity_Ptr[lDD::getIndex(indexLocal, iD)] = collision.getHydrodynamicVelocity()[iD];
        localForce_Ptr[lDD::getIndex(indexLocal, iD)] = collision.getForce()[iD];
      }

    }

  public:
    double getCommunicationTime() {
      return dtCommunication.count();
    }

    double getComputationTime() {
      return dtComputation.count();
    }

    double getTotalTime() {
      return dtTotal.count();
    }
  };

  template<class T, Architecture architecture>
    class Algorithm<T, AlgorithmType::Pull, architecture>
    : public Algorithm<T, AlgorithmType::Generic, architecture> {
  private:
    using Algorithm<T, AlgorithmType::Generic, architecture>::localDensity_Ptr;
    using Algorithm<T, AlgorithmType::Generic, architecture>::localVelocity_Ptr;
    using Algorithm<T, AlgorithmType::Generic, architecture>::localForce_Ptr;
    using Algorithm<T, AlgorithmType::Generic, architecture>::localAlpha_Ptr;

    using Algorithm<T, AlgorithmType::Generic, architecture>::haloDistribution_Previous_Ptr;
    using Algorithm<T, AlgorithmType::Generic, architecture>::haloDistribution_Next_Ptr;

    using Algorithm<T, AlgorithmType::Generic, architecture>::communication;
    using Algorithm<T, AlgorithmType::Generic, architecture>::collision;
    using Algorithm<T, AlgorithmType::Generic, architecture>::moment;
    using Algorithm<T, AlgorithmType::Generic, architecture>::computationLocal;
    using Algorithm<T, AlgorithmType::Generic, architecture>::computationHalo;

    using Algorithm<T, AlgorithmType::Generic, architecture>::dtComputation;
    using Algorithm<T, AlgorithmType::Generic, architecture>::dtCommunication;
    using Algorithm<T, AlgorithmType::Generic, architecture>::dtTotal;

    using Algorithm<T, AlgorithmType::Generic, architecture>::isWritten;

    using Algorithm<T, AlgorithmType::Generic, architecture>::storeLocalFields;

    Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                     partitionningT, L::dimD> periodicBoundary;

  public:
    Algorithm(Field<T, 1, architecture, writeDensity>& densityField_in,
              Field<T, L::dimD, architecture, writeVelocity>& velocityField_in,
              Field<T, L::dimD, architecture, writeDensity>& forceField_in,
              Field<T, 1, architecture, writeDensity>& alphaField_in,
              Distribution<T, architecture>& f_Previous_in,
              Distribution<T, architecture>& f_Next_in,
              Communication<T, latticeT, algorithmT, memoryL,
              partitionningT, L::dimD>& communication_in)
      : Algorithm<T, AlgorithmType::Generic, architecture>(densityField_in, velocityField_in,
                                                           forceField_in, alphaField_in,
                                                           f_Previous_in, f_Next_in,
                                                           communication_in)
    {}

    HOST DEVICE
    void operator()(const MathVector<unsigned int, 3>& iP) {
      moment.calculateMoments(haloDistribution_Previous_Ptr, iP);

      collision.setForce(iP+gD::offset(communication.getRankMPI()));
      collision.setVariables(haloDistribution_Previous_Ptr, iP,
                             moment.getDensity(), moment.getVelocity());

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        haloDistribution_Next_Ptr[hD::getIndex(iP, iQ)] =
                            collision.calculate(haloDistribution_Previous_Ptr,
                                                iP-uiL::celerity()[iQ], iQ);
      }

      storeLocalFields(iP);
    }

    void setIsWritten(bool isWritten_in) {
      isWritten = isWritten_in;
    }

    HOST
    void iterate() {
      INSTRUMENT_ON("Algorithm<T, AlgorithmType::Pull>::iterate",2)

      std::swap(haloDistribution_Previous_Ptr, haloDistribution_Next_Ptr);

      //force.update(iteration);

      auto t0 = std::chrono::high_resolution_clock::now();

      communication.communicateHalos(haloDistribution_Previous_Ptr);
      // computationHalo.Do((&PeriodicBoundary::apply),
      //                   haloDistribution_Previous_Ptr);
      computationHalo.Do(periodicBoundary,
                         haloDistribution_Previous_Ptr);

      //boundary.apply(f_Previous.haloData());

      auto t1 = std::chrono::high_resolution_clock::now();

      computationLocal.Do(*this);

      auto t2 = std::chrono::high_resolution_clock::now();

      dtCommunication = (t1 - t0);
      dtComputation = (t2 - t1);
      dtTotal = (t2 - t0);
    }

    using Algorithm<T, AlgorithmType::Generic, architecture>::getCommunicationTime;
    using Algorithm<T, AlgorithmType::Generic, architecture>::getComputationTime;
    using Algorithm<T, AlgorithmType::Generic, architecture>::getTotalTime;
  };

}




#endif // ALGORITHM_H
