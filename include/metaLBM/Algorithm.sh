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

  template<class T, DomainType initDomainType>
  class Algorithm<T, AlgorithmType::Pull, initDomainType, Architecture::NVSHMEM>
    : public Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic> {
  private:
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::localDensity_Ptr;
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::localVelocity_Ptr;
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::localForce_Ptr;
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::localAlpha_Ptr;

    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::haloDistribution_Previous_Ptr;
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::haloDistribution_Next_Ptr;

    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::communication;
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::collision;
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::moment;
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::computationLocal;
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::computationHalo;

    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::dtComputation;
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::dtCommunication;
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::dtTotal;
    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::isWritten;

    using Algorithm<T, AlgorithmType::Generic, initDomainType, Architecture::Generic>::storeLocalFields;

    Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                     partitionningT, L::dimD> periodicBoundary;

  public:
    Algorithm(Field<T, 1, initDomainType, Architecture::GPU, writeDensity>& densityField_in,
              Field<T, L::dimD, initDomainType, Architecture::GPU, writeVelocity>& velocityField_in,
              Field<T, L::dimD, initDomainType, Architecture::GPU, writeDensity>& forceField_in,
              Field<T, 1, initDomainType, Architecture::GPU, writeDensity>& alphaField_in,
              Distribution<T, initDomainType, Architecture::GPU>& f_Previous_in,
              Distribution<T, initDomainType, Architecture::GPU>& f_Next_in,
              Communication<T, latticeT, algorithmT, memoryL,
              partitionningT, L::dimD>& communication_in)
      : Algorithm<T, AlgorithmType::Generic, initDomainType,
                  Architecture::Generic>(densityField_in, velocityField_in,
                                forceField_in, alphaField_in,
                                f_Previous_in, f_Next_in,
                                communication_in)
    {}

    HOST DEVICE
    void operator()(const MathVector<unsigned int, 3>& iP) {
      moment.calculateMoments(haloDistribution_Previous_Ptr, iP);

      collision.setForce(localForce_Ptr, iP, gSD::offset(communication.getRankMPI()));
      collision.setVariables(haloDistribution_Previous_Ptr, iP,
                             moment.getDensity(), moment.getVelocity());

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        haloDistribution_Next_Ptr[hSD::getIndex(iP, iQ)] =
                            collision.calculate(haloDistribution_Previous_Ptr,
                                                iP-uiL::celerity()[iQ], iQ);
      }

      storeLocalFields(iP);
    }

    void setIsWritten(bool isWritten_in) {
      isWritten = isWritten_in;
    }

    HOST
    void iterate(const unsigned int iteration) {
      INSTRUMENT_ON("Algorithm<T, AlgorithmType::Pull>::iterate",2)

      std::swap(haloDistribution_Previous_Ptr, haloDistribution_Next_Ptr);

      collision.update(iteration);

      auto t0 = std::chrono::high_resolution_clock::now();

      communication.communicateHalos(haloDistribution_Previous_Ptr);
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

    using Algorithm<T, AlgorithmType::Generic, initDomainType,
                    Architecture::Generic>::getCommunicationTime;
    using Algorithm<T, AlgorithmType::Generic, initDomainType,
                    Architecture::Generic>::getComputationTime;
    using Algorithm<T, AlgorithmType::Generic, initDomainType,
                    Architecture::Generic>::getTotalTime;
  };

}

#endif // ALGORITHM_H
