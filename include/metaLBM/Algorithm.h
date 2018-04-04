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

  template<class T, AlgorithmType algorithmType,
    Architecture architecture, Implementation implementation>
  class Algorithm {
  public:
    HOST DEVICE
      void operator()(const Position iP);
  };


  template<class T, Architecture architecture, Implementation implementation>
  class Algorithm<T, AlgorithmType::Generic, architecture, implementation> {
  protected:
    T * localDensity_Ptr;
    T * localVelocity_Ptr;
    T * localForce_Ptr;
    T * localAlpha_Ptr;

    T * localDistribution_Ptr;
    T * haloDistributionPrevious_Ptr;
    T * haloDistributionNext_Ptr;
    const unsigned int numberElements;

    Packer<T> packer;
    Unpacker<T> unpacker;
    Communication<T, latticeT, algorithmT, memoryL,
                  partitionningT, implementation, L::dimD> communication;
    Collision_ collision;

    Computation<architecture, L::dimD> computationLocal;
    Computation<architecture, L::dimD> computationBottom;
    Computation<architecture, L::dimD> computationTop;
    Computation<architecture, L::dimD> computationFront;
    Computation<architecture, L::dimD> computationBack;

    std::chrono::duration<double> dtComputation;
    std::chrono::duration<double> dtCommunication;

  public:
    bool isStored;

    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in,
              const unsigned int numberElements_in,
              Communication<T, latticeT, algorithmT, memoryL,
              partitionningT, implementation, L::dimD>& communication_in)
      : localDensity_Ptr(fieldList_in.density.getLocalData())
      , localVelocity_Ptr(fieldList_in.velocity.getLocalData())
      , localForce_Ptr(fieldList_in.force.getLocalData())
      , localAlpha_Ptr(fieldList_in.alpha.getLocalData())
      , localDistribution_Ptr(distribution_in.getLocalData())
      , haloDistributionPrevious_Ptr(distribution_in.getHaloDataPrevious())
      , haloDistributionNext_Ptr(distribution_in.getHaloDataNext())
      , numberElements(numberElements_in)
      , communication(communication_in)
      , collision(relaxationTime, forceAmplitude, forceWaveLength, forcekMin, forcekMax)
      , computationLocal(lSD::sStart()+L::halo(),
                         lSD::sEnd()+L::halo())
      , computationBottom({hSD::start()[d::X], L::halo()[d::Y], hSD::start()[d::Z]},
                          {hSD::end()[d::X], 2*L::halo()[d::Y], hSD::end()[d::Z]})
      , computationTop({hSD::start()[d::X], L::halo()[d::Y]+lSD::sLength()[d::Y] - 1,
                        hSD::start()[d::Z]},
                       {hSD::end()[d::X],  2*L::halo()[d::Y]+ lSD::sLength()[d::Y] - 1,
                        hSD::end()[d::Z]})
      , computationFront({hSD::start()[d::X], hSD::start()[d::Y], L::halo()[d::Z]},
                         {hSD::end()[d::X], hSD::end()[d::Y], 2*L::halo()[d::Z]})
      , computationBack({hSD::start()[d::X], hSD::start()[d::Y],
                         L::halo()[d::Z]+lSD::sLength()[d::Z] - 1},
                        {hSD::end()[d::X], hSD::end()[d::Y],
                         2*L::halo()[d::Z] + lSD::sLength()[d::Z] - 1})
      , dtComputation()
      , dtCommunication()
      , isStored(false)
    {}

    HOST
    void pack() {
      computationLocal.Do(packer, localDistribution_Ptr,
                          haloDistributionNext_Ptr, numberElements);
    }

    HOST
    void unpack() {
      computationLocal.Do(unpacker, haloDistributionNext_Ptr,
                          localDistribution_Ptr, numberElements);
    }

    double getCommunicationTime() {
      return dtCommunication.count();
    }

    double getComputationTime() {
      return dtComputation.count();
    }

  protected:
    DEVICE HOST
    void storeLocalFields(const Position& iP) {
      INSTRUMENT_OFF("Algorithm<T, AlgorithmType::Pull>::storeLocalFields",4)

      const auto indexLocal = hSD::getIndexLocal(iP);

      localDensity_Ptr[indexLocal] = collision.getDensity();
      for(auto iD = 0; iD < L::dimD; ++iD) {
        (localVelocity_Ptr+iD*numberElements)[indexLocal]
          = collision.getHydrodynamicVelocity()[iD];
      }

      if(writeAlpha) localAlpha_Ptr[indexLocal] = collision.getAlpha();

      if(writeForce) {
        for(auto iD = 0; iD < L::dimD; ++iD) {
          (localForce_Ptr+iD*numberElements)[indexLocal] = collision.getForce()[iD];
        }
      }
    }

  };


  template<class T, Architecture architecture, Implementation implementation>
  class Algorithm<T, AlgorithmType::Pull, architecture, implementation>
    : public Algorithm<T, AlgorithmType::Generic, architecture, implementation> {
  private:
    using Base = Algorithm<T, AlgorithmType::Generic, architecture, implementation>;

    using Base::localDensity_Ptr;
    using Base::localVelocity_Ptr;
    using Base::localForce_Ptr;
    using Base::localAlpha_Ptr;

    using Base::localDistribution_Ptr;
    using Base::haloDistributionPrevious_Ptr;
    using Base::haloDistributionNext_Ptr;
    using Base::numberElements;

    using Base::communication;
    using Base::collision;

    Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
             partitionningT, implementation, L::dimD> periodicBoundary;

  public:
    using Base::Algorithm;

    HOST DEVICE
    void operator()(const Position& iP) {
      collision.calculateMoments(haloDistributionPrevious_Ptr, iP);

      collision.setForce(localForce_Ptr, iP,
                         numberElements, gSD::sOffset(communication.rankMPI));
      collision.calculateRelaxationTime(haloDistributionNext_Ptr,
                                        haloDistributionPrevious_Ptr, iP);

      #pragma unroll
      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        collision.collideAndStream(haloDistributionNext_Ptr,
                                   haloDistributionPrevious_Ptr, iP, iQ);
      }

      if(Base::isStored) {
        Base::storeLocalFields(iP);
      }
    }

    HOST
    void iterate(const unsigned int iteration) {
      INSTRUMENT_ON("Algorithm<T, AlgorithmType::Pull>::iterate",2)

      std::swap(haloDistributionPrevious_Ptr, haloDistributionNext_Ptr);

      collision.update(iteration);

      auto t0 = std::chrono::high_resolution_clock::now();

      communication.communicateHalos(haloDistributionPrevious_Ptr);

      // TODO: Run only at the boundaries
      Base::computationBottom.Do(periodicBoundary.applyYBottom,
                                 haloDistributionPrevious_Ptr);
      Base::computationTop.Do(periodicBoundary.applyYTop,
                              haloDistributionPrevious_Ptr);
      Base::computationFront.Do(periodicBoundary.applyZFront,
                                haloDistributionPrevious_Ptr);
      Base::computationBack.Do(periodicBoundary.applyZBack,
                               haloDistributionPrevious_Ptr);

      //boundary.apply(f_Previous.haloData());

      auto t1 = std::chrono::high_resolution_clock::now();

      Base::computationLocal.Do(*this);

      auto t2 = std::chrono::high_resolution_clock::now();

      Base::dtCommunication = (t1 - t0);
      Base::dtComputation = (t2 - t1);
    }

    using Base::pack;
    using Base::unpack;
    using Base::getCommunicationTime;
    using Base::getComputationTime;
  };


}

#endif // ALGORITHM_H
