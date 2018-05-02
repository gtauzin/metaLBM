#pragma once

#include <chrono>
#include <cstdio>
#include <utility>

#include "Boundary.h"
#include "Collision.h"
#include "Commons.h"
#include "Communication.h"
#include "Computation.h"
#include "Event.h"
#include "Distribution.h"
#include "Domain.h"
#include "Field.h"
#include "FieldList.h"
#include "Lattice.h"
#include "Moment.h"
#include "Options.h"

namespace lbm {

  template <class T, AlgorithmType algorithmType, Architecture architecture,
            Implementation implementation, Overlapping overlapping>
  class Algorithm {
  public:
    LBM_HOST LBM_DEVICE void operator()(const Position iP);
  };

  template <class T, Architecture architecture, Implementation implementation,
            Overlapping overlapping>
  class Algorithm<T, AlgorithmType::Generic, architecture, implementation, overlapping> {
  protected:
    T* localDensity_Ptr;
    T* localVelocity_Ptr;
    T* localForce_Ptr;
    T* localAlpha_Ptr;

    T* localDistribution_Ptr;
    T* haloDistributionPrevious_Ptr;
    T* haloDistributionNext_Ptr;
    const unsigned int numberElements;

    Packer<T> packer;
    Unpacker<T> unpacker;
    Communication<T, latticeT, algorithmT, memoryL, partitionningT,
                  implementation, L::dimD> communication;
    Collision_ collision;

    std::chrono::duration<double> dtComputation;
    std::chrono::duration<double> dtCommunication;

  public:
    bool isStored;

    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in,
              const unsigned int numberElements_in,
              Communication<T, latticeT, algorithmT, memoryL, partitionningT,
              implementation, L::dimD>& communication_in)
      : localDensity_Ptr(fieldList_in.density.getLocalData())
      , localVelocity_Ptr(fieldList_in.velocity.getLocalData())
      , localForce_Ptr(fieldList_in.force.getLocalData())
      , localAlpha_Ptr(fieldList_in.alpha.getLocalData())
      , localDistribution_Ptr(distribution_in.getLocalData())
      , haloDistributionPrevious_Ptr(distribution_in.getHaloDataPrevious())
      , haloDistributionNext_Ptr(distribution_in.getHaloDataNext())
      , numberElements(numberElements_in)
      , communication(communication_in)
      , collision(relaxationTime,
                  forceAmplitude,
                  forceWaveLength,
                  forcekMin,
                  forcekMax)
      , dtComputation()
      , dtCommunication()
      , isStored(false)
    {}

    double getCommunicationTime() { return dtCommunication.count(); }

    double getComputationTime() { return dtComputation.count(); }

  protected:
    LBM_DEVICE LBM_HOST void storeLocalFields(const Position& iP) {
      LBM_INSTRUMENT_OFF(
                         "Algorithm<T, AlgorithmType::Pull>::storeLocalFields", 4)

        const auto indexLocal = hSD::getIndexLocal(iP);

      localDensity_Ptr[indexLocal] = collision.getDensity();
      for (auto iD = 0; iD < L::dimD; ++iD) {
        (localVelocity_Ptr + iD * numberElements)[indexLocal] =
          collision.getHydrodynamicVelocity()[iD];
      }

      if (writeAlpha)
        localAlpha_Ptr[indexLocal] = collision.getAlpha();

      if (writeForce) {
        for (auto iD = 0; iD < L::dimD; ++iD) {
          (localForce_Ptr + iD * numberElements)[indexLocal] =
            collision.getForce()[iD];
        }
      }
    }
  };

  template <class T, Architecture architecture, Implementation implementation>
  class Algorithm<T, AlgorithmType::Pull, architecture, implementation, Overlapping::Off>
    : public Algorithm<T, AlgorithmType::Generic, architecture,
                       implementation, Overlapping::Off> {
  private:
    using Base =
      Algorithm<T, AlgorithmType::Generic, architecture, implementation, Overlapping::Off>;
    using Clock = std::chrono::high_resolution_clock;

  protected:
    Computation<architecture, L::dimD> computationLocal;
    Computation<architecture, L::dimD> computationBottom;
    Computation<architecture, L::dimD> computationTop;
    Computation<architecture, L::dimD> computationFront;
    Computation<architecture, L::dimD> computationBack;

    BottomBoundary<T, BoundaryType::Periodic, AlgorithmType::Pull, partitionningT,
                   implementation, L::dimD> bottomBoundary;
    TopBoundary<T,BoundaryType::Periodic, AlgorithmType::Pull, partitionningT,
                implementation, L::dimD> topBoundary;
    FrontBoundary<T, BoundaryType::Periodic, AlgorithmType::Pull, partitionningT,
                  implementation, L::dimD> frontBoundary;
    BackBoundary<T, BoundaryType::Periodic, AlgorithmType::Pull, partitionningT,
                 implementation, L::dimD> backBoundary;

  public:
    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in,
              const unsigned int numberElements_in,
              Communication<T, latticeT, algorithmT, memoryL, partitionningT,
              implementation, L::dimD>& communication_in)
      : Base(fieldList_in, distribution_in, numberElements_in, communication_in)
      , computationLocal(lSD::sStart() + L::halo(),
                         lSD::sEnd() + L::halo(),
                         {d::X, d::Y, d::Z})
      , computationBottom({hSD::start()[d::X], L::halo()[d::Y], hSD::start()[d::Z]},
                          {hSD::end()[d::X], 2 * L::halo()[d::Y], hSD::end()[d::Z]},
                          {d::X, d::Z, d::Y})
      , computationTop({hSD::start()[d::X], L::halo()[d::Y] + lSD::sLength()[d::Y] - 1,
                           hSD::start()[d::Z]},
                       {hSD::end()[d::X], 2 * L::halo()[d::Y] + lSD::sLength()[d::Y] - 1,
                           hSD::end()[d::Z]},
                       {d::X, d::Z, d::Y})
      , computationFront({hSD::start()[d::X], hSD::start()[d::Y], L::halo()[d::Z]},
                         {hSD::end()[d::X], hSD::end()[d::Y], 2 * L::halo()[d::Z]},
                         {d::X, d::Y, d::Z})
      , computationBack({hSD::start()[d::X], hSD::start()[d::Y],
                            L::halo()[d::Z] + lSD::sLength()[d::Z] - 1},
                        {hSD::end()[d::X], hSD::end()[d::Y],
                            2 * L::halo()[d::Z] + lSD::sLength()[d::Z] - 1},
                        {d::X, d::Y, d::Z})
    {}

    LBM_DEVICE void operator()(const Position& iP) {
      Base::collision.calculateMoments(Base::haloDistributionPrevious_Ptr, iP);

      Base::collision.setForce(Base::localForce_Ptr, iP, Base::numberElements,
                               gSD::sOffset(Base::communication.rankMPI));
      Base::collision.calculateRelaxationTime(Base::haloDistributionNext_Ptr,
                                              Base::haloDistributionPrevious_Ptr, iP);

      #pragma unroll
      for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
        Base::collision.collideAndStream(Base::haloDistributionNext_Ptr,
                                         Base::haloDistributionPrevious_Ptr, iP,
                                         iQ);
      }

      if (Base::isStored) {
        Base::storeLocalFields(iP);
      }
    }

    LBM_HOST
    void iterate(const unsigned int iteration,
                 Stream<architecture>& defaultStream,
                 Stream<architecture>& bulkStream,
                 Stream<architecture>& leftStream,
                 Stream<architecture>& rightStream,
                 Event<architecture>& leftEvent,
                 Event<architecture>& rightEvent) {
      {LBM_INSTRUMENT_ON("Algorithm<T, AlgorithmType::Pull>::iterate", 2)}

      std::swap(Base::haloDistributionPrevious_Ptr,
                Base::haloDistributionNext_Ptr);

      Base::collision.update(iteration);

      auto t0 = Clock::now();
      Base::communication.communicateHalos(Base::haloDistributionPrevious_Ptr);
      computationBottom.Do(defaultStream, bottomBoundary,
                           Base::haloDistributionPrevious_Ptr);
      computationTop.Do(defaultStream, topBoundary,
                        Base::haloDistributionPrevious_Ptr);
      computationFront.Do(defaultStream, frontBoundary,
                          Base::haloDistributionPrevious_Ptr);
      computationBack.Do(defaultStream, backBoundary,
                         Base::haloDistributionPrevious_Ptr);
      auto t1 = Clock::now();
      Base::dtCommunication = (t1 - t0);

      computationLocal.Do(defaultStream, *this);
      t0 = Clock::now();
      Base::dtComputation = (t0 - t1);
    }

    LBM_HOST
    void pack(const Stream<architecture>& stream) {
      computationLocal.Do(stream, Base::packer, Base::localDistribution_Ptr,
                          Base::haloDistributionNext_Ptr, Base::numberElements);
    }

    LBM_HOST
    void unpack(const Stream<architecture>& stream) {
      computationLocal.Do(stream, Base::unpacker, Base::haloDistributionNext_Ptr,
                          Base::localDistribution_Ptr, Base::numberElements);
    }
  };


  template <class T, Architecture architecture, Implementation implementation>
  class Algorithm<T, AlgorithmType::Pull, architecture, implementation, Overlapping::On>
    : public Algorithm<T, AlgorithmType::Pull, architecture,
                       implementation, Overlapping::Off> {
  private:
    using Base =
      Algorithm<T, AlgorithmType::Pull, architecture, implementation, Overlapping::Off>;
    using Clock = std::chrono::high_resolution_clock;

    Computation<architecture, L::dimD> computationBulk;
    Computation<architecture, L::dimD> computationRight;
    Computation<architecture, L::dimD> computationLeft;

  public:
    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in,
              const unsigned int numberElements_in,
              Communication<T, latticeT, algorithmT, memoryL, partitionningT,
              implementation, L::dimD>& communication_in)
      : Base(fieldList_in, distribution_in, numberElements_in, communication_in)
      , computationBulk({lSD::sStart()[d::X]+2*L::halo()[d::X],
                          lSD::sStart()[d::Y]+L::halo()[d::Y],
                          lSD::sStart()[d::Z]+L::halo()[d::Z]},
                        {lSD::sEnd()[d::X], lSD::sEnd()[d::Y]+L::halo()[d::Y],
                          lSD::sEnd()[d::Z]+L::halo()[d::Z]},
                        {d::X, d::Y, d::Z})
      , computationLeft({lSD::sStart()[d::X]+L::halo()[d::X],
                          lSD::sStart()[d::Y]+L::halo()[d::Y],
                          lSD::sStart()[d::Z]+L::halo()[d::Z]},
                        {lSD::sStart()[d::X]+2*L::halo()[d::X],
                          lSD::sEnd()[d::Y]+L::halo()[d::Y],
                          lSD::sEnd()[d::Z]+L::halo()[d::Z]},
                        {d::X, d::Y, d::Z})
      , computationRight({lSD::sEnd()[d::X],
                          lSD::sStart()[d::Y]+L::halo()[d::Y],
                          lSD::sStart()[d::Z]+L::halo()[d::Z]},
                         {lSD::sEnd()[d::X]+L::halo()[d::X],
                          lSD::sEnd()[d::Y]+L::halo()[d::Y],
                          lSD::sEnd()[d::Z]+L::halo()[d::Z]},
                         {d::X, d::Y, d::Z})
    {}

    LBM_HOST
    void iterate(const unsigned int iteration,
                 Stream<architecture>& defaultStream,
                 Stream<architecture>& bulkStream,
                 Stream<architecture>& leftStream,
                 Stream<architecture>& rightStream,
                 Event<architecture>& leftEvent,
                 Event<architecture>& rightEvent) {
      { LBM_INSTRUMENT_ON("Algorithm<T, AlgorithmType::Pull>::iterate", 2) }

      std::swap(Base::haloDistributionPrevious_Ptr,
                Base::haloDistributionNext_Ptr);

      Base::collision.update(iteration);

      auto t0 = Clock::now();
      Base::computationBottom.Do(bulkStream, Base::bottomBoundary,
                                 Base::haloDistributionPrevious_Ptr);
      Base::computationTop.Do(bulkStream, Base::topBoundary,
                              Base::haloDistributionPrevious_Ptr);
      Base::computationFront.Do(bulkStream, Base::frontBoundary,
                                Base::haloDistributionPrevious_Ptr);
      Base::computationBack.Do(bulkStream, Base::backBoundary,
                               Base::haloDistributionPrevious_Ptr);

      leftEvent.record(bulkStream);
      rightEvent.record(bulkStream);

      auto t1 = Clock::now();
      Base::dtCommunication = (t1 - t0);

      computationBulk.Do(bulkStream, *this);
      t0 = Clock::now();
      Base::dtComputation = (t0 - t1);

      // TODO: if only 1 GPU, use leftBoundary and rightBoundary instead of MPI
      Base::communication.communicateHalos(Base::haloDistributionPrevious_Ptr);

      leftEvent.wait(leftStream);
      rightEvent.wait(rightStream);

      t1 = Clock::now();
      Base::dtCommunication += (t1 - t0);

      computationLeft.Do(leftStream, *this);
      computationRight.Do(rightStream, *this);

      bulkStream.synchronize();
      leftStream.synchronize();
      rightStream.synchronize();

      t0 = Clock::now();
      Base::dtComputation += (t0 - t1);

    }

    using Base::unpack;

  };

}  // namespace lbm
