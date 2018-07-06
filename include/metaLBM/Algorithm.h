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
  public:
    T* densityPtr;
    T* velocityPtr;
    T* forcePtr;
    T* alphaPtr;
    T* T2Ptr;
    T* T3Ptr;
    T* T4Ptr;
    T* distributionPtr;

  protected:
    T* haloDistributionPreviousPtr;
    T* haloDistributionNextPtr;

    Packer<T> packer;
    Unpacker<T> unpacker;
    Communication<T, latticeT, algorithmT, memoryL, partitionningT,
                  implementation, L::dimD> communication;
    Collision_<architecture> collision;

    std::chrono::duration<double> dtComputation;
    std::chrono::duration<double> dtCommunication;

  public:
    bool isStored;

    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in,
              Communication<T, latticeT, algorithmT, memoryL, partitionningT,
                            implementation, L::dimD>& communication_in)
      : densityPtr(fieldList_in.density.getData(FFTWInit::numberElements))
      , velocityPtr(fieldList_in.velocity.getData(FFTWInit::numberElements))
      , forcePtr(fieldList_in.force.getData(FFTWInit::numberElements))
      , alphaPtr(fieldList_in.alpha.getData(FFTWInit::numberElements))
      , T2Ptr(fieldList_in.T2.getData(FFTWInit::numberElements))
      , T3Ptr(fieldList_in.T3.getData(FFTWInit::numberElements))
      , T4Ptr(fieldList_in.T4.getData(FFTWInit::numberElements))
      , distributionPtr(distribution_in.getData(FFTWInit::numberElements))
      , haloDistributionPreviousPtr(distribution_in.getHaloDataPrevious())
      , haloDistributionNextPtr(distribution_in.getHaloDataNext())
      , communication(communication_in)
      , collision(relaxationTime, fieldList_in, forceAmplitude, forceWaveLength,
                  forcekMin, forcekMax)
      , dtComputation()
      , dtCommunication()
      , isStored(false)
    {}

    double getCommunicationTime() { return dtCommunication.count(); }

    double getComputationTime() { return dtComputation.count(); }

  protected:
    LBM_DEVICE LBM_HOST void storeFields(const Position& iP,
                                         const unsigned int numberElements) {
      LBM_INSTRUMENT_OFF("Algorithm<T, AlgorithmType::Pull>::storeFields", 4)

        const auto indexLocal = hSD::getIndexLocal(iP);

      densityPtr[indexLocal] = collision.getDensity();
      for(auto iD = 0; iD < L::dimD; ++iD) {
        (velocityPtr + iD * numberElements)[indexLocal] =
          collision.getHydrodynamicVelocity()[iD];
      }

      if(writeAlpha) {
        alphaPtr[indexLocal] = collision.getAlpha();
      }

      if(writeT) {
        T2Ptr[indexLocal] = collision.getT2();
        T3Ptr[indexLocal] = collision.getT3();
        T4Ptr[indexLocal] = collision.getT4();
      }

      if(writeForce) {
        for(auto iD = 0; iD < L::dimD; ++iD) {
          (forcePtr + iD * numberElements)[indexLocal] =
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
              Communication<T, latticeT, algorithmT, memoryL, partitionningT,
              implementation, L::dimD>& communication_in)
      : Base(fieldList_in, distribution_in, communication_in)
      , computationLocal(L::halo(), lSD::sEnd() + L::halo(), {d::X, d::Y, d::Z})
      , computationBottom({hSD::start()[d::X], L::halo()[d::Y], hSD::start()[d::Z]},
          {hSD::end()[d::X], 2 * L::halo()[d::Y], hSD::end()[d::Z]},
          {d::X, d::Y, d::Z})
      , computationTop({hSD::start()[d::X], lSD::sLength()[d::Y], hSD::start()[d::Z]},
          {hSD::end()[d::X], L::halo()[d::Y] + lSD::sLength()[d::Y], hSD::end()[d::Z]},
          {d::X, d::Y, d::Z})
      , computationFront({hSD::start()[d::X], hSD::start()[d::Y], L::halo()[d::Z]},
          {hSD::end()[d::X], hSD::end()[d::Y], 2 * L::halo()[d::Z]},
          {d::X, d::Y, d::Z})
      , computationBack({hSD::start()[d::X], hSD::start()[d::Y], lSD::sLength()[d::Z]},
          {hSD::end()[d::X], hSD::end()[d::Y], L::halo()[d::Z] + lSD::sLength()[d::Z]},
          {d::X, d::Y, d::Z})
    {}

    LBM_DEVICE void operator()(const Position& iP, const unsigned int numberElements,
                               const MathVector<int, 3> rank) {
      Base::collision.calculateMoments(Base::haloDistributionPreviousPtr, iP);

      Base::collision.setForce(Base::forcePtr, iP, gSD::sOffset(rank), numberElements);
      Base::collision.calculateRelaxationTime(Base::haloDistributionNextPtr,
                                              Base::haloDistributionPreviousPtr, iP);

      if (Base::isStored) {
        Base::collision.calculateTs(Base::haloDistributionPreviousPtr,
                                    Base::haloDistributionNextPtr, iP);
      }

      #pragma unroll
      for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
        Base::collision.collideAndStream(Base::haloDistributionNextPtr,
                                         Base::haloDistributionPreviousPtr, iP,
                                         iQ);
      }

      if (Base::isStored) {
        Base::storeFields(iP, numberElements);
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
      LBM_INSTRUMENT_ON("Algorithm<T, AlgorithmType::Pull>::iterate", 2)

      std::swap(Base::haloDistributionPreviousPtr, Base::haloDistributionNextPtr);

      Base::collision.update(iteration, FFTWInit::numberElements);

      auto t0 = Clock::now();
      Base::communication.communicateHalos(Base::haloDistributionPreviousPtr);
      computationBottom.Do(defaultStream, bottomBoundary,
                           Base::haloDistributionPreviousPtr);
      computationTop.Do(defaultStream, topBoundary,
                        Base::haloDistributionPreviousPtr);
      computationFront.Do(defaultStream, frontBoundary,
                          Base::haloDistributionPreviousPtr);
      computationBack.Do(defaultStream, backBoundary,
                         Base::haloDistributionPreviousPtr);
      auto t1 = Clock::now();
      Base::dtCommunication = (t1 - t0);

      computationLocal.Do(defaultStream, *this, FFTWInit::numberElements, MPIInit::rank);
      computationLocal.synchronize();
      t0 = Clock::now();
      Base::dtComputation = (t0 - t1);
    }

    LBM_HOST
    void pack(const Stream<architecture>& stream) {
      computationLocal.Do(stream, Base::packer, Base::distributionPtr,
                          Base::haloDistributionNextPtr, FFTWInit::numberElements);
      computationLocal.synchronize();
    }

    LBM_HOST
    void unpack(const Stream<architecture>& stream) {
      computationLocal.Do(stream, Base::unpacker, Base::haloDistributionNextPtr,
                          Base::distributionPtr, FFTWInit::numberElements);
      computationLocal.synchronize();
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
              Communication<T, latticeT, algorithmT, memoryL, partitionningT,
              implementation, L::dimD>& communication_in)
      : Base(fieldList_in, distribution_in, communication_in)
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
      LBM_INSTRUMENT_ON("Algorithm<T, AlgorithmType::Pull>::iterate", 2)

      std::swap(Base::haloDistributionPreviousPtr,
                Base::haloDistributionNextPtr);

      Base::collision.update(iteration);

      auto t0 = Clock::now();
      Base::computationBottom.Do(bulkStream, Base::bottomBoundary,
                                 Base::haloDistributionPreviousPtr);
      Base::computationTop.Do(bulkStream, Base::topBoundary,
                              Base::haloDistributionPreviousPtr);
      Base::computationFront.Do(bulkStream, Base::frontBoundary,
                                Base::haloDistributionPreviousPtr);
      Base::computationBack.Do(bulkStream, Base::backBoundary,
                               Base::haloDistributionPreviousPtr);

      bulkStream.synchronize();
      //leftEvent.record(bulkStream);
      //rightEvent.record(bulkStream);

      auto t1 = Clock::now();
      Base::dtCommunication = (t1 - t0);

      computationBulk.Do(bulkStream, *this, FFTWInit::numberElements, MPIInit::rank);
      t0 = Clock::now();
      Base::dtComputation = (t0 - t1);

      // TODO: if only 1 GPU, use leftBoundary and rightBoundary instead of MPI
      Base::communication.communicateHalos(Base::haloDistributionPreviousPtr);

      //leftEvent.wait(leftStream);
      //rightEvent.wait(rightStream);

      t1 = Clock::now();
      Base::dtCommunication += (t1 - t0);

      computationLeft.Do(leftStream, *this, FFTWInit::numberElements, MPIInit::rank);
      computationRight.Do(rightStream, *this, FFTWInit::numberElements, MPIInit::rank);

      bulkStream.synchronize();
      leftStream.synchronize();
      rightStream.synchronize();

      t0 = Clock::now();
      Base::dtComputation += (t0 - t1);

    }

    using Base::unpack;

  };

#ifdef USE_NVSHMEM

  template <class T, Overlapping overlapping>
  class Algorithm<T, AlgorithmType::Pull, Architecture::GPU,
                  Implementation::NVSHMEM_IN, overlapping>
    : public Algorithm<T, AlgorithmType::Generic, Architecture::GPU,
                       Implementation::MPI, Overlapping::Off> {
  private:
    using Base =
      Algorithm<T, AlgorithmType::Generic, Architecture::GPU,
                Implementation::MPI, Overlapping::Off>;
    using Clock = std::chrono::high_resolution_clock;

  public:
    using Base::Algorithm;

    LBM_DEVICE void operator()(const Position& iP, const unsigned int numberElements,
                               const MathVector<int, 3> rank) {
      if(iP[d::X] == Base::computationLocal.start[Base::computationLocal.dir[d::X]]) {
        for (auto iQ = 1; iQ < L::faceQ + 1; ++iQ) {
          Base::haloDistributionPtr[hSD::getIndex(iP_Destination, iQ)] =
            haloDistributionLeftPtr[hSD::getIndex(iP, iQ)];
        }
      }
      else if(iP[d::X] == Base::computationLocal.end[Base::computationLocal.dir[d::X]]) {
        for (auto iQ = 1; iQ < L::faceQ + 1; ++iQ) {
          Base::haloDistributionPtr[hSD::getIndex(iP_Destination, iQ)] =
            haloDistributionRightPtr[hSD::getIndex(iP, iQ)];
        }
      }

      else if(iP[d::Y] == Base::computationLocal.start[Base::computationLocal.dir[d::Y]]) {

      }

      else if(iP[d::Y] == Base::computationLocal.end[Base::computationLocal.dir[d::Y]]) {
      }

      else if(iP[d::Z] == Base::computationLocal.start[Base::computationLocal.dir[d::Z]]) {
      }

      else if(iP[d::Z] == Base::computationLocal.end[Base::computationLocal.dir[d::Z]]) {
      }

      // handle corners!!!!!

      Base::collision.calculateMoments(Base::haloDistributionPreviousPtr, iP);

      Base::collision.setForce(Base::forcePtr, iP, gSD::sOffset(rank), numberElements);
      Base::collision.calculateRelaxationTime(Base::haloDistributionNextPtr,
                                              Base::haloDistributionPreviousPtr, iP);

#pragma unroll
      for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
        Base::collision.collideAndStream(Base::haloDistributionNextPtr,
                                         Base::haloDistributionPreviousPtr, iP,
                                         iQ);
      }

      if (Base::isStored) {
        Base::storeFields(iP);
      }
    }

    LBM_HOST
    void iterate(const unsigned int iteration,
                 Stream<Architecture::GPU>& defaultStream,
                 Stream<Architecture::GPU>& bulkStream,
                 Stream<Architecture::GPU>& leftStream,
                 Stream<Architecture::GPU>& rightStream,
                 Event<Architecture::GPU>& leftEvent,
                 Event<Architecture::GPU>& rightEvent) {
      LBM_INSTRUMENT_ON("Algorithm<T, AlgorithmType::Pull>::iterate", 2)

        std::swap(Base::haloDistributionPreviousPtr,
                  Base::haloDistributionNextPtr);

      Base::collision.update(iteration, FFTWInit::numberElements);

      auto t0 = Clock::now();
      Base::computationBottom.Do(defaultStream, Base::bottomBoundary,
                           Base::haloDistributionPreviousPtr);
      Base::computationTop.Do(defaultStream, Base::topBoundary,
                        Base::haloDistributionPreviousPtr);
      Base::computationFront.Do(defaultStream, Base::frontBoundary,
                          Base::haloDistributionPreviousPtr);
      Base::computationBack.Do(defaultStream, Base::backBoundary,
                         Base::haloDistributionPreviousPtr);
      auto t1 = Clock::now();
      Base::dtCommunication = (t1 - t0);

      Base::computationLocal.Do(defaultStream, *this, FFTWInit::numberElements,
                                MPIInit::rank);
      Base::computationLocal::synchronize();

      t0 = Clock::now();
      Base::dtComputation = (t0 - t1);
    }

  };
#endif  // USE_NVSHMEM

}  // namespace lbm
