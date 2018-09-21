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
    T* T2_approxPtr;
    T* T3_approxPtr;
    T* T4_approxPtr;
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
      , T2_approxPtr(fieldList_in.T2_approx.getData(FFTWInit::numberElements))
      , T3_approxPtr(fieldList_in.T3_approx.getData(FFTWInit::numberElements))
      , T4_approxPtr(fieldList_in.T4_approx.getData(FFTWInit::numberElements))
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

        T2_approxPtr[indexLocal] = collision.getT2_approx();
        T3_approxPtr[indexLocal] = collision.getT3_approx();
        T4_approxPtr[indexLocal] = collision.getT4_approx();
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
      , computationBottom({hSD::start()[d::X], hSD::start()[d::Y], hSD::start()[d::Z]},
          {hSD::end()[d::X], L::halo()[d::Y], hSD::end()[d::Z]},
          {d::X, d::Y, d::Z})
      , computationTop({hSD::start()[d::X], L::halo()[d::Y] + lSD::sLength()[d::Y], hSD::start()[d::Z]},
          {hSD::end()[d::X], 2 * L::halo()[d::Y] + lSD::sLength()[d::Y], hSD::end()[d::Z]},
          {d::X, d::Y, d::Z})
      , computationFront({hSD::start()[d::X], hSD::start()[d::Y], hSD::start()[d::Z]},
          {hSD::end()[d::X], hSD::end()[d::Y], L::halo()[d::Z]},
          {d::X, d::Y, d::Z})
      , computationBack({hSD::start()[d::X], hSD::start()[d::Y],  L::halo()[d::Z] + lSD::sLength()[d::Z]},
          {hSD::end()[d::X], hSD::end()[d::Y], 2 * L::halo()[d::Z] + lSD::sLength()[d::Z]},
          {d::X, d::Y, d::Z})
    {}

    LBM_DEVICE void operator()(const Position& iP, const unsigned int numberElements,
                               const MathVector<int, 3> rank) {
      Base::collision.calculateMoments(Base::haloDistributionPreviousPtr, iP);

      Base::collision.setForce(Base::forcePtr, iP, gSD::sOffset(rank), numberElements);
      Base::collision.calculateRelaxationTime(Base::haloDistributionNextPtr,
                                              Base::haloDistributionPreviousPtr, iP,
					      Base::alphaPtr[hSD::getIndexLocal(iP)]);
      Base::alphaPtr[hSD::getIndexLocal(iP)] = Base::collision.getAlpha();

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

      Base::collision.update(iteration, FFTWInit::numberElements);

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

  template <class T>
  class Algorithm<T, AlgorithmType::Pull, Architecture::GPU,
    Implementation::NVSHMEM_IN, Overlapping::Off>
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
      if(iP[d::X] < L::dimH) {
        Position iP_Destination, iP_Source;

        #pragma unroll
        for(auto iQ = L::faceQ + 1; iQ < L::faceQ + 1 + L::level()[0]; ++iQ) {
          iP_Destination = iP - uiL::celerity()[iQ];

          iP_Source{iP_Destination[d::X] < L::halo()[d::X] ? iP_Destination[d::X] + lSD::sLength()[d::X]
                                                           : iP_Destination[d::X],
                    iP_Destination[d::Y] < L::halo()[d::Y] ? iP_Destination[d::Y] + lSD::sLength()[d::Y]
                                                           : iP_Destination[d::Y],
                    iP_Destination[d::Z] < L::halo()[d::Z] ? iP_Destination[d::Z] + lSD::sLength()[d::Z]
                                                           : iP_Destination[d::Z]};

          haloDistributionPtr[hSD::getIndex(iP_Destination, iQ)] =
            haloDistributionPtr[hSD::getIndex(iP_Source, iQ)];
        }
      }

      for(auto iH = 1; iH < L::dimH - iP[d::X]; ++iQ) {

          #pragma unroll
          for(auto iQ = L::faceQ + 1 + L::level()[iH-1]; iQ < L::faceQ + 1 + L::level()[iH]; ++iQ) {
              iP_Destination = iP - uiL::celerity()[iQ];

              iP_Source{iP_Destination[d::X] < L::halo()[d::X] ? iP_Destination[d::X] + lSD::sLength()[d::X]
                                                               : iP_Destination[d::X],
                        iP_Destination[d::Y] < L::halo()[d::Y] ? iP_Destination[d::Y] + lSD::sLength()[d::Y]
                                                               : iP_Destination[d::Y],
                        iP_Destination[d::Z] < L::halo()[d::Z] ? iP_Destination[d::Z] + lSD::sLength()[d::Z]
                                                               : iP_Destination[d::Z]};

            haloDistributionPtr[hSD::getIndex(iP_Destination, iQ)] =
              haloDistributionPtr[hSD::getIndex(iP_Source, iQ)];
          }
        }
      }
    }
    else if(iP[d::X] >  lSD::sLength()[d::X] - 1 - L::dimH) {
        Position iP_Destination, iP_Source;

        #pragma unroll
        for(auto iQ = 1; iQ < 1 + L::level()[0]; ++iQ) {
          iP_Destination = iP - uiL::celerity()[L::iQ_Top()[iQ]];

          iP_Source{iP_Destination[d::X] <= L::halo()[d::X] ? iP_Destination[d::X]
                                                            : iP_Destination[d::X] - lSD::sLength()[d::X],
                    iP_Destination[d::Y] <= L::halo()[d::Y] ? iP_Destination[d::Y]
                                                            : iP_Destination[d::Y] - lSD::sLength()[d::Y],
                    iP_Destination[d::Z] <= L::halo()[d::Z] ? iP_Destination[d::Z]
                                                            : iP_Destination[d::Z] - lSD::sLength()[d::Z]};

          haloDistributionPtr[hSD::getIndex(iP_Destination, iQ)] =
            haloDistributionPtr[hSD::getIndex(iP_Source, iQ)];
        }

        for(auto iH = 1; iH < lSD::sLength()[d::X] - 1 - L::dimH - iP[d::X]; ++iQ) {

          #pragma unroll
          for(auto iQ = 1 + L::level()[iH-1]; iQ < 1 + L::level()[iH]; ++iQ) {
            iP_Destination = iP - uiL::celerity()[L::iQ_Top()[iQ]];

            iP_Source{iP_Destination[d::X] <= L::halo()[d::X] ? iP_Destination[d::X]
                                                              : iP_Destination[d::X] - lSD::sLength()[d::X],
                      iP_Destination[d::Y] <= L::halo()[d::Y] ? iP_Destination[d::Y]
                                                              : iP_Destination[d::Y] - lSD::sLength()[d::Y],
                      iP_Destination[d::Z] <= L::halo()[d::Z] ? iP_Destination[d::Z]
                                                              : iP_Destination[d::Z] - lSD::sLength()[d::Z]};

            haloDistributionPtr[hSD::getIndex(iP_Destination, iQ)] =
              haloDistributionPtr[hSD::getIndex(iP_Source, iQ)];
          }

    }

      else {
        if(iP[d::X] < L::halo()[d::Y]) {
          Base::bottomBoundary(iP, haloDistributionPtr);
        }

        else if(iP[d::X] >  lSD::sLength()[d::Y] - 1 - L::halo()[d::Y]) {
          Base::topBoundary(iP, haloDistributionPtr);
        }

        else {
          if((iP[d::X] < L::halo()[d::Z])) {
            Base::frontBoundary(iP, haloDistributionPtr);
          }
          else if(iP[d::X] >  lSD::sLength()[d::Z] - 1 - L::halo()[d::Z]) {
            Base::backBoundary(iP, haloDistributionPtr);
          }
        }
      }

    Base(iP, numberElements, rank);

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

      auto t1 = Clock::now();

      Base::computationLocal.Do(defaultStream, *this, FFTWInit::numberElements,
                                MPIInit::rank);
      Base::computationLocal::synchronize();

      t0 = Clock::now();
      Base::dtComputation = (t0 - t1);
    }

  };
#endif  // USE_NVSHMEM

}  // namespace lbm
