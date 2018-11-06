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
            MemoryLayout memoryLayout, PartitionningType partitionningType,
            CommunicationType communicationType, Overlapping overlapping>
  class Algorithm {
  public:
    LBM_HOST LBM_DEVICE void operator()(const Position iP);
  };

  template <class T, Architecture architecture, Overlapping overlapping,
            MemoryLayout memoryLayout>
  class Algorithm<T, AlgorithmType::Generic, architecture, memoryLayout,
                  PartitionningType::Generic, CommunicationType::Generic, overlapping> {
  public:
    T* densityPtr;
    T* velocityPtr;
    T* forcePtr;
    T* alphaPtr;
    T* numberIterationsPtr;
    T* T2Ptr;
    T* T3Ptr;
    T* T4Ptr;
    T* T2_approxPtr;
    T* T3_approxPtr;
    T* T4_approxPtr;
    T* fNeq_5_6_8Ptr;
    T* f1_5_6_8Ptr;
    T* piNeqDiagonalPtr;
    T* piNeqSymmetricPtr;
    T* pi1DiagonalPtr;
    T* pi1SymmetricPtr;
    T* squaredQContractedPiNeqPtr;
    T* cubedQContractedPiNeqPtr;
    T* squaredQContractedPi1Ptr;
    T* cubedQContractedPi1Ptr;
    T* distributionPtr;

  protected:
    T* haloDistributionPreviousPtr;
    T* haloDistributionNextPtr;

    Computation<architecture, L::dimD> computationLocal;

    Packer<T> packer;
    Unpacker<T> unpacker;
    Collision_<architecture> collision;

    std::chrono::duration<double> dtComputation;
    std::chrono::duration<double> dtCommunication;

  public:
    bool isStored;

    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in)
      : densityPtr(fieldList_in.density.getData(FFTWInit::numberElements))
      , velocityPtr(fieldList_in.velocity.getData(FFTWInit::numberElements))
      , forcePtr(fieldList_in.force.getData(FFTWInit::numberElements))
      , alphaPtr(fieldList_in.alpha.getData(FFTWInit::numberElements))
      , numberIterationsPtr(fieldList_in.numberIterations.getData(FFTWInit::numberElements))
      , T2Ptr(fieldList_in.T2.getData(FFTWInit::numberElements))
      , T3Ptr(fieldList_in.T3.getData(FFTWInit::numberElements))
      , T4Ptr(fieldList_in.T4.getData(FFTWInit::numberElements))
      , T2_approxPtr(fieldList_in.T2_approx.getData(FFTWInit::numberElements))
      , T3_approxPtr(fieldList_in.T3_approx.getData(FFTWInit::numberElements))
      , T4_approxPtr(fieldList_in.T4_approx.getData(FFTWInit::numberElements))
      , fNeq_5_6_8Ptr(fieldList_in.fNeq_5_6_8.getData(FFTWInit::numberElements))
      , f1_5_6_8Ptr(fieldList_in.f1_5_6_8.getData(FFTWInit::numberElements))
      , piNeqDiagonalPtr(fieldList_in.piNeqDiagonal.getData(FFTWInit::numberElements))
      , piNeqSymmetricPtr(fieldList_in.piNeqSymmetric.getData(FFTWInit::numberElements))
      , pi1DiagonalPtr(fieldList_in.pi1Diagonal.getData(FFTWInit::numberElements))
      , pi1SymmetricPtr(fieldList_in.pi1Symmetric.getData(FFTWInit::numberElements))
      , squaredQContractedPiNeqPtr(fieldList_in.squaredQContractedPiNeq.getData(FFTWInit::numberElements))
      , cubedQContractedPiNeqPtr(fieldList_in.cubedQContractedPiNeq.getData(FFTWInit::numberElements))
      , squaredQContractedPi1Ptr(fieldList_in.squaredQContractedPi1.getData(FFTWInit::numberElements))
      , cubedQContractedPi1Ptr(fieldList_in.cubedQContractedPi1.getData(FFTWInit::numberElements))
      , distributionPtr(distribution_in.getData(FFTWInit::numberElements))
      , haloDistributionPreviousPtr(distribution_in.getHaloDataPrevious())
      , haloDistributionNextPtr(distribution_in.getHaloDataNext())
      , computationLocal(L::halo(), lSD::sEnd() + L::halo(), {d::X, d::Y, d::Z})
      , collision(relaxationTime, fieldList_in, forceAmplitude, forceWaveLength,
                  forcekMin, forcekMax)
      , dtComputation()
      , dtCommunication()
      , isStored(false)
    {}

    LBM_DEVICE
    void operator()(const Position& iP, const unsigned int numberElements,
                    const MathVector<int, 3> rank) {
      collision.calculateMoments(haloDistributionPreviousPtr, iP);

      collision.setForce(forcePtr, iP, gSD::sOffset(rank), numberElements);
      collision.calculateRelaxationTime(haloDistributionNextPtr,
                                        haloDistributionPreviousPtr, iP,
                                        alphaPtr[hSD::getIndexLocal(iP)]);
      alphaPtr[hSD::getIndexLocal(iP)] = collision.getAlpha();

      if (isStored) {
        collision.calculateObservables(haloDistributionPreviousPtr,
                                       haloDistributionNextPtr, iP);
      }

#pragma unroll
      for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
        collision.collideAndStream(haloDistributionNextPtr,
                                   haloDistributionPreviousPtr, iP,
                                   iQ);
      }

      if (isStored) {
        storeFields(iP, numberElements);
      }
    }

    double getCommunicationTime() { return dtCommunication.count(); }

    double getComputationTime() { return dtComputation.count(); }

    LBM_HOST
    void pack(const Stream<architecture>& stream) {
      LBM_INSTRUMENT_OFF("Algorithm<T, AlgorithmType::Pull>::pack", 4)

        computationLocal.Do(stream, packer, distributionPtr,
                            haloDistributionNextPtr, FFTWInit::numberElements);
      computationLocal.synchronize();
    }

    LBM_HOST
    void unpack(const Stream<architecture>& stream) {
      LBM_INSTRUMENT_OFF("Algorithm<T, AlgorithmType::Pull>::unpack", 4)
        computationLocal.Do(stream, unpacker, haloDistributionNextPtr,
                            distributionPtr, FFTWInit::numberElements);
      computationLocal.synchronize();
    }

  protected:
    LBM_DEVICE LBM_HOST
    void storeFields(const Position& iP, const unsigned int numberElements) {
      LBM_INSTRUMENT_OFF("Algorithm<T, AlgorithmType::Pull>::storeFields", 4)

        const auto indexLocal = hSD::getIndexLocal(iP);

      densityPtr[indexLocal] = collision.getDensity();
      for(auto iD = 0; iD < L::dimD; ++iD) {
        (velocityPtr + iD * numberElements)[indexLocal] =
          collision.getHydrodynamicVelocity()[iD];
      }

      if(writeAlpha) {
        alphaPtr[indexLocal] = collision.getAlpha();
        numberIterationsPtr[indexLocal] = collision.getNumberIterations();
      }

      if(writeKinetics) {
        T2Ptr[indexLocal] = collision.getT2();
        T3Ptr[indexLocal] = collision.getT3();
        T4Ptr[indexLocal] = collision.getT4();

        T2_approxPtr[indexLocal] = collision.getT2_approx();
        T3_approxPtr[indexLocal] = collision.getT3_approx();
        T4_approxPtr[indexLocal] = collision.getT4_approx();

        for(auto iD = 0; iD < 3; ++iD) {
          (fNeq_5_6_8Ptr + iD * numberElements)[indexLocal] =
            collision.getFNeq_5_6_8()[iD];
          (f1_5_6_8Ptr + iD * numberElements)[indexLocal] =
            collision.getF1_5_6_8()[iD];
        }
        for(auto iD = 0; iD < L::dimD; ++iD) {
          (piNeqDiagonalPtr + iD * numberElements)[indexLocal] =
            collision.getPiNeqDiagonal()[iD];
          (pi1DiagonalPtr + iD * numberElements)[indexLocal] =
            collision.getPi1Diagonal()[iD];
        }
        for(auto iD = 0; iD < 2 * L::dimD - 3; ++iD) {
          (piNeqSymmetricPtr + iD * numberElements)[indexLocal] =
            collision.getPiNeqSymmetric()[iD];
          (pi1SymmetricPtr + iD * numberElements)[indexLocal] =
            collision.getPi1Symmetric()[iD];
        }

        squaredQContractedPiNeqPtr[indexLocal] = collision.getSquaredQContractedPiNeq();
        cubedQContractedPiNeqPtr[indexLocal] = collision.getCubedQContractedPiNeq();

        squaredQContractedPi1Ptr[indexLocal] = collision.getSquaredQContractedPi1();
        cubedQContractedPi1Ptr[indexLocal] = collision.getCubedQContractedPi1();
      }

      if(writeForce) {
        for(auto iD = 0; iD < L::dimD; ++iD) {
          (forcePtr + iD * numberElements)[indexLocal] =
            collision.getForce()[iD];
        }
      }
    }
  };

  template <class T, Architecture architecture, MemoryLayout memoryLayout>
  class Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                  PartitionningType::OneD, CommunicationType::Generic, Overlapping::Off>
    : public Algorithm<T, AlgorithmType::Generic, architecture, memoryLayout,
                       PartitionningType::Generic, CommunicationType::Generic,
                       Overlapping::Off> {
  private:
    using Base =
      Algorithm<T, AlgorithmType::Generic, architecture, memoryLayout,
                PartitionningType::Generic, CommunicationType::Generic, Overlapping::Off>;
    using Clock = std::chrono::high_resolution_clock;

  protected:
    Computation<architecture, L::dimD> computationBottom;
    Computation<architecture, L::dimD> computationTop;
    Computation<architecture, L::dimD> computationFront;
    Computation<architecture, L::dimD> computationBack;
    using Base::computationLocal;

    BottomBoundary<T, BoundaryType::Periodic, AlgorithmType::Pull, partitionningT,
                   CommunicationType::MPI, L::dimD> bottomBoundary;
    TopBoundary<T,BoundaryType::Periodic, AlgorithmType::Pull, partitionningT,
                CommunicationType::MPI, L::dimD> topBoundary;
    FrontBoundary<T, BoundaryType::Periodic, AlgorithmType::Pull, partitionningT,
                  CommunicationType::MPI, L::dimD> frontBoundary;
    BackBoundary<T, BoundaryType::Periodic, AlgorithmType::Pull, partitionningT,
                 CommunicationType::MPI, L::dimD> backBoundary;

  public:
    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in)
      : Base(fieldList_in, distribution_in)
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

    using Base::pack;
    using Base::unpack;

  };

  template <class T, Architecture architecture, MemoryLayout memoryLayout>
  class Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                  PartitionningType::OneD, CommunicationType::Generic, Overlapping::On>
    : public Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                       PartitionningType::OneD, CommunicationType::Generic,
                       Overlapping::Off> {
  private:
    using Base =
      Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                PartitionningType::OneD, CommunicationType::Generic, Overlapping::Off>;
    using Clock = std::chrono::high_resolution_clock;

  protected:
    Computation<architecture, L::dimD> computationBulk;
    Computation<architecture, L::dimD> computationRight;
    Computation<architecture, L::dimD> computationLeft;

    using Base::computationLocal;
    using Base::computationBottom;
    using Base::computationTop;
    using Base::computationFront;
    using Base::computationBack;

    using Base::bottomBoundary;
    using Base::topBoundary;
    using Base::frontBoundary;
    using Base::backBoundary;

  public:
    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in)
      : Base(fieldList_in, distribution_in)
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

    using Base::pack;
    using Base::unpack;
  };


  template <class T, Architecture architecture, MemoryLayout memoryLayout>
  class Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                  PartitionningType::OneD, CommunicationType::MPI, Overlapping::Off>
    : public Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                       PartitionningType::OneD, CommunicationType::Generic,
                       Overlapping::Off> {
  private:
    using Base =
      Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                PartitionningType::OneD, CommunicationType::Generic, Overlapping::Off>;
    using Clock = std::chrono::high_resolution_clock;

  protected:
    Communication<T, L::Type, AlgorithmType::Pull, memoryLayout, PartitionningType::OneD,
                  CommunicationType::MPI, L::dimD> communication;

  public:
    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in,
              Communication<T, L::Type, AlgorithmType::Pull, memoryLayout,
              PartitionningType::OneD, CommunicationType::MPI,
              L::dimD>& communication_in)
      : Base(fieldList_in, distribution_in)
      , communication(communication_in)
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

      std::swap(Base::haloDistributionPreviousPtr, Base::haloDistributionNextPtr);

      Base::collision.update(iteration, FFTWInit::numberElements);

      auto t0 = Clock::now();
      communication.communicateHalos(Base::haloDistributionPreviousPtr);

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
      Base::computationLocal.Do(defaultStream, *this, FFTWInit::numberElements, MPIInit::rank);
      Base::computationLocal.synchronize();

      t0 = Clock::now();
      Base::dtComputation = (t0 - t1);

    }

    using Base::pack;
    using Base::unpack;

  };


  template <class T, Architecture architecture, MemoryLayout memoryLayout>
  class Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                  PartitionningType::OneD, CommunicationType::MPI, Overlapping::On>
    : public Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                       PartitionningType::OneD, CommunicationType::Generic,
                       Overlapping::On> {
  private:
    using Base =
      Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                PartitionningType::OneD, CommunicationType::Generic, Overlapping::On>;
    using Clock = std::chrono::high_resolution_clock;

  protected:
    Communication<T, L::Type, AlgorithmType::Pull, memoryLayout, PartitionningType::OneD,
                  CommunicationType::MPI, L::dimD> communication;

  public:
    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in,
              Communication<T, L::Type, AlgorithmType::Pull, memoryLayout,
              PartitionningType::OneD, CommunicationType::MPI,
              L::dimD>& communication_in)
      : Base(fieldList_in, distribution_in)
      , communication(communication_in)
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

      Base::computationBulk.Do(bulkStream, *this, FFTWInit::numberElements, MPIInit::rank);
      t0 = Clock::now();
      Base::dtComputation = (t0 - t1);

      // TODO: if only 1 GPU, use leftBoundary and rightBoundary instead of MPI
      communication.communicateHalos(Base::haloDistributionPreviousPtr);

      //leftEvent.wait(leftStream);
      //rightEvent.wait(rightStream);

      t1 = Clock::now();
      Base::dtCommunication += (t1 - t0);

      Base::computationLeft.Do(leftStream, *this, FFTWInit::numberElements, MPIInit::rank);
      Base::computationRight.Do(rightStream, *this, FFTWInit::numberElements, MPIInit::rank);

      bulkStream.synchronize();
      leftStream.synchronize();
      rightStream.synchronize();

      t0 = Clock::now();
      Base::dtComputation += (t0 - t1);

    }

    using Base::pack;
    using Base::unpack;

  };

#ifdef USE_NVSHMEM

  template <class T, Architecture architecture, MemoryLayout memoryLayout>
  class Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                  PartitionningType::OneD, CommunicationType::NVSHMEM_OUT, Overlapping::Off>
    : public Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                       PartitionningType::OneD, CommunicationType::Generic,
                       Overlapping::Off> {
  private:
    using Base =
      Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                PartitionningType::OneD, CommunicationType::Generic, Overlapping::Off>;
    using Clock = std::chrono::high_resolution_clock;

  protected:
    Communication<T, L::Type, AlgorithmType::Pull, memoryLayout, PartitionningType::OneD,
                  CommunicationType::NVSHMEM_OUT, L::dimD> communication;

  public:
    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in,
              Communication<T, L::Type, AlgorithmType::Pull, memoryLayout,
              PartitionningType::OneD, CommunicationType::NVSHMEM_OUT,
              L::dimD>& communication_in)
      : Base(fieldList_in, distribution_in)
      , communication(communication_in)
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

        std::swap(Base::haloDistributionPreviousPtr, Base::haloDistributionNextPtr);

      Base::collision.update(iteration, FFTWInit::numberElements);

      auto t0 = Clock::now();
      communication.communicateHalos(Base::haloDistributionPreviousPtr);

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

      Base::computationLocal.Do(defaultStream, *this, FFTWInit::numberElements, MPIInit::rank);
      Base::computationLocal.synchronize();
      t0 = Clock::now();
      Base::dtComputation = (t0 - t1);
    }

    using Base::pack;
    using Base::unpack;

  };


  template <class T, Architecture architecture, MemoryLayout memoryLayout>
  class Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                  PartitionningType::OneD, CommunicationType::NVSHMEM_OUT, Overlapping::On>
    : public Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                       PartitionningType::OneD, CommunicationType::Generic,
                       Overlapping::On> {
  private:
    using Base =
      Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                PartitionningType::OneD, CommunicationType::Generic, Overlapping::Off>;
    using Clock = std::chrono::high_resolution_clock;

  protected:
    Communication<T, L::Type, AlgorithmType::Pull, memoryLayout, PartitionningType::OneD,
                  CommunicationType::NVSHMEM_OUT, L::dimD> communication;

  public:
    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in,
              Communication<T, L::Type, AlgorithmType::Pull, memoryLayout,
              PartitionningType::OneD, CommunicationType::NVSHMEM_OUT,
              L::dimD>& communication_in)
      : Base(fieldList_in, distribution_in)
      , communication(communication_in)
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

      Base::computationBulk.Do(bulkStream, *this, FFTWInit::numberElements, MPIInit::rank);
      t0 = Clock::now();
      Base::dtComputation = (t0 - t1);

      // TODO: if only 1 GPU, use leftBoundary and rightBoundary instead of MPI
      communication.communicateHalos(Base::haloDistributionPreviousPtr);

      //leftEvent.wait(leftStream);
      //rightEvent.wait(rightStream);

      t1 = Clock::now();
      Base::dtCommunication += (t1 - t0);

      Base::computationLeft.Do(leftStream, *this, FFTWInit::numberElements, MPIInit::rank);
      Base::computationRight.Do(rightStream, *this, FFTWInit::numberElements, MPIInit::rank);

      bulkStream.synchronize();
      leftStream.synchronize();
      rightStream.synchronize();

      t0 = Clock::now();
      Base::dtComputation += (t0 - t1);

    }

    using Base::pack;
    using Base::unpack;
  };


  template <class T, Architecture architecture, MemoryLayout memoryLayout>
  class Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                  PartitionningType::OneD, CommunicationType::NVSHMEM_IN, Overlapping::Off>
    : public Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                       PartitionningType::OneD, CommunicationType::Generic,
                       Overlapping::Off> {
  private:
    using Base =
      Algorithm<T, AlgorithmType::Pull, architecture, memoryLayout,
                PartitionningType::OneD, CommunicationType::Generic, Overlapping::Off>;
    using Clock = std::chrono::high_resolution_clock;

  protected:
    Communication<T, L::Type, AlgorithmType::Pull, memoryLayout, PartitionningType::OneD,
                  CommunicationType::NVSHMEM_IN, L::dimD> communication;

  public:
    Algorithm(FieldList<T, architecture>& fieldList_in,
              Distribution<T, architecture>& distribution_in,
              Communication<T, L::Type, AlgorithmType::Pull, memoryLayout,
              PartitionningType::OneD, CommunicationType::NVSHMEM_IN,
              L::dimD>& communication_in)
      : Base(fieldList_in, distribution_in)
      , communication(communication_in)
    {}

    LBM_DEVICE
    void operator()(const Position& iP, const unsigned int numberElements,
                    const MathVector<int, 3> rank, const int rankLeft, const int rankRight) {
      if(iP[d::X] < L::dimH) {
        Position iP_Destination, iP_Source;

        #pragma unroll
        for(auto iQ = L::faceQ + 1; iQ < L::faceQ + 1 + L::level()[0]; ++iQ) {
          iP_Destination = iP - uiL::celerity()[iQ];

          iP_Source =
            Position{{iP_Destination[d::X] < L::halo()[d::X]
                      ? iP_Destination[d::X] + lSD::sLength()[d::X]
                      : iP_Destination[d::X],
                      iP_Destination[d::Y] < L::halo()[d::Y]
                      ? iP_Destination[d::Y] + lSD::sLength()[d::Y]
                      : iP_Destination[d::Y],
                      iP_Destination[d::Z] < L::halo()[d::Z]
                      ? iP_Destination[d::Z] + lSD::sLength()[d::Z]
                      : iP_Destination[d::Z]}};

          Base::haloDistributionPtr[hSD::getIndex(iP_Destination, iQ)] =
            shmem_double_g(Base::haloDistributionPtr[hSD::getIndex(iP_Source, iQ)], rankLeft);
        }

        for(auto iH = 1; iH < L::dimH - iP[d::X]; ++iH) {
          #pragma unroll
          for(auto iQ = L::faceQ + 1 + L::level()[iH-1]; iQ < L::faceQ + 1 + L::level()[iH]; ++iQ) {
            iP_Destination = iP - uiL::celerity()[iQ];

            iP_Source =
              Position{{iP_Destination[d::X] < L::halo()[d::X]
                        ? iP_Destination[d::X] + lSD::sLength()[d::X]
                        : iP_Destination[d::X],
                        iP_Destination[d::Y] < L::halo()[d::Y]
                        ? iP_Destination[d::Y] + lSD::sLength()[d::Y]
                        : iP_Destination[d::Y],
                        iP_Destination[d::Z] < L::halo()[d::Z]
                        ? iP_Destination[d::Z] + lSD::sLength()[d::Z]
                        : iP_Destination[d::Z]}};

            Base::haloDistributionPtr[hSD::getIndex(iP_Destination, iQ)] =
              shmem_double_g(Base::haloDistributionPtr[hSD::getIndex(iP_Source, iQ)], rankLeft);
          }
        }
      }
      else if(iP[d::X] >  lSD::sLength()[d::X] - 1 - L::dimH) {
        Position iP_Destination, iP_Source;

        #pragma unroll
        for(auto iQ = 1; iQ < 1 + L::level()[0]; ++iQ) {
          iP_Destination = iP - uiL::celerity()[L::iQ_Top()[iQ]];

          iP_Source =
            Position{{iP_Destination[d::X] > lSD::sLength()[d::X] - 1
                      ? iP_Destination[d::X] - lSD::sLength()[d::X]
                      : iP_Destination[d::X],
                      iP_Destination[d::Y] > lSD::sLength()[d::Y] - 1
                      ? iP_Destination[d::Y] - lSD::sLength()[d::Y]
                      : iP_Destination[d::Y],
                      iP_Destination[d::Z] > lSD::sLength()[d::Z] - 1
                      ? iP_Destination[d::Z] - lSD::sLength()[d::Z]
                      : iP_Destination[d::Z]}};

          Base::haloDistributionPtr[hSD::getIndex(iP_Destination, iQ)] =
            shmem_double_g(Base::haloDistributionPtr[hSD::getIndex(iP_Source, iQ)], rankRight);
        }

        for(auto iH = 1; iH < lSD::sLength()[d::X] - 1 - L::dimH - iP[d::X]; ++iH) {

          #pragma unroll
          for(auto iQ = 1 + L::level()[iH-1]; iQ < 1 + L::level()[iH]; ++iQ) {
            iP_Destination = iP - uiL::celerity()[L::iQ_Top()[iQ]];

            iP_Source =
              Position{{iP_Destination[d::X] > lSD::sLength()[d::X] - 1
                        ? iP_Destination[d::X] - lSD::sLength()[d::X]
                        : iP_Destination[d::X],
                        iP_Destination[d::Y] > lSD::sLength()[d::Y] - 1
                        ? iP_Destination[d::Y] - lSD::sLength()[d::Y]
                        : iP_Destination[d::Y],
                        iP_Destination[d::Z] > lSD::sLength()[d::Z] - 1
                        ? iP_Destination[d::Z] - lSD::sLength()[d::Z]
                        : iP_Destination[d::Z]}};

            Base::haloDistributionPtr[hSD::getIndex(iP_Destination, iQ)] =
              shmem_double_g(Base::haloDistributionPtr[hSD::getIndex(iP_Source, iQ)], rankRight);
          }
        }
      }
      else {
        if(iP[d::X] < L::halo()[d::Y]) {
          Base::bottomBoundary(iP, Base::haloDistributionPtr);
        }
        else if(iP[d::X] >  lSD::sLength()[d::Y] - 1 - L::halo()[d::Y]) {
          Base::topBoundary(iP, Base::haloDistributionPtr);
        }
        else {
          if((iP[d::X] < L::halo()[d::Z])) {
            Base::frontBoundary(iP, Base::haloDistributionPtr);
          }
          else if(iP[d::X] >  lSD::sLength()[d::Z] - 1 - L::halo()[d::Z]) {
            Base::backBoundary(iP, Base::haloDistributionPtr);
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
                                MPIInit::rank, MPIInit::rankLeft, MPIInit::rankRight);
      Base::computationLocal::synchronize();

      auto t0 = Clock::now();
      Base::dtComputation = (t0 - t1);
    }

  };
#endif  // USE_NVSHMEM

}  // namespace lbm
