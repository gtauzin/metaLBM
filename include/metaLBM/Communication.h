#pragma once

#include <mpi.h>
#include <string>

#ifdef USE_NVSHMEM
#include <shmem.h>
#include <shmemx.h>
#endif

#include "Boundary.h"
#include "Commons.h"
#include "Computation.h"
#include "Domain.h"
#include "DynamicArray.cuh"
#include "Lattice.h"
#include "Options.h"

namespace lbm {

template <class T,
          LatticeType latticeType,
          AlgorithmType algorithmType,
          MemoryLayout memoryLayout,
          PartitionningType partitionningType,
          Implementation implementation,
          unsigned int Dimension>
class Communication {};

template <class T, LatticeType latticeType>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    MemoryLayout::Generic,
                    PartitionningType::Generic,
                    Implementation::MPI,
                    0> {
 protected:
  const MathVector<int, 3> rankMPI;
  const MathVector<int, 3> sizeMPI;
  const std::string processorName;

  Computation<Architecture::CPU, L::dimD> computationLocal;
  const unsigned int rightXRankMPI;
  const unsigned int leftXRankMPI;
  MPI_Status statusXRightMPI[2];
  MPI_Request requestXRightMPI[2];
  MPI_Status statusXLeftMPI[2];
  MPI_Request requestXLeftMPI[2];

  const unsigned int rightYRankMPI;
  const unsigned int leftYRankMPI;
  MPI_Status statusYMPI[4];
  MPI_Request requestYMPI[4];

  const unsigned int rightZRankMPI;
  const unsigned int leftZRankMPI;
  MPI_Status statusZMPI[4];
  MPI_Request requestZMPI[4];

  LBM_HOST
  Communication(const MathVector<int, 3>& rankMPI_in,
                const MathVector<int, 3>& sizeMPI_in,
                const std::string& processorName_in)
      : rankMPI(rankMPI_in),
        sizeMPI(sizeMPI_in),
        processorName(processorName_in),
        computationLocal(lSD::sStart(), lSD::sEnd()),
        leftXRankMPI((rankMPI_in[d::X] + sizeMPI_in[d::X] - 1) %
                     sizeMPI_in[d::X]),
        rightXRankMPI((rankMPI_in[d::X] + 1) % sizeMPI_in[d::X]),
        statusXRightMPI(),
        requestXRightMPI(),
        statusXLeftMPI(),
        requestXLeftMPI(),
        leftYRankMPI((rankMPI_in[d::Y] + sizeMPI_in[d::Y] - 1) %
                     sizeMPI_in[d::Y]),
        rightYRankMPI((rankMPI_in[d::Y] + 1) % sizeMPI_in[d::Y]),
        statusYMPI(),
        requestYMPI(),
        leftZRankMPI((rankMPI_in[d::Z] + sizeMPI_in[d::Z] - 1) %
                     sizeMPI_in[d::Z]),
        rightZRankMPI((rankMPI_in[d::Z] + 1) % sizeMPI_in[d::Z]),
        statusZMPI(),
        requestZMPI() {}

  LBM_HOST
  void sendGlobalToLocal(T* globalPtr,
                         T* localPtr,
                         unsigned int numberComponents) {
    LBM_SCOREP_INSTRUMENT_ON("Communication<6>::sendGlobalToLocal", 3)

    MPI_Scatter(globalPtr, numberComponents * lSD::pVolume(), MPI_DOUBLE,
                localPtr, numberComponents * lSD::pVolume(), MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
  }

  LBM_HOST
  void sendLocalToGlobal(T* localPtr,
                         T* globalPtr,
                         unsigned int numberComponents) {
    LBM_SCOREP_INSTRUMENT_ON("Communication<6>::sendLocalToGlobal", 3)

    MPI_Gather(localPtr, numberComponents * lSD::pVolume(), MPI_DOUBLE,
               globalPtr, numberComponents * lSD::pVolume(), MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
  }

  LBM_HOST
  void reduce(T* localSumPtr, unsigned int numberComponents) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (rankMPI[d::X] == 0) {
      MPI_Reduce(MPI_IN_PLACE, localSumPtr, numberComponents, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(localSumPtr, localSumPtr, numberComponents, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }

  LBM_HOST
  T reduce(T* localPtr) {
    T localSum = (T)0;
    computationLocal.Do([&] LBM_HOST(const Position& iP) {
      localSum += localPtr[lSD::getIndex(iP)];
    });

    reduce(&localSum, 1);
    return localSum;
  }
};

template <class T, LatticeType latticeType>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    MemoryLayout::SoA,
                    PartitionningType::Generic,
                    Implementation::MPI,
                    0> : public Communication<T,
                                              latticeType,
                                              AlgorithmType::Pull,
                                              MemoryLayout::Generic,
                                              PartitionningType::Generic,
                                              Implementation::MPI,
                                              0> {
 private:
  using Base = Communication<T,
                             latticeType,
                             AlgorithmType::Pull,
                             MemoryLayout::Generic,
                             PartitionningType::Generic,
                             Implementation::MPI,
                             0>;
  using Base::leftXRankMPI;
  using Base::requestXLeftMPI;
  using Base::requestXRightMPI;
  using Base::rightXRankMPI;
  using Base::statusXLeftMPI;
  using Base::statusXRightMPI;

  using Base::leftYRankMPI;
  using Base::requestYMPI;
  using Base::rightYRankMPI;
  using Base::statusYMPI;

  using Base::leftZRankMPI;
  using Base::requestZMPI;
  using Base::rightZRankMPI;
  using Base::statusZMPI;

  typedef Domain<DomainType::HaloSpace,
                 PartitionningType::Generic,
                 MemoryLayout::SoA,
                 L::dimQ>
      hMLSD;

  unsigned int sizeStripeX;
  unsigned int sendToRightBeginX;
  unsigned int receivedFromLeftBeginX;
  unsigned int sendToLeftBeginX;
  unsigned int receivedFromRightBeginX;

 protected:
  LBM_HOST
  void sendAndReceiveHaloXRight(T* haloDistributionPtr) {
    {
      LBM_SCOREP_INSTRUMENT_ON(
          "Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloXRight", 4)
    }

    for (auto iQ = L::faceQ + 1; iQ < 2 * L::faceQ + 1; ++iQ) {
      sendToRightBeginX = hMLSD::getIndex(
          Position({L::halo()[d::X] + lSD::sLength()[d::X] - 1,
                    hMLSD::start()[d::Y], hMLSD::start()[d::Z]}),
          iQ);
      receivedFromLeftBeginX = hMLSD::getIndex(
          Position({0, hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

      MPI_Irecv(haloDistributionPtr + receivedFromLeftBeginX, sizeStripeX,
                MPI_DOUBLE, leftXRankMPI, 17, MPI_COMM_WORLD,
                &requestXRightMPI[0]);

      MPI_Isend(haloDistributionPtr + sendToRightBeginX, sizeStripeX,
                MPI_DOUBLE, rightXRankMPI, 17, MPI_COMM_WORLD,
                &requestXRightMPI[1]);

      MPI_Waitall(2, requestXRightMPI, statusXRightMPI);
    }
  }

  LBM_HOST
  void sendAndReceiveHaloXLeft(T* haloDistributionPtr) {
    {
      LBM_SCOREP_INSTRUMENT_ON(
          "Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloXLeft", 4)
    }

    for (auto iQ = 1; iQ < L::faceQ + 1; ++iQ) {
      sendToLeftBeginX =
          hMLSD::getIndex(Position({L::halo()[d::X], hMLSD::start()[d::Y],
                                    hMLSD::start()[d::Z]}),
                          iQ);

      receivedFromRightBeginX = hMLSD::getIndex(
          Position({L::halo()[d::X] + lSD::sLength()[d::X],
                    hMLSD::start()[d::Y], hMLSD::start()[d::Z]}),
          iQ);

      MPI_Irecv(haloDistributionPtr + receivedFromRightBeginX, sizeStripeX,
                MPI_DOUBLE, rightXRankMPI, 23, MPI_COMM_WORLD,
                &requestXLeftMPI[0]);

      MPI_Isend(haloDistributionPtr + sendToLeftBeginX, sizeStripeX, MPI_DOUBLE,
                leftXRankMPI, 23, MPI_COMM_WORLD, &requestXLeftMPI[1]);

      MPI_Waitall(2, requestXLeftMPI, statusXLeftMPI);
    }
  }

  LBM_HOST
  void sendAndReceiveHaloYBottom(T* haloDistributionPtr) {
    {
      LBM_SCOREP_INSTRUMENT_ON(
          "Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloYBottom", 4)
    }
    // TODO - PACK AND UNPACK
  }

  LBM_HOST
  void sendAndReceiveHaloYTop(T* haloDistributionPtr) {
    {
      LBM_SCOREP_INSTRUMENT_ON(
          "Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloYTop", 4)
    }
    // TODO - PACK AND UNPACK
  }

  LBM_HOST
  void sendAndReceiveHaloZFront(T* haloDistributionPtr) {
    {
      LBM_SCOREP_INSTRUMENT_ON(
          "Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloZFront", 4)
    }
    // TODO: PACK AND UNPACK
  }

  LBM_HOST
  void sendAndReceiveHaloZBack(T* haloDistributionPtr) {
    {
      LBM_SCOREP_INSTRUMENT_ON(
          "Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloZBack", 4)
    }
    // TODO: PACK AND UNPACK
  }

 public:
  LBM_HOST
  Communication(const MathVector<int, 3>& rankMPI_in,
                const MathVector<int, 3>& sizeMPI_in,
                const std::string& processorName_in)
      : Communication<T,
                      latticeType,
                      AlgorithmType::Pull,
                      MemoryLayout::Generic,
                      PartitionningType::Generic,
                      Implementation::MPI,
                      0>(rankMPI_in, sizeMPI_in, processorName_in),
        sizeStripeX(hMLSD::volume() / hMLSD::length()[d::X]),
        sendToRightBeginX(0),
        receivedFromLeftBeginX(0),
        sendToLeftBeginX(0),
        receivedFromRightBeginX(0) {}

  using Base::rankMPI;
  using Base::reduce;
  using Base::sendGlobalToLocal;
  using Base::sendLocalToGlobal;
  using Base::sizeMPI;
};

template <class T, LatticeType latticeType>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    MemoryLayout::AoS,
                    PartitionningType::Generic,
                    Implementation::MPI,
                    0> : public Communication<T,
                                              latticeType,
                                              AlgorithmType::Pull,
                                              MemoryLayout::Generic,
                                              PartitionningType::Generic,
                                              Implementation::MPI,
                                              0> {
 private:
  using Base = Communication<T,
                             latticeType,
                             AlgorithmType::Pull,
                             MemoryLayout::Generic,
                             PartitionningType::Generic,
                             Implementation::MPI,
                             0>;

  using Base::leftXRankMPI;
  using Base::requestXLeftMPI;
  using Base::requestXRightMPI;
  using Base::rightXRankMPI;
  using Base::statusXLeftMPI;
  using Base::statusXRightMPI;

  using Base::leftYRankMPI;
  using Base::requestYMPI;
  using Base::rightYRankMPI;
  using Base::statusYMPI;

  using Base::leftZRankMPI;
  using Base::requestZMPI;
  using Base::rightZRankMPI;
  using Base::statusZMPI;

  typedef Domain<DomainType::HaloSpace,
                 PartitionningType::Generic,
                 MemoryLayout::AoS,
                 L::dimQ>
      hMLSD;

  unsigned int sizeStripeX;
  unsigned int sendToRightBeginX;
  unsigned int receivedFromLeftBeginX;
  unsigned int sendToLeftBeginX;
  unsigned int receivedFromRightBeginX;

 public:
  LBM_HOST
  Communication(const MathVector<int, 3>& rankMPI_in,
                const MathVector<int, 3>& sizeMPI_in,
                const std::string& processorName_in)
      : Communication<T,
                      latticeType,
                      AlgorithmType::Pull,
                      MemoryLayout::Generic,
                      PartitionningType::Generic,
                      Implementation::MPI,
                      0>(rankMPI_in, sizeMPI_in, processorName_in),
        sizeStripeX(L::dimQ * hMLSD::volume() * L::halo()[d::X] /
                    hMLSD::length()[d::X]),
        sendToRightBeginX(hMLSD::getIndex(
            Position({L::halo()[d::X] + lSD::sLength()[d::X] - 1,
                      hMLSD::start()[d::Y], hMLSD::start()[d::Z]}),
            0)),
        receivedFromLeftBeginX(hMLSD::getIndex(
            Position({0, hMLSD::start()[d::Y], hMLSD::start()[d::Z]}),
            0)),
        sendToLeftBeginX(
            hMLSD::getIndex(Position({L::halo()[d::X], hMLSD::start()[d::Y],
                                      hMLSD::start()[d::Z]}),
                            0)),
        receivedFromRightBeginX(hMLSD::getIndex(
            Position({L::halo()[d::X] + lSD::sLength()[d::X],
                      hMLSD::start()[d::Y], hMLSD::start()[d::Z]}),
            0)) {}

  using Base::rankMPI;
  using Base::reduce;
  using Base::sendGlobalToLocal;
  using Base::sendLocalToGlobal;
  using Base::sizeMPI;

 protected:
  LBM_HOST
  void sendAndReceiveHaloXRight(T* haloDistributionPtr) {
    LBM_SCOREP_INSTRUMENT_ON(
        "Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloXRight", 4)

    MPI_Irecv(haloDistributionPtr + receivedFromLeftBeginX, sizeStripeX,
              MPI_DOUBLE, leftXRankMPI, 17, MPI_COMM_WORLD,
              &requestXRightMPI[0]);

    MPI_Isend(haloDistributionPtr + sendToRightBeginX, sizeStripeX, MPI_DOUBLE,
              rightXRankMPI, 17, MPI_COMM_WORLD, &requestXRightMPI[1]);

    MPI_Waitall(2, requestXRightMPI, statusXRightMPI);
  }

  LBM_HOST
  void sendAndReceiveHaloXLeft(T* haloDistributionPtr) {
    LBM_SCOREP_INSTRUMENT_ON(
        "Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloXLeft", 4)

    MPI_Irecv(haloDistributionPtr + receivedFromRightBeginX, sizeStripeX,
              MPI_DOUBLE, rightXRankMPI, 23, MPI_COMM_WORLD,
              &requestXLeftMPI[0]);

    MPI_Isend(haloDistributionPtr + sendToLeftBeginX, sizeStripeX, MPI_DOUBLE,
              leftXRankMPI, 23, MPI_COMM_WORLD, &requestXLeftMPI[1]);

    MPI_Waitall(2, requestXLeftMPI, statusXLeftMPI);
  }

  LBM_HOST
  void sendAndReceiveHaloYBottom(T* haloDistributionPtr) {
    {
      LBM_SCOREP_INSTRUMENT_ON(
          "Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloYBottom", 4)
    }
    // TODO - PACK AND UNPACK
  }

  LBM_HOST
  void sendAndReceiveHaloYTop(T* haloDistributionPtr) {
    {
      LBM_SCOREP_INSTRUMENT_ON(
          "Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloYTop", 4)
    }
    // TODO - PACK AND UNPACK
  }

  LBM_HOST
  void sendAndReceiveHaloZFront(T* haloDistributionPtr) {
    {
      LBM_SCOREP_INSTRUMENT_ON(
          "Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloZFront", 4)
    }
    // TODO: PACK AND UNPACK
  }

  LBM_HOST
  void sendAndReceiveHaloZBack(T* haloDistributionPtr) {
    {
      LBM_SCOREP_INSTRUMENT_ON(
          "Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloZBack", 4)
    }
    // TODO: PACK AND UNPACK
  }
};

#ifdef USE_NVSHMEM

template <class T, LatticeType latticeType>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    MemoryLayout::SoA,
                    PartitionningType::Generic,
                    Implementation::NVSHMEM_OUT,
                    0> : public Communication<T,
                                              latticeType,
                                              AlgorithmType::Pull,
                                              MemoryLayout::SoA,
                                              PartitionningType::Generic,
                                              Implementation::MPI,
                                              0> {
 private:
  using Base = Communication<T,
                             latticeType,
                             AlgorithmType::Pull,
                             MemoryLayout::SoA,
                             PartitionningType::Generic,
                             Implementation::MPI,
                             0>;

  typedef Domain<DomainType::HaloSpace,
                 PartitionningType::Generic,
                 MemoryLayout::SoA,
                 L::dimQ>
      hMLSD;

 protected:
  HOST void sendAndReceiveHaloXRight(T* haloDistributionPtr) {
    {
      INSTRUMENT_ON(
          "Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloXRight", 4)
    }

    for (auto iQ = faceQ + 1; iQ < 2 * L::faceQ + 1; ++iQ) {
      Base::sendToRightBeginX = hMLSD::getIndex(
          Position({L::halo()[d::X] + lSD::sLength()[d::X] - 1,
                    hMLSD::start()[d::Y], hMLSD::start()[d::Z]}),
          iQ);
      Base::receivedFromLeftBeginX = hMLSD::getIndex(
          Position({0, hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

      shmem_double_put(haloDistributionPtr + Base::receivedFromLeftBeginX,
                       haloDistributionPtr + Base::sendToRightBeginX,
                       Base::sizeStripeX, Base::rightXRankMPI);
      shmem_barrier_all();
    }
  }

  void sendAndReceiveHaloXLeft(T* haloDistributionPtr) {
    {
      INSTRUMENT_ON(
          "Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloXLeft", 4)
    }

    for (auto iQ = 1; iQ < L::faceQ + 1; ++iQ) {
      Base::sendToLeftBeginX =
          hMLSD::getIndex(Position({L::halo()[d::X], hMLSD::start()[d::Y],
                                    hMLSD::start()[d::Z]}),
                          iQ);

      Base::receivedFromRightBeginX = hMLSD::getIndex(
          Position({L::halo()[d::X] + lSD::sLength()[d::X],
                    hMLSD::start()[d::Y], hMLSD::start()[d::Z]}),
          iQ);

      shmem_double_put(haloDistributionPtr + Base::sendToLeftBeginX,
                       haloDistributionPtr + Base::receivedFromRightBeginX,
                       Base::sizeStripeX, Base::leftXRankMPI);
      shmem_barrier_all();
    }
  }

  using Base::sendAndReceiveHaloYBottom;
  using Base::sendAndReceiveHaloYTop;
  using Base::sendAndReceiveHaloZBack;
  using Base::sendAndReceiveHaloZFront;

 public:
  using Base::Communication;

  using Base::rankMPI;
  using Base::reduce;
  using Base::sendGlobalToLocal;
  using Base::sendLocalToGlobal;
  using Base::sizeMPI;
};

template <class T, LatticeType latticeType>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    MemoryLayout::AoS,
                    PartitionningType::Generic,
                    Implementation::NVSHMEM_OUT,
                    0> : public Communication<T,
                                              latticeType,
                                              AlgorithmType::Pull,
                                              MemoryLayout::AoS,
                                              PartitionningType::Generic,
                                              Implementation::MPI,
                                              0> {
 private:
  using Base = Communication<T,
                             latticeType,
                             AlgorithmType::Pull,
                             MemoryLayout::AoS,
                             PartitionningType::Generic,
                             Implementation::MPI,
                             0>;

  typedef Domain<DomainType::HaloSpace,
                 PartitionningType::Generic,
                 MemoryLayout::AoS,
                 L::dimQ>
      hMLSD;

  using Base::receivedFromLeftBeginX;
  using Base::receivedFromRightBeginX;
  using Base::sendToLeftBeginX;
  using Base::sendToRightBeginX;
  using Base::sizeStripeX;

  using Base::Communication;

  using Base::rankMPI;
  using Base::reduce;
  using Base::sendGlobalToLocal;
  using Base::sendLocalToGlobal;
  using Base::sizeMPI;

 protected:
  HOST void sendAndReceiveHaloXRight(T* haloDistributionPtr) {
    {INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloX",
                   4)}

    shmem_double_put(haloDistributionPtr + Base::receivedFromLeftBeginX,
                     haloDistributionPtr + Base::sendToRightBeginX,
                     Base::sizeStripeX, Base::rightXRankMPI);
    shmem_barrier_all();
  }

  HOST void sendAndReceiveHaloXLeft(T* haloDistributionPtr) {
    {INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloX",
                   4)}

    shmem_double_put(haloDistributionPtr + Base::sendToLeftBeginX,
                     haloDistributionPtr + Base::receivedFromRightBeginX,
                     Base::sizeStripeX, Base::leftXRankMPI);
    shmem_barrier_all();
  }

  using Base::sendAndReceiveHaloBackZ;
  using Base::sendAndReceiveHaloFrontZ;
  using Base::sendAndReceiveHaloYBottom;
  using Base::sendAndReceiveHaloYTop;
};

template <class T, LatticeType latticeType, MemoryLayout memoryLayout>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    memoryLayout,
                    PartitionningType::Generic,
                    Implementation::NVSHMEM_IN,
                    0> : public Communication<T,
                                              latticeType,
                                              AlgorithmType::Pull,
                                              MemoryLayout::Generic,
                                              PartitionningType::Generic,
                                              Implementation::MPI,
                                              0> {
 private:
  using Base = Communication<T,
                             latticeType,
                             AlgorithmType::Pull,
                             MemoryLayout::Generic,
                             PartitionningType::Generic,
                             Implementation::MPI,
                             0>;

 public:
  using Base::Communication;

  using Base::rankMPI;
  using Base::reduce;
  using Base::sendGlobalToLocal;
  using Base::sendLocalToGlobal;
  using Base::sizeMPI;

 protected:
  /* HOST */
  /* void sendAndReceiveHaloX(T * haloDistributionPtr) { */
  /* } */

  /* HOST */
  /* void sendAndReceiveHaloY(T * haloDistributionPtr) { */
  /* } */

  /* HOST */
  /* void sendAndReceiveHaloZ(T * haloDistributionPtr) { */
  /* } */
};
#endif  // USE_NVSHMEM

template <class T,
          LatticeType latticeType,
          MemoryLayout memoryLayout,
          unsigned int Dimension>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    memoryLayout,
                    PartitionningType::OneD,
                    Implementation::MPI,
                    Dimension>
    : public Communication<T,
                           latticeType,
                           AlgorithmType::Pull,
                           memoryLayout,
                           PartitionningType::Generic,
                           Implementation::MPI,
                           0> {
 private:
  using Base = Communication<T,
                             latticeType,
                             AlgorithmType::Pull,
                             memoryLayout,
                             PartitionningType::Generic,
                             Implementation::MPI,
                             0>;

 public:
  LBM_HOST
  Communication(const MathVector<int, 3>& rankMPI_in,
                const MathVector<int, 3>& sizeMPI_in,
                const std::string& processorName_in)
      : Communication<T,
                      latticeType,
                      AlgorithmType::Pull,
                      memoryLayout,
                      PartitionningType::Generic,
                      Implementation::MPI,
                      0>(rankMPI_in, sizeMPI_in, processorName_in) {}

  using Base::rankMPI;
  using Base::reduce;
  using Base::sendGlobalToLocal;
  using Base::sendLocalToGlobal;
  using Base::sizeMPI;

  LBM_HOST
  inline void communicateHalos(T* haloDistributionPtr) {
    {LBM_SCOREP_INSTRUMENT_ON(
        "Communication<6>::communicateHalos",
        3)} Base::sendAndReceiveHaloXRight(haloDistributionPtr);
    Base::sendAndReceiveHaloXLeft(haloDistributionPtr);
  }
};

template <class T,
          LatticeType latticeType,
          MemoryLayout memoryLayout,
          unsigned int Dimension>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    memoryLayout,
                    PartitionningType::TwoD,
                    Implementation::MPI,
                    Dimension>
    : public Communication<T,
                           latticeType,
                           AlgorithmType::Pull,
                           memoryLayout,
                           PartitionningType::Generic,
                           Implementation::MPI,
                           0> {
 private:
  using Base = Communication<T,
                             latticeType,
                             AlgorithmType::Pull,
                             memoryLayout,
                             PartitionningType::Generic,
                             Implementation::MPI,
                             0>;

 public:
  using Base::Communication;

  using Base::rankMPI;
  using Base::reduce;
  using Base::sendGlobalToLocal;
  using Base::sendLocalToGlobal;
  using Base::sizeMPI;
};

template <class T, LatticeType latticeType, MemoryLayout memoryLayout>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    memoryLayout,
                    PartitionningType::TwoD,
                    Implementation::MPI,
                    2> : public Communication<T,
                                              latticeType,
                                              AlgorithmType::Pull,
                                              memoryLayout,
                                              PartitionningType::Generic,
                                              Implementation::MPI,
                                              0> {
 private:
  using Base = Communication<T,
                             latticeType,
                             AlgorithmType::Pull,
                             memoryLayout,
                             PartitionningType::Generic,
                             Implementation::MPI,
                             0>;

 public:
  using Base::Communication;

  using Base::rankMPI;
  using Base::reduce;
  using Base::sendGlobalToLocal;
  using Base::sendLocalToGlobal;
  using Base::sizeMPI;

  LBM_HOST
  inline void communicateHalos(T* haloDistributionPtr) {
    LBM_SCOREP_INSTRUMENT_ON("Communication<6>::communicateHalos", 3)

    Base::sendAndReceiveHaloXRight(haloDistributionPtr);
    Base::sendAndReceiveHaloXLeft(haloDistributionPtr);
    Base::sendAndReceiveHaloYBottom(haloDistributionPtr);
    Base::sendAndReceiveHaloYTop(haloDistributionPtr);
  }
};

template <class T, LatticeType latticeType, MemoryLayout memoryLayout>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    memoryLayout,
                    PartitionningType::TwoD,
                    Implementation::MPI,
                    3> : public Communication<T,
                                              latticeType,
                                              AlgorithmType::Pull,
                                              memoryLayout,
                                              PartitionningType::TwoD,
                                              Implementation::MPI,
                                              2> {
 private:
  using Base = Communication<T,
                             latticeType,
                             AlgorithmType::Pull,
                             memoryLayout,
                             PartitionningType::TwoD,
                             Implementation::MPI,
                             2>;

 public:
  using Base::Communication;

  using Base::rankMPI;
  using Base::reduce;
  using Base::sendGlobalToLocal;
  using Base::sendLocalToGlobal;
  using Base::sizeMPI;

  using Base::communicateHalos;
};

template <class T,
          LatticeType latticeType,
          MemoryLayout memoryLayout,
          unsigned int Dimension>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    memoryLayout,
                    PartitionningType::ThreeD,
                    Implementation::MPI,
                    Dimension>
    : public Communication<T,
                           latticeType,
                           AlgorithmType::Pull,
                           memoryLayout,
                           PartitionningType::Generic,
                           Implementation::MPI,
                           0> {
 private:
  using Base = Communication<T,
                             latticeType,
                             AlgorithmType::Pull,
                             memoryLayout,
                             PartitionningType::Generic,
                             Implementation::MPI,
                             0>;

  using Base::Communication;

  using Base::rankMPI;
  using Base::reduce;
  using Base::sendGlobalToLocal;
  using Base::sendLocalToGlobal;
  using Base::sizeMPI;

 private:
  using Base::sendAndReceiveHaloX;
};

template <class T, LatticeType latticeType, MemoryLayout memoryLayout>
class Communication<T,
                    latticeType,
                    AlgorithmType::Pull,
                    memoryLayout,
                    PartitionningType::ThreeD,
                    Implementation::MPI,
                    3> : public Communication<T,
                                              latticeType,
                                              AlgorithmType::Pull,
                                              memoryLayout,
                                              PartitionningType::Generic,
                                              Implementation::MPI,
                                              0> {
 private:
  using Base = Communication<T,
                             latticeType,
                             AlgorithmType::Pull,
                             memoryLayout,
                             PartitionningType::Generic,
                             Implementation::MPI,
                             0>;

 public:
  using Base::Communication;

  using Base::rankMPI;
  using Base::reduce;
  using Base::sendGlobalToLocal;
  using Base::sendLocalToGlobal;
  using Base::sizeMPI;

  LBM_HOST
  inline void communicateHalos(T* haloDistributionPtr) {
    LBM_SCOREP_INSTRUMENT_ON("Communication<6>::communicateHalos", 3)

    Base::sendAndReceiveHaloZFront(haloDistributionPtr);
    Base::sendAndReceiveHaloZBack(haloDistributionPtr);
    Base::sendAndReceiveHaloYBottom(haloDistributionPtr);
    Base::sendAndReceiveHaloYTop(haloDistributionPtr);
    Base::sendAndReceiveHaloXRight(haloDistributionPtr);
    Base::sendAndReceiveHaloXLeft(haloDistributionPtr);
  }
};

typedef Communication<dataT,
                      latticeT,
                      algorithmT,
                      memoryL,
                      partitionningT,
                      implementationT,
                      L::dimD>
    Communication_;

}  // namespace lbm
