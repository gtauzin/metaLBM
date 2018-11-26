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

  template <class T, LatticeType latticeType, AlgorithmType algorithmType,
            MemoryLayout memoryLayout, PartitionningType partitionningType,
            CommunicationType communicationType, unsigned int Dimension>
  class Communication {};

  template <class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::Generic,
                      PartitionningType::Generic, CommunicationType::MPI, 0> {
  protected:
    Computation<Architecture::CPU, L::dimD> computationLocal;

    MPI_Status statusXRightMPI[2];
    MPI_Request requestXRightMPI[2];
    MPI_Status statusXLeftMPI[2];
    MPI_Request requestXLeftMPI[2];

    MPI_Status statusYMPI[4];
    MPI_Request requestYMPI[4];

    MPI_Status statusZMPI[4];
    MPI_Request requestZMPI[4];

    LBM_HOST
    Communication()
      : computationLocal(lSD::sStart(), lSD::sEnd())
      , statusXRightMPI()
      , requestXRightMPI()
      , statusXLeftMPI()
      , requestXLeftMPI()
      , statusYMPI()
      , requestYMPI()
      , statusZMPI()
      , requestZMPI()
    {}

    LBM_HOST
    void sendGlobalToLocal(T* globalPtr, T* localPtr,
                           unsigned int numberComponents) {
      LBM_INSTRUMENT_ON("Communication<6>::sendGlobalToLocal", 3)

        MPI_Scatter(globalPtr, numberComponents * lSD::pVolume(), MPI_DOUBLE,
                    localPtr, numberComponents * lSD::pVolume(), MPI_DOUBLE, 0,
                    MPI_COMM_WORLD);
    }

    LBM_HOST
    void sendLocalToGlobal(T* localPtr, T* globalPtr,
                           unsigned int numberComponents) {
      LBM_INSTRUMENT_ON("Communication<6>::sendLocalToGlobal", 3)

        MPI_Gather(localPtr, numberComponents * lSD::pVolume(), MPI_DOUBLE,
                   globalPtr, numberComponents * lSD::pVolume(), MPI_DOUBLE, 0,
                   MPI_COMM_WORLD);
    }

    LBM_HOST
    void reduce(T* localSumPtr, unsigned int numberComponents) {
      MPI_Barrier(MPI_COMM_WORLD);

      if (MPIInit::rank[d::X] == 0) {
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
      computationLocal.synchronize();

      reduce(&localSum, 1);
      return localSum;
    }
  };

  template <class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::SoA,
                      PartitionningType::Generic, CommunicationType::MPI, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::Generic,
                           PartitionningType::Generic, CommunicationType::MPI, 0> {
  protected:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::Generic,
                               PartitionningType::Generic, CommunicationType::MPI, 0>;

    using Base::requestXLeftMPI;
    using Base::requestXRightMPI;
    using Base::statusXLeftMPI;
    using Base::statusXRightMPI;

    using Base::requestYMPI;
    using Base::statusYMPI;

    using Base::requestZMPI;
    using Base::statusZMPI;

    typedef Domain<DomainType::HaloSpace, PartitionningType::Generic,
                   MemoryLayout::SoA, L::dimQ> hMLSD;

    unsigned int sizeStripeX;
    unsigned int sendToRightBeginX;
    unsigned int receivedFromLeftBeginX;
    unsigned int sendToLeftBeginX;
    unsigned int receivedFromRightBeginX;

  protected:
    LBM_HOST
    void sendAndReceiveHaloXRight(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON( "Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloXRight", 4)

      for(auto iQ = L::faceQ + 1; iQ < 2 * L::faceQ + 1; ++iQ) {
        sendToRightBeginX = hMLSD::getIndex(
          Position({lSD::sLength()[d::X],
                hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);
        receivedFromLeftBeginX = hMLSD::getIndex(
          Position({0, hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

        LBM_MPI_CALL(MPI_Irecv(haloDistributionPtr + receivedFromLeftBeginX, sizeStripeX,
                               MPI_DOUBLE, MPIInit::rankLeft, 17, MPI_COMM_WORLD,
                               &requestXRightMPI[0]));

        LBM_MPI_CALL(MPI_Isend(haloDistributionPtr + sendToRightBeginX, sizeStripeX,
                               MPI_DOUBLE, MPIInit::rankRight, 17, MPI_COMM_WORLD,
                               &requestXRightMPI[1]));

        MPI_Waitall(2, requestXRightMPI, statusXRightMPI);
      }
    }

    LBM_HOST
    void sendAndReceiveHaloXLeft(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloXLeft", 4)

      for(auto iQ = 1; iQ < L::faceQ + 1; ++iQ) {
        sendToLeftBeginX = hMLSD::getIndex(
          Position({L::halo()[d::X], hMLSD::start()[d::Y],
                hMLSD::start()[d::Z]}), iQ);

        receivedFromRightBeginX = hMLSD::getIndex(
          Position({L::halo()[d::X] + lSD::sLength()[d::X],
                    hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

        LBM_MPI_CALL(MPI_Irecv(haloDistributionPtr + receivedFromRightBeginX, sizeStripeX,
                               MPI_DOUBLE, MPIInit::rankRight, 23, MPI_COMM_WORLD,
                               &requestXLeftMPI[0]));

        LBM_MPI_CALL(MPI_Isend(haloDistributionPtr + sendToLeftBeginX, sizeStripeX,
                               MPI_DOUBLE, MPIInit::rankLeft, 23, MPI_COMM_WORLD,
                               &requestXLeftMPI[1]));

        MPI_Waitall(2, requestXLeftMPI, statusXLeftMPI);
      }
    }

    LBM_HOST
    void sendAndReceiveHaloYBottom(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloYBottom", 4)
        // TODO - PACK AND UNPACK
        }

    LBM_HOST
    void sendAndReceiveHaloYTop(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloYTop", 4)
        // TODO - PACK AND UNPACK
        }

    LBM_HOST
    void sendAndReceiveHaloZFront(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloZFront", 4)

        // TODO: PACK AND UNPACK
        }

    LBM_HOST
    void sendAndReceiveHaloZBack(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloZBack", 4)

        // TODO: PACK AND UNPACK
        }

  public:
    LBM_HOST
    Communication()
      : Base()
      , sizeStripeX(L::halo()[d::X] * hMLSD::volume() / hMLSD::length()[d::X])
      , sendToRightBeginX(0)
      , receivedFromLeftBeginX(0)
      , sendToLeftBeginX(0)
      , receivedFromRightBeginX(0)
    {}

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
  };

  template <class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::AoS,
                      PartitionningType::Generic, CommunicationType::MPI, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::Generic,
                           PartitionningType::Generic, CommunicationType::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::Generic,
                               PartitionningType::Generic, CommunicationType::MPI, 0>;

    using Base::requestXLeftMPI;
    using Base::requestXRightMPI;
    using Base::statusXLeftMPI;
    using Base::statusXRightMPI;

    using Base::requestYMPI;
    using Base::statusYMPI;

    using Base::requestZMPI;
    using Base::statusZMPI;

    typedef Domain<DomainType::HaloSpace, PartitionningType::Generic,
                   MemoryLayout::AoS, L::dimQ> hMLSD;

  public:
    LBM_HOST
    Communication()
      : Base()
      , sizeStripeX(L::dimQ * hMLSD::volume() * L::halo()[d::X] /
                  hMLSD::length()[d::X])
      , sendToRightBeginX(hMLSD::getIndex(
          Position({lSD::sLength()[d::X],
                    hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), 0))
      , receivedFromLeftBeginX(hMLSD::getIndex(
          Position({0, hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), 0))
      , sendToLeftBeginX( hMLSD::getIndex(
          Position({L::halo()[d::X], hMLSD::start()[d::Y],
                    hMLSD::start()[d::Z]}), 0))
      , receivedFromRightBeginX(hMLSD::getIndex(
          Position({L::halo()[d::X] + lSD::sLength()[d::X],
                    hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), 0))
    {}

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;

  protected:
    unsigned int sizeStripeX;
    unsigned int sendToRightBeginX;
    unsigned int receivedFromLeftBeginX;
    unsigned int sendToLeftBeginX;
    unsigned int receivedFromRightBeginX;

    LBM_HOST
    void sendAndReceiveHaloXRight(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloXRight", 4)

      MPI_Irecv(haloDistributionPtr + receivedFromLeftBeginX, sizeStripeX,
                MPI_DOUBLE, MPIInit::rankLeft, 17, MPI_COMM_WORLD,
                &requestXRightMPI[0]);

      MPI_Isend(haloDistributionPtr + sendToRightBeginX, sizeStripeX, MPI_DOUBLE,
                MPIInit::rankRight, 17, MPI_COMM_WORLD, &requestXRightMPI[1]);

      MPI_Waitall(2, requestXRightMPI, statusXRightMPI);
    }

    LBM_HOST
    void sendAndReceiveHaloXLeft(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloXLeft", 4)

      MPI_Irecv(haloDistributionPtr + receivedFromRightBeginX, sizeStripeX,
                MPI_DOUBLE, MPIInit::rankRight, 23, MPI_COMM_WORLD,
                &requestXLeftMPI[0]);

      MPI_Isend(haloDistributionPtr + sendToLeftBeginX, sizeStripeX, MPI_DOUBLE,
                MPIInit::rankLeft, 23, MPI_COMM_WORLD, &requestXLeftMPI[1]);

      MPI_Waitall(2, requestXLeftMPI, statusXLeftMPI);
    }

    LBM_HOST
    void sendAndReceiveHaloYBottom(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloYBottom", 4)

        // TODO - PACK AND UNPACK
        }

    LBM_HOST
    void sendAndReceiveHaloYTop(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloYTop", 4)

        // TODO - PACK AND UNPACK
        }

    LBM_HOST
    void sendAndReceiveHaloZFront(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloZFront", 4)

        // TODO: PACK AND UNPACK
        }

    LBM_HOST
    void sendAndReceiveHaloZBack(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloZBack", 4)

        // TODO: PACK AND UNPACK
        }
  };

  #ifdef USE_NVSHMEM

  template <class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::SoA,
                      PartitionningType::Generic, CommunicationType::NVSHMEM_OUT, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::SoA,
                           PartitionningType::Generic, CommunicationType::MPI, 0> {

  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::SoA,
                               PartitionningType::Generic, CommunicationType::MPI, 0>;

    typedef Domain<DomainType::HaloSpace, PartitionningType::Generic,
                   MemoryLayout::SoA, L::dimQ> hMLSD;

  protected:
    LBM_HOST void sendAndReceiveHaloXRight(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloXRight", 4)

        for (auto iQ = L::faceQ + 1; iQ < 2 * L::faceQ + 1; ++iQ) {
          Base::sendToRightBeginX = hMLSD::getIndex(
            Position({L::halo()[d::X] + lSD::sLength()[d::X] - 1,
                      hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);
          Base::receivedFromLeftBeginX = hMLSD::getIndex(
            Position({0, hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

          shmem_putmem(haloDistributionPtr + Base::receivedFromLeftBeginX,
                       haloDistributionPtr + Base::sendToRightBeginX,
                       Base::sizeStripeX * sizeof(T), MPIInit::rankRight);
          //shmem_quiet();
          shmem_barrier_all();
        }
    }

    LBM_HOST
    void sendAndReceiveHaloXLeft(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloXLeft", 4)

        for (auto iQ = 1; iQ < L::faceQ + 1; ++iQ) {
          Base::sendToLeftBeginX = hMLSD::getIndex(
            Position({L::halo()[d::X], hMLSD::start()[d::Y],
                      hMLSD::start()[d::Z]}), iQ);

          Base::receivedFromRightBeginX = hMLSD::getIndex(
            Position({L::halo()[d::X] + lSD::sLength()[d::X],
                      hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

          shmem_putmem(haloDistributionPtr + Base::receivedFromRightBeginX,
                       haloDistributionPtr + Base::sendToLeftBeginX,
                       Base::sizeStripeX * sizeof(T), MPIInit::rankLeft);
          //shmem_quiet();
          shmem_barrier_all();
        }
    }

    using Base::sendAndReceiveHaloYBottom;
    using Base::sendAndReceiveHaloYTop;
    using Base::sendAndReceiveHaloZBack;
    using Base::sendAndReceiveHaloZFront;

  public:
    using Base::Communication;

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
  };

  template <class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::AoS,
                      PartitionningType::Generic, CommunicationType::NVSHMEM_OUT, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::AoS,
                               PartitionningType::Generic, CommunicationType::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, MemoryLayout::AoS,
                               PartitionningType::Generic, CommunicationType::MPI, 0>;

    typedef Domain<DomainType::HaloSpace, PartitionningType::Generic,
                   MemoryLayout::AoS, L::dimQ> hMLSD;

    using Base::receivedFromLeftBeginX;
    using Base::receivedFromRightBeginX;
    using Base::sendToLeftBeginX;
    using Base::sendToRightBeginX;
    using Base::sizeStripeX;

  public:
    using Base::Communication;

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;

  protected:
    LBM_HOST void sendAndReceiveHaloXRight(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloX", 4)

      shmem_double_put(haloDistributionPtr + Base::receivedFromLeftBeginX,
                       haloDistributionPtr + Base::sendToRightBeginX,
                       Base::sizeStripeX, MPIInit::rankRight);
      shmem_barrier_all();
    }

    LBM_HOST void sendAndReceiveHaloXLeft(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloX", 4)

      shmem_double_put(haloDistributionPtr + Base::receivedFromRightBeginX,
                       haloDistributionPtr + Base::sendToLeftBeginX,
                       Base::sizeStripeX, MPIInit::rankLeft);
      shmem_barrier_all();
    }

    using Base::sendAndReceiveHaloYBottom;
    using Base::sendAndReceiveHaloYTop;
    using Base::sendAndReceiveHaloZBack;
    using Base::sendAndReceiveHaloZFront;
  };

  template <class T, LatticeType latticeType, MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                      PartitionningType::Generic, CommunicationType::NVSHMEM_IN, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                           PartitionningType::Generic, CommunicationType::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                               PartitionningType::Generic, CommunicationType::MPI, 0>;

  public:
    using Base::Communication;

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;

  protected:
    LBM_HOST void sendAndReceiveHaloXRight(T* haloDistributionPtr) {
    }

    LBM_HOST void sendAndReceiveHaloXLeft(T* haloDistributionPtr) {
    }
  };


  template <class T, LatticeType latticeType, MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                      PartitionningType::Generic, CommunicationType::Persistent_NVSHMEM_IN, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                           PartitionningType::Generic, CommunicationType::NVSHMEM_IN, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                               PartitionningType::Generic, CommunicationType::NVSHMEM_IN, 0>;

  public:
    using Base::Communication;

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;

  protected:
    LBM_HOST void sendAndReceiveHaloXRight(T* haloDistributionPtr) {
    }

    LBM_HOST void sendAndReceiveHaloXLeft(T* haloDistributionPtr) {
    }
  };

  #endif  // USE_NVSHMEM


  template <class T, LatticeType latticeType, MemoryLayout memoryLayout,
            CommunicationType communicationType, unsigned int Dimension>
  class Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                      PartitionningType::OneD, communicationType, Dimension>
    : public Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                           PartitionningType::Generic, communicationType, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                               PartitionningType::Generic, communicationType, 0>;

  public:
    using Base::Communication;

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;

    LBM_HOST
    inline void communicateHalos(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<6>::communicateHalos",3)

      Base::sendAndReceiveHaloXRight(haloDistributionPtr);
      Base::sendAndReceiveHaloXLeft(haloDistributionPtr);
    }
  };

  template <class T, LatticeType latticeType, MemoryLayout memoryLayout,
            CommunicationType communicationType, unsigned int Dimension>
  class Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                      PartitionningType::TwoD, communicationType, Dimension>
    : public Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                           PartitionningType::Generic, communicationType, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                               PartitionningType::Generic, communicationType, 0>;

  public:
    using Base::Communication;

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
  };

  template <class T, LatticeType latticeType, MemoryLayout memoryLayout,
            CommunicationType communicationType>
  class Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                      PartitionningType::TwoD, communicationType, 2>
    : public Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                           PartitionningType::Generic, communicationType, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                               PartitionningType::Generic, communicationType, 0>;

  public:
    using Base::Communication;

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;

    LBM_HOST
    inline void communicateHalos(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<6>::communicateHalos", 3)

      Base::sendAndReceiveHaloXRight(haloDistributionPtr);
      Base::sendAndReceiveHaloXLeft(haloDistributionPtr);
      Base::sendAndReceiveHaloYBottom(haloDistributionPtr);
      Base::sendAndReceiveHaloYTop(haloDistributionPtr);
    }
  };

  template <class T, LatticeType latticeType, MemoryLayout memoryLayout,
            CommunicationType communicationType>
  class Communication<T,  latticeType, AlgorithmType::Pull, memoryLayout,
                      PartitionningType::TwoD, communicationType, 3>
    : public Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                           PartitionningType::TwoD, communicationType, 2> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                               PartitionningType::TwoD, communicationType, 2>;

  public:
    using Base::Communication;

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;

    using Base::communicateHalos;
  };

  template <class T, LatticeType latticeType, MemoryLayout memoryLayout,
            CommunicationType communicationType, unsigned int Dimension>
  class Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                      PartitionningType::ThreeD, communicationType, Dimension>
    : public Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                           PartitionningType::Generic, communicationType, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                           PartitionningType::Generic, communicationType, 0>;

    using Base::Communication;

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;

  private:
    using Base::sendAndReceiveHaloX;
  };

  template <class T, LatticeType latticeType, MemoryLayout memoryLayout,
            CommunicationType communicationType>
  class Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                      PartitionningType::ThreeD, communicationType, 3>
    : public Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                           PartitionningType::Generic, communicationType, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull, memoryLayout,
                               PartitionningType::Generic, communicationType, 0>;

  public:
    using Base::Communication;

    using Base::reduce;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;

    LBM_HOST
    inline void communicateHalos(T* haloDistributionPtr) {
      LBM_INSTRUMENT_ON("Communication<6>::communicateHalos", 3)

      Base::sendAndReceiveHaloZFront(haloDistributionPtr);
      Base::sendAndReceiveHaloZBack(haloDistributionPtr);
      Base::sendAndReceiveHaloYBottom(haloDistributionPtr);
      Base::sendAndReceiveHaloYTop(haloDistributionPtr);
      Base::sendAndReceiveHaloXRight(haloDistributionPtr);
      Base::sendAndReceiveHaloXLeft(haloDistributionPtr);
    }
  };

  typedef Communication<dataT, latticeT, algorithmT, memoryL, partitionningT,
                        communicationT, L::dimD> Communication_;

}  // namespace lbm
