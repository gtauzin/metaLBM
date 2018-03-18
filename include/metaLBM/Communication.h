#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <mpi.h>
#include <string>

#ifdef USE_NVSHMEM
  #include<shmem.h>
  #include<shmemx.h>
#endif

#include "Commons.h"
#include "Options.h"
#include "DynamicArray.cuh"
#include "Lattice.h"
#include "Domain.h"
#include "Boundary.h"

#ifdef USE_NVSHMEM
  #include <shmem.h>
  #include <shmemx.h>
#endif
#include "Computation.h"


namespace lbm {

  template<class T, LatticeType latticeType, AlgorithmType algorithmType,
           MemoryLayout memoryLayout, PartitionningType partitionningType,
           Implementation implementation, unsigned int Dimension>
  class Communication {};


  template<class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::Generic, PartitionningType::Generic,
                      Implementation::MPI, 0> {
  protected:
    const MathVector<int, 3> rankMPI;
    const MathVector<int, 3> sizeMPI;
    const std::string processorName;

    const unsigned int rightXRankMPI;
    const unsigned int leftXRankMPI;
    MPI_Status statusXMPI[4];
    MPI_Request requestXMPI[4];

    const unsigned int rightYRankMPI;
    const unsigned int leftYRankMPI;
    MPI_Status statusYMPI[4];
    MPI_Request requestYMPI[4];

    const unsigned int rightZRankMPI;
    const unsigned int leftZRankMPI;
    MPI_Status statusZMPI[4];
    MPI_Request requestZMPI[4];

    HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : rankMPI(rankMPI_in)
      , sizeMPI(sizeMPI_in)
      , processorName(processorName_in)
      , leftXRankMPI((rankMPI_in[d::X] + sizeMPI_in[d::X] - 1) % sizeMPI_in[d::X])
      , rightXRankMPI((rankMPI_in[d::X] + 1) % sizeMPI_in[d::X])
      , statusXMPI()
      , requestXMPI()
      , leftYRankMPI((rankMPI_in[d::Y] + sizeMPI_in[d::Y] - 1) % sizeMPI_in[d::Y])
      , rightYRankMPI((rankMPI_in[d::Y] + 1) % sizeMPI_in[d::Y])
      , statusYMPI()
      , requestYMPI()
      , leftZRankMPI((rankMPI_in[d::Z] + sizeMPI_in[d::Z] - 1) % sizeMPI_in[d::Z])
      , rightZRankMPI((rankMPI_in[d::Z] + 1) % sizeMPI_in[d::Z])
      , statusZMPI()
      , requestZMPI()
    {}

    void printInputs() {
      std::cout << "MPI #" << rankMPI << " of " << sizeMPI
                << " running on host " << processorName << std::endl
                << "Right MPI #" << rightXRankMPI
                << ", left MPI #" << leftXRankMPI << std::endl;
    }

    DEVICE HOST
    MathVector<int, 3> getRankMPI() {
      return rankMPI;
    }

    HOST
    void sendGlobalToLocal(T * globalPtr,
                           T * localPtr,
                           unsigned int numberComponents) {
      INSTRUMENT_ON("Communication<6>::sendGlobalToLocal",3)

        MPI_Scatter(globalPtr, numberComponents*lSD::pVolume(), MPI_DOUBLE,
                    localPtr, numberComponents*lSD::pVolume(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

    }

    HOST
    void sendLocalToGlobal(T * localPtr,
                           T * globalPtr,
                           unsigned int numberComponents) {
      INSTRUMENT_ON("Communication<6>::sendLocalToGlobal",3)

      MPI_Gather(localPtr, numberComponents*lSD::pVolume(), MPI_DOUBLE,
                 globalPtr, numberComponents*lSD::pVolume(), MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    }

    HOST
    T reduce(T * localPtr) {
      T localSum = (T) 0;
      for(auto iZ = lSD::sStart()[d::Z]; iZ < lSD::sEnd()[d::Z]; ++iZ) {
        for(auto iY = lSD::sStart()[d::Y]; iY < lSD::sEnd()[d::Y]; ++iY) {
          for(auto iX = lSD::sStart()[d::X]; iX < lSD::sEnd()[d::X]; ++iX) {
            localSum
              += localPtr[lSD::getIndex(Position({iX, iY, iZ}))];
          }
        }
      }

      MPI_Barrier(MPI_COMM_WORLD);

      T globalSum = (T) 0;
      MPI_Reduce(&localSum, &globalSum, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);

      MPI_Barrier(MPI_COMM_WORLD);

      return globalSum;
    }

  };


  template<class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::SoA, PartitionningType::Generic,
                      Implementation::MPI, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::Generic, PartitionningType::Generic,
                           Implementation::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull,
                               MemoryLayout::Generic, PartitionningType::Generic,
                               Implementation::MPI, 0>;
    using Base::rightXRankMPI;
    using Base::leftXRankMPI;
    using Base::statusXMPI;
    using Base::requestXMPI;

    using Base::rightYRankMPI;
    using Base::leftYRankMPI;
    using Base::statusYMPI;
    using Base::requestYMPI;

    using Base::rightZRankMPI;
    using Base::leftZRankMPI;
    using Base::statusZMPI;
    using Base::requestZMPI;


    typedef Domain<DomainType::HaloSpace, PartitionningType::Generic,
                   MemoryLayout::SoA, L::dimQ> hMLSD;

    unsigned int sizeStripeX;
    unsigned int sendToRightBeginX;
    unsigned int receivedFromLeftBeginX;
    unsigned int sendToLeftBeginX;
    unsigned int receivedFromRightBeginX;

    unsigned int sizeStripeY;
    unsigned int sendToRightBeginY;
    unsigned int receivedFromLeftBeginY;
    unsigned int sendToLeftBeginY;
    unsigned int receivedFromRightBeginY;

    unsigned int sizeStripeZ;
    unsigned int sendToRightBeginZ;
    unsigned int receivedFromLeftBeginZ;
    unsigned int sendToLeftBeginZ;
    unsigned int receivedFromRightBeginZ;

  protected:
    HOST
    void sendAndReceiveHaloX(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloX",4)

      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        sendToRightBeginX = hMLSD::getIndex(Position({L::halo()[d::X]+lSD::sLength()[d::X]-1,
                hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);
        receivedFromLeftBeginX = hMLSD::getIndex(Position({0,
                hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

        sendToLeftBeginX = hMLSD::getIndex(Position({L::halo()[d::X],
                hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

        receivedFromRightBeginX = hMLSD::getIndex(Position({L::halo()[d::X]+lSD::sLength()[d::X],
                hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

        MPI_Irecv(haloDistributionPtr+receivedFromLeftBeginX, sizeStripeX,
                  MPI_DOUBLE, leftXRankMPI, 17, MPI_COMM_WORLD, &requestXMPI[0]);

        MPI_Irecv(haloDistributionPtr+receivedFromRightBeginX, sizeStripeX,
                  MPI_DOUBLE, rightXRankMPI, 23, MPI_COMM_WORLD, &requestXMPI[1]);

        MPI_Isend(haloDistributionPtr+sendToRightBeginX, sizeStripeX,
                  MPI_DOUBLE, rightXRankMPI, 17, MPI_COMM_WORLD, &requestXMPI[2]);

        MPI_Isend(haloDistributionPtr+sendToLeftBeginX, sizeStripeX,
                  MPI_DOUBLE, leftXRankMPI, 23, MPI_COMM_WORLD, &requestXMPI[3]);

        MPI_Waitall(4, requestXMPI, statusXMPI);
      }
    }

    HOST
    void sendAndReceiveHaloY(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloY",4)

        //TODO - PACK AND UNPACK


        // for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        //     sendToRightBeginY = hMLSD::getIndex(Position({hMLSD::start()[d::X],
        //             L::halo()[d::Y]+lSD::sLength()[d::Y]-1,
        //             hMLSD::start()[d::Z]}), iQ);
        //     receivedFromLeftBeginY = hMLSD::getIndex(Position({hMLSD::start()[d::X], 0, hMLSD::start()[d::Z]}), iQ);

        //     sendToLeftBeginY = hMLSD::getIndex(Position({hMLSD::start()[d::X],
        //             L::halo()[d::Y], hMLSD::start()[d::Z]}), iQ);

        //     receivedFromRightBeginY = hMLSD::getIndex(Position({hMLSD::start()[d::X], L::halo()[d::Y]+lSD::sLength()[d::Y], hMLSD::start()[d::Z]}), iQ);


        //     MPI_Irecv(haloDistributionPtr+receivedFromLeftBeginY, sizeStripeY,
        //               MPI_DOUBLE, leftYRankMPI, 17, MPI_COMM_WORLD, &requestYMPI[0]);

        //     MPI_Irecv(haloDistributionPtr+receivedFromRightBeginY, sizeStripeY,
        //               MPI_DOUBLE, rightYRankMPI, 23, MPI_COMM_WORLD, &requestYMPI[1]);

        //     MPI_Isend(haloDistributionPtr+sendToRightBeginY, sizeStripeY,
        //               MPI_DOUBLE, rightYRankMPI, 17, MPI_COMM_WORLD, &requestYMPI[2]);

        //     MPI_Isend(haloDistributionPtr+sendToLeftBeginY, sizeStripeY,
        //               MPI_DOUBLE, leftYRankMPI, 23, MPI_COMM_WORLD, &requestYMPI[3]);

        //     MPI_Waitall(4, requestYMPI, statusYMPI);
        //   }
        }

    HOST
    void sendAndReceiveHaloZ(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloZ",4)

        //TODO: PACK AND UNPACK

        // for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        //     sendToRightBeginZ = hMLSD::getIndex(Position({hMLSD::start()[d::X],
        //             hMLSD::start()[d::Y],
        //             L::halo()[d::Z]+lSD::sLength()[d::Z]-1}), iQ);

        //     receivedFromLeftBeginZ = hMLSD::getIndex(Position({hMLSD::start()[d::X],
        //             hMLSD::start()[d::Y], 0}), iQ);

        //     sendToLeftBeginZ = hMLSD::getIndex(Position({hMLSD::start()[d::X],
        //             hMLSD::start()[d::Y], L::halo()[d::Z]}), iQ);

        //     receivedFromRightBeginZ = hMLSD::getIndex(Position({hMLSD::start()[d::X], hMLSD::start()[d::Y], L::halo()[d::Z]+lSD::sLength()[d::Z]}), iQ);


        //     MPI_Irecv(haloDistributionPtr+receivedFromLeftBeginZ, sizeStripeZ,
        //               MPI_DOUBLE, leftZRankMPI, 17, MPI_COMM_WORLD, &requestZMPI[0]);

        //     MPI_Irecv(haloDistributionPtr+receivedFromRightBeginZ, sizeStripeZ,
        //               MPI_DOUBLE, rightZRankMPI, 23, MPI_COMM_WORLD, &requestZMPI[1]);

        //     MPI_Isend(haloDistributionPtr+sendToRightBeginZ, sizeStripeZ,
        //               MPI_DOUBLE, rightZRankMPI, 17, MPI_COMM_WORLD, &requestZMPI[2]);

        //     MPI_Isend(haloDistributionPtr+sendToLeftBeginZ, sizeStripeZ,
        //               MPI_DOUBLE, leftZRankMPI, 23, MPI_COMM_WORLD, &requestZMPI[3]);

        //     MPI_Waitall(4, requestZMPI, statusZMPI);
        //   }
        }


  public:
    HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::Generic, PartitionningType::Generic,
                      Implementation::MPI, 0>(rankMPI_in, sizeMPI_in,
                                         processorName_in)
      , sizeStripeX(hMLSD::volume()/hMLSD::length()[d::X])
      , sendToRightBeginX(0)
      , receivedFromLeftBeginX(0)
      , sendToLeftBeginX(0)
      , receivedFromRightBeginX(0)
      , sizeStripeY(hMLSD::volume()/hMLSD::length()[d::Y])
      , sendToRightBeginY(0)
      , receivedFromLeftBeginY(0)
      , sendToLeftBeginY(0)
      , receivedFromRightBeginY(0)
      , sizeStripeZ(hMLSD::volume()/hMLSD::length()[d::Z])
      , sendToRightBeginZ(0)
      , receivedFromLeftBeginZ(0)
      , sendToLeftBeginZ(0)
      , receivedFromRightBeginZ(0)

    {}

    using Base::printInputs;
    using Base::getRankMPI;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
    using Base::reduce;
  };


  template<class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::AoS, PartitionningType::Generic, Implementation::MPI, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::Generic, PartitionningType::Generic,
                           Implementation::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, PartitionningType::Generic,
                               Implementation::MPI, 0>;

    using Base::rightXRankMPI;
    using Base::leftXRankMPI;
    using Base::statusXMPI;
    using Base::requestXMPI;

    using Base::rightYRankMPI;
    using Base::leftYRankMPI;
    using Base::statusYMPI;
    using Base::requestYMPI;

    using Base::rightZRankMPI;
    using Base::leftZRankMPI;
    using Base::statusZMPI;
    using Base::requestZMPI;


    typedef Domain<DomainType::HaloSpace, PartitionningType::Generic,
                   MemoryLayout::AoS, L::dimQ> hMLSD;

    unsigned int sizeStripeX;
    unsigned int sendToRightBeginX;
    unsigned int receivedFromLeftBeginX;
    unsigned int sendToLeftBeginX;
    unsigned int receivedFromRightBeginX;

    unsigned int sizeStripeY;
    unsigned int sendToRightBeginY;
    unsigned int receivedFromLeftBeginY;
    unsigned int sendToLeftBeginY;
    unsigned int receivedFromRightBeginY;

    unsigned int sizeStripeZ;
    unsigned int sendToRightBeginZ;
    unsigned int receivedFromLeftBeginZ;
    unsigned int sendToLeftBeginZ;
    unsigned int receivedFromRightBeginZ;

  public:
    HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
      MemoryLayout::Generic, PartitionningType::Generic,
                        Implementation::MPI, 0>(rankMPI_in, sizeMPI_in,
                                         processorName_in)
      , sizeStripeX(L::dimQ*hMLSD::volume()*L::halo()[d::X]/hMLSD::length()[d::X])
      , sendToRightBeginX(hMLSD::getIndex(Position({L::halo()[d::X]+lSD::sLength()[d::X]-1,
                hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), 0))
      , receivedFromLeftBeginX(hMLSD::getIndex(Position({0,
                hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), 0))
      , sendToLeftBeginX(hMLSD::getIndex(Position({L::halo()[d::X],
                hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), 0))
      , receivedFromRightBeginX(hMLSD::getIndex(Position({L::halo()[d::X]+lSD::sLength()[d::X],
                hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), 0))
      , sizeStripeY(L::dimQ*hMLSD::volume()*L::halo()[d::Y]/hMLSD::length()[d::Y])

      , sendToRightBeginY(hMLSD::getIndex(Position({hMLSD::start()[d::X], L::halo()[d::Y]+lSD::sLength()[d::Y]-1, hMLSD::start()[d::Z]}), 0))

      , receivedFromLeftBeginY(hMLSD::getIndex(Position({hMLSD::start()[d::X],
                0, hMLSD::start()[d::Z]}), 0))
      , sendToLeftBeginY(hMLSD::getIndex(Position({hMLSD::start()[d::X],
                L::halo()[d::Y], hMLSD::start()[d::Z]}), 0))
      , receivedFromRightBeginY(hMLSD::getIndex(Position({hMLSD::start()[d::X], L::halo()[d::Y]+lSD::sLength()[d::Y], hMLSD::start()[d::Z]}), 0))


      , sizeStripeZ(L::dimQ*hMLSD::volume()*L::halo()[d::Z]/hMLSD::length()[d::Z])

      , sendToRightBeginZ(hMLSD::getIndex(Position({hMLSD::start()[d::X], hMLSD::start()[d::Y], L::halo()[d::Z]+lSD::sLength()[d::Z]-1}), 0))

      , receivedFromLeftBeginZ(hMLSD::getIndex(Position({hMLSD::start()[d::X],
                hMLSD::start()[d::Y], 0}), 0))
      , sendToLeftBeginZ(hMLSD::getIndex(Position({hMLSD::start()[d::X],
                hMLSD::start()[d::Y], L::halo()[d::Z]}), 0))
      , receivedFromRightBeginZ(hMLSD::getIndex(Position({hMLSD::start()[d::X], hMLSD::start()[d::Y], L::halo()[d::Z]+lSD::sLength()[d::Z]}), 0))

    {}

    using Base::printInputs;
    using Base::getRankMPI;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
    using Base::reduce;

  protected:
    HOST
    void sendAndReceiveHaloX(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloX",4)

      MPI_Irecv(haloDistributionPtr+receivedFromLeftBeginX, sizeStripeX,
                  MPI_DOUBLE, leftXRankMPI, 17, MPI_COMM_WORLD, &requestXMPI[0]);

      MPI_Irecv(haloDistributionPtr+receivedFromRightBeginX, sizeStripeX,
                MPI_DOUBLE, rightXRankMPI, 23, MPI_COMM_WORLD, &requestXMPI[1]);

      MPI_Isend(haloDistributionPtr+sendToRightBeginX, sizeStripeX,
                MPI_DOUBLE, rightXRankMPI, 17, MPI_COMM_WORLD, &requestXMPI[2]);

      MPI_Isend(haloDistributionPtr+sendToLeftBeginX, sizeStripeX,
                MPI_DOUBLE, leftXRankMPI, 23, MPI_COMM_WORLD, &requestXMPI[3]);

      MPI_Waitall(4, requestXMPI, statusXMPI);
    }

    HOST
    void sendAndReceiveHaloY(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloY",4)

      MPI_Irecv(haloDistributionPtr+receivedFromLeftBeginY, sizeStripeY,
                  MPI_DOUBLE, leftYRankMPI, 17, MPI_COMM_WORLD, &requestYMPI[0]);

      MPI_Irecv(haloDistributionPtr+receivedFromRightBeginY, sizeStripeY,
                MPI_DOUBLE, rightYRankMPI, 23, MPI_COMM_WORLD, &requestYMPI[1]);

      MPI_Isend(haloDistributionPtr+sendToRightBeginY, sizeStripeY,
                MPI_DOUBLE, rightYRankMPI, 17, MPI_COMM_WORLD, &requestYMPI[2]);

      MPI_Isend(haloDistributionPtr+sendToLeftBeginY, sizeStripeY,
                MPI_DOUBLE, leftYRankMPI, 23, MPI_COMM_WORLD, &requestYMPI[3]);

      MPI_Waitall(4, requestYMPI, statusYMPI);
    }

    HOST
    void sendAndReceiveHaloZ(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloZ",4)

        MPI_Irecv(haloDistributionPtr+receivedFromLeftBeginZ, sizeStripeZ,
                  MPI_DOUBLE, leftZRankMPI, 17, MPI_COMM_WORLD, &requestZMPI[0]);

      MPI_Irecv(haloDistributionPtr+receivedFromRightBeginZ, sizeStripeZ,
                MPI_DOUBLE, rightZRankMPI, 23, MPI_COMM_WORLD, &requestZMPI[1]);

      MPI_Isend(haloDistributionPtr+sendToRightBeginZ, sizeStripeZ,
                MPI_DOUBLE, rightZRankMPI, 17, MPI_COMM_WORLD, &requestZMPI[2]);

      MPI_Isend(haloDistributionPtr+sendToLeftBeginZ, sizeStripeZ,
                MPI_DOUBLE, leftZRankMPI, 23, MPI_COMM_WORLD, &requestZMPI[3]);

      MPI_Waitall(4, requestZMPI, statusZMPI);
    }
  };

#ifdef USE_NVSHMEM

  template<class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::SoA, PartitionningType::Generic,
                      Implementation::NVSHMEM_OUT, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::SoA, PartitionningType::Generic,
                           Implementation::MPI, 0> {
  private:
   using Base = Communication<T, latticeType, AlgorithmType::Pull,
                              MemoryLayout::SoA, PartitionningType::Generic,
                              Implementation::MPI, 0>;

    using Base::rightXRankMPI;
    using Base::leftXRankMPI;
    using Base::statusXMPI;
    using Base::requestXMPI;

    using Base::rightYRankMPI;
    using Base::leftYRankMPI;
    using Base::statusYMPI;
    using Base::requestYMPI;

    using Base::rightZRankMPI;
    using Base::leftZRankMPI;
    using Base::statusZMPI;
    using Base::requestZMPI;


    typedef Domain<DomainType::HaloSpace, PartitionningType::Generic,
                   MemoryLayout::SoA, L::dimQ> hMLSD;

    using Base::sizeStripeX;
    using Base::sendToRightBeginX;
    using Base::receivedFromLeftBeginX;
    using Base::sendToLeftBeginX;
    using Base::receivedFromRightBeginX;

    using Base::sizeStripeY;
    using Base::sendToRightBeginY;
    using Base::receivedFromLeftBeginY;
    using Base::sendToLeftBeginY;
    using Base::receivedFromRightBeginY;

    using Base::sizeStripeZ;
    using Base::sendToRightBeginZ;
    using Base::receivedFromLeftBeginZ;
    using Base::sendToLeftBeginZ;
    using Base::receivedFromRightBeginZ;

  protected:
    HOST
    void sendAndReceiveHaloX(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloX",4)
        //TODO: Replace with NVSHMEM directives
        for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
          sendToRightBeginX = hMLSD::getIndex(Position({L::halo()[d::X]+lSD::sLength()[d::X]-1,
                  hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);
          receivedFromLeftBeginX = hMLSD::getIndex(Position({0,
                  hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

          sendToLeftBeginX = hMLSD::getIndex(Position({L::halo()[d::X],
                  hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

          receivedFromRightBeginX = hMLSD::getIndex(Position({L::halo()[d::X]+lSD::sLength()[d::X],
                  hMLSD::start()[d::Y], hMLSD::start()[d::Z]}), iQ);

          /* MPI_Irecv(haloDistributionPtr+receivedFromLeftBeginX, sizeStripeX, */
          /*           MPI_DOUBLE, leftXRankMPI, 17, MPI_COMM_WORLD, &requestXMPI[0]); */

          /* MPI_Irecv(haloDistributionPtr+receivedFromRightBeginX, sizeStripeX, */
          /*           MPI_DOUBLE, rightXRankMPI, 23, MPI_COMM_WORLD, &requestXMPI[1]); */

          /* MPI_Isend(haloDistributionPtr+sendToRightBeginX, sizeStripeX, */
          /*           MPI_DOUBLE, rightXRankMPI, 17, MPI_COMM_WORLD, &requestXMPI[2]); */

          /* MPI_Isend(haloDistributionPtr+sendToLeftBeginX, sizeStripeX, */
          /*           MPI_DOUBLE, leftXRankMPI, 23, MPI_COMM_WORLD, &requestXMPI[3]); */

          /* MPI_Waitall(4, requestXMPI, statusXMPI); */

          shmem_double_put(haloDistributionPtr+receivedFromLeftBeginX,
                           haloDistributionPtr+sendToRightBeginX,
                           sizeStripeX, rightXRankMPI);
          shmem_double_put(haloDistributionPtr+sendToLeftBeginX,
                           haloDistributionPtr+receivedFromRightBeginX,
                           sizeStripeX, leftXRankMPI);
          shmem_barrier_all();

        }
    }

    using Base::sendAndReceiveHaloY;

    using Base::sendAndReceiveHaloZ;



  public:
    using Base::Communication;

    using Base::printInputs;
    using Base::getRankMPI;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
    using Base::reduce;
  };


  template<class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::AoS, PartitionningType::Generic, Implementation::NVSHMEM_OUT, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::AoS, PartitionningType::Generic,
                           Implementation::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::AoS, PartitionningType::Generic,
                               Implementation::MPI, 0>;

    using Base::rightXRankMPI;
    using Base::leftXRankMPI;
    using Base::statusXMPI;
    using Base::requestXMPI;

    using Base::rightYRankMPI;
    using Base::leftYRankMPI;
    using Base::statusYMPI;
    using Base::requestYMPI;

    using Base::rightZRankMPI;
    using Base::leftZRankMPI;
    using Base::statusZMPI;
    using Base::requestZMPI;


    typedef Domain<DomainType::HaloSpace, PartitionningType::Generic,
                   MemoryLayout::AoS, L::dimQ> hMLSD;

    using Base::sizeStripeX;
    using Base::sendToRightBeginX;
    using Base::receivedFromLeftBeginX;
    using Base::sendToLeftBeginX;
    using Base::receivedFromRightBeginX;


    using Base::sizeStripeY;
    using Base::sendToRightBeginY;
    using Base::receivedFromLeftBeginY;
    using Base::sendToLeftBeginY;
    using Base::receivedFromRightBeginY;


    using Base::sizeStripeZ;
    using Base::sendToRightBeginZ;
    using Base::receivedFromLeftBeginZ;
    using Base::sendToLeftBeginZ;
    using Base::receivedFromRightBeginZ;

  using Base::Communication;

    using Base::printInputs;
    using Base::getRankMPI;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
    using Base::reduce;

  protected:
    HOST
    void sendAndReceiveHaloX(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloX",4)
        //TODO: Replace with NVSHMEM directives
      /* MPI_Irecv(haloDistributionPtr+receivedFromLeftBeginX, sizeStripeX, */
      /*             MPI_DOUBLE, leftXRankMPI, 17, MPI_COMM_WORLD, &requestXMPI[0]); */

      /* MPI_Irecv(haloDistributionPtr+receivedFromRightBeginX, sizeStripeX, */
      /*           MPI_DOUBLE, rightXRankMPI, 23, MPI_COMM_WORLD, &requestXMPI[1]); */

      /* MPI_Isend(haloDistributionPtr+sendToRightBeginX, sizeStripeX, */
      /*           MPI_DOUBLE, rightXRankMPI, 17, MPI_COMM_WORLD, &requestXMPI[2]); */

      /* MPI_Isend(haloDistributionPtr+sendToLeftBeginX, sizeStripeX, */
      /*           MPI_DOUBLE, leftXRankMPI, 23, MPI_COMM_WORLD, &requestXMPI[3]); */

      /* MPI_Waitall(4, requestXMPI, statusXMPI); */

      shmem_double_put(haloDistributionPtr+receivedFromLeftBeginX,
                       haloDistributionPtr+sendToRightBeginX,
                       sizeStripeX, rightXRankMPI);
      shmem_double_put(haloDistributionPtr+sendToLeftBeginX,
                       haloDistributionPtr+receivedFromRightBeginX,
                       sizeStripeX, leftXRankMPI);
      shmem_barrier_all();

    }

    using Base::sendAndReceiveHaloY;

    using Base::sendAndReceiveHaloZ;

  };


  template<class T, LatticeType latticeType, MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::Generic, Implementation::NVSHMEM_IN, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::Generic, PartitionningType::Generic,
                           Implementation::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull,
                               MemoryLayout::Generic, PartitionningType::Generic,
                               Implementation::MPI, 0>;
  public:
    using Base::Communication;

    using Base::printInputs;
    using Base::getRankMPI;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
    using Base::reduce;

  protected:
    HOST
    void sendAndReceiveHaloX(T * haloDistributionPtr) {
    }

    HOST
    void sendAndReceiveHaloY(T * haloDistributionPtr) {
    }

    HOST
    void sendAndReceiveHaloZ(T * haloDistributionPtr) {
    }
  };
#endif // USE_NVSHMEM


  template<class T, LatticeType latticeType,
           MemoryLayout memoryLayout, unsigned int Dimension>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::OneD,  Implementation::MPI, Dimension>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic,
                           Implementation::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull,
                               memoryLayout, PartitionningType::Generic,
                               Implementation::MPI, 0>;

  public:
    HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::Generic,
                      Implementation::MPI, 0>(rankMPI_in, sizeMPI_in,
                                       processorName_in)
    {}

    using Base::printInputs;
    using Base::getRankMPI;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
    using Base::reduce;

    HOST
    inline void communicateHalos(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<6>::communicateHalos",3)
        sendAndReceiveHaloX(haloDistributionPtr);
    }

  private:
    using Base::sendAndReceiveHaloX;
  };

  template<class T, LatticeType latticeType, MemoryLayout memoryLayout,
           unsigned int Dimension>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::TwoD, Implementation::MPI, Dimension>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic, Implementation::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic,
                               Implementation::MPI, 0>;

  public:
    using Base::Communication;
    using Base::printInputs;
    using Base::getRankMPI;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
    using Base::reduce;

  private:
    using Base::sendAndReceiveHaloX;
  };

  template<class T, LatticeType latticeType, MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::TwoD, Implementation::MPI, 2>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic,
                           Implementation::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic,
                               Implementation::MPI, 0>;

  public:
    using Base::Communication;
    using Base::printInputs;
    using Base::getRankMPI;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
    using Base::reduce;

    HOST
    inline void communicateHalos(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<6>::communicateHalos",3)

      sendAndReceiveHaloX(haloDistributionPtr);
      sendAndReceiveHaloY(haloDistributionPtr);
    }

  private:
    using Base::sendAndReceiveHaloX;
    using Base::sendAndReceiveHaloY;
  };

  template<class T, LatticeType latticeType, MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::TwoD, Implementation::MPI, 3>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::TwoD, Implementation::MPI, 2> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull,
                               memoryLayout, PartitionningType::TwoD,
                               Implementation::MPI, 2>;

  public:
    using Base::Communication;
    using Base::printInputs;
    using Base::getRankMPI;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
    using Base::reduce;
    using Base::communicateHalos;

  private:
    using Base::sendAndReceiveHaloX;
    using Base::sendAndReceiveHaloY;
  };


  template<class T, LatticeType latticeType,
           MemoryLayout memoryLayout, unsigned int Dimension>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::ThreeD, Implementation::MPI, Dimension>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic, Implementation::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull,
                               memoryLayout, PartitionningType::Generic,
                               Implementation::MPI, 0>;

    using Base::Communication;

    using Base::printInputs;
    using Base::getRankMPI;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
    using Base::reduce;

  private:
    using Base::sendAndReceiveHaloX;
  };

  template<class T, LatticeType latticeType, MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::ThreeD, Implementation::MPI, 3>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic, Implementation::MPI, 0> {
  private:
    using Base = Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic,
                               Implementation::MPI, 0>;

  public:
    using Base::Communication;

    using Base::printInputs;
    using Base::getRankMPI;
    using Base::sendGlobalToLocal;
    using Base::sendLocalToGlobal;
    using Base::reduce;

    HOST
    inline void communicateHalos(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<6>::communicateHalos",3)

      sendAndReceiveHaloZ(haloDistributionPtr);
      sendAndReceiveHaloY(haloDistributionPtr);
      sendAndReceiveHaloX(haloDistributionPtr);
    }

  private:
    using Base::sendAndReceiveHaloX;
    using Base::sendAndReceiveHaloY;
    using Base::sendAndReceiveHaloZ;
  };

}

#endif // COMMUNICATION_H
