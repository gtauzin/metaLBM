#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <mpi.h>
#include <string>

#include "Commons.h"
#include "Options.h"
#include "DynamicArray.cuh"
#include "Lattice.h"
#include "Domain.h"
#include "Boundary.h"
#include "Computation.h"

namespace lbm {

  template<class T, LatticeType latticeType, AlgorithmType algorithmType,
           MemoryLayout memoryLayout, PartitionningType partitionningType,
           unsigned int Dimension>
  class Communication {};


  template<class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::Generic, PartitionningType::Generic, 0> {
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
      , rightXRankMPI((rankMPI_in[d::X] + sizeMPI_in[d::X] - 1) % sizeMPI_in[d::X])
      , leftXRankMPI((rankMPI_in[d::X] + 1) % sizeMPI_in[d::X])
      , statusXMPI()
      , requestXMPI()
      , rightYRankMPI((rankMPI_in[d::Y] + sizeMPI_in[d::Y] - 1) % sizeMPI_in[d::Y])
      , leftYRankMPI((rankMPI_in[d::Y] + 1) % sizeMPI_in[d::Y])
      , statusYMPI()
      , requestYMPI()
      , rightZRankMPI((rankMPI_in[d::Z] + sizeMPI_in[d::Z] - 1) % sizeMPI_in[d::Z])
      , leftZRankMPI((rankMPI_in[d::Z] + 1) % sizeMPI_in[d::Z])
      , statusZMPI()
      , requestZMPI()
    {}

    void printInputs() {
      std::cout << "MPI #" << rankMPI << " of " << sizeMPI
                << " running on host " << processorName << std::endl
                << "Right MPI #" << rightYRankMPI
                << ", left MPI #" << leftYRankMPI << std::endl;
    }

    DEVICE HOST
    MathVector<int, 3> getRankMPI() {
      return rankMPI;
    }

    HOST
    void sendGlobalToLocal(DynamicArray<T, Architecture::CPU>& globalArray,
                           DynamicArray<T, Architecture::CPU>& localHostArray,
                           DynamicArray<T, Architecture::GPU>& localDeviceArray,
                           unsigned int numberComponents) {
      INSTRUMENT_ON("Communication<6>::sendGlobalToLocal",3)

      MPI_Scatter(globalArray.data(), numberComponents*lD::volume(), MPI_DOUBLE,
                  localHostArray.data(), numberComponents*lD::volume(), MPI_DOUBLE,
                  0, MPI_COMM_WORLD);

      localDeviceArray.copyFrom(localHostArray);
    }

    HOST
    void sendLocalToGlobal(DynamicArray<T, Architecture::GPU>& localDeviceArray,
                           DynamicArray<T, Architecture::CPU>& localHostArray,
                           DynamicArray<T, Architecture::CPU>& globalArray,
                           unsigned int numberComponents) {
      INSTRUMENT_ON("Communication<6>::sendLocalToGlobal",3)

      localDeviceArray.copyTo(localHostArray);

      MPI_Gather(localHostArray.data(), numberComponents*lD::volume(), MPI_DOUBLE,
                 globalArray.data(), numberComponents*lD::volume(), MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    }

    HOST
    T reduce(const DynamicArray<T, Architecture::GPU>& localDeviceArray,
             DynamicArray<T, Architecture::CPU>& localHostArray) {

      localDeviceArray.copyTo(localHostArray);

      T localSum = (T) 0;
      for(unsigned int iZ = lD::start()[d::Z]; iZ < lD::end()[d::Z]; ++iZ) {
        for(unsigned int iY = lD::start()[d::Y]; iY < lD::end()[d::Y]; ++iY) {
          for(unsigned int iX = lD::start()[d::X]; iX < lD::end()[d::X]; ++iX) {
            localSum
              += localHostArray[lD::getIndex(MathVector<unsigned int, 3>({iX, iY, iZ}))];
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



  template<class T, LatticeType latticeType, MemoryLayout memoryLayout,
           unsigned int Dimension>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::Serial, Dimension> {
  public:
    HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
    {}

    void printInputs() {
      std::cout << "Serial run" << std::endl;
    }

    DEVICE HOST
    MathVector<int, 3> getRankMPI() {
      return MathVector<int, 3>{{0}};
    }

    HOST
    void sendGlobalToLocal(DynamicArray<T, Architecture::CPU>& globalArray,
                           DynamicArray<T, Architecture::CPU>& localHostArray,
                           DynamicArray<T, Architecture::GPU>& localDeviceArray,
                           unsigned int numberComponents) {
      INSTRUMENT_ON("Communication<6>::sendGlobalToLocal",3)

      localHostArray.copyFrom(globalArray);
      localDeviceArray.copyFrom(localHostArray);
    }

    HOST
    void sendLocalToGlobal(DynamicArray<T, Architecture::GPU>& localDeviceArray,
                           DynamicArray<T, Architecture::CPU>& localHostArray,
                           DynamicArray<T, Architecture::CPU>& globalArray,
                           unsigned int numberComponents) {
      INSTRUMENT_ON("Communication<6>::sendLocalToGlobal",3)

      localDeviceArray.copyTo(localHostArray);
      localHostArray.copyTo(globalArray);
    }

    HOST
    inline void communicateHalos(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<6>::communicateHalos",3)
    }

    HOST
    T reduce(const DynamicArray<T, Architecture::GPU>& localDeviceArray,
             DynamicArray<T, Architecture::CPU>& localHostArray) {

      localDeviceArray.copyTo(localHostArray);

      T localSum = (T) 0;
      for(unsigned int iZ = lD::start()[d::Z]; iZ < lD::end()[d::Z]; ++iZ) {
        for(unsigned int iY = lD::start()[d::Y]; iY < lD::end()[d::Y]; ++iY) {
          for(unsigned int iX = lD::start()[d::X]; iX < lD::end()[d::X]; ++iX) {
            localSum
              += localHostArray[lD::getIndex(MathVector<unsigned int, 3>({iX, iY, iZ}))];
          }
        }
      }

      return localSum;
    }
  };

  template<class T, LatticeType latticeType>
    class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::SoA, PartitionningType::Generic, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::Generic, PartitionningType::Generic, 0> {
      private:
  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::rightXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::leftXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::statusXMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::requestXMPI;

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::rightYRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::leftYRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::statusYMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::requestYMPI;

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::rightZRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::leftZRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::statusZMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::requestZMPI;


    typedef Domain<DomainType::Halo, PartitionningType::Generic,
                   MemoryLayout::SoA, L::dimQ> hMLD;

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
    void sendAndReceiveHaloX(T * RESTRICT haloDistributionPtr) {
      INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloX",4)

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
          sendToRightBeginX = hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X]+lD::length()[d::X]-1,
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), iQ);
          receivedFromLeftBeginX = hMLD::getIndex(MathVector<unsigned int, 3>({0,
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), iQ);

          sendToLeftBeginX = hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X],
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), iQ);

          receivedFromRightBeginX = hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X]+lD::length()[d::X],
                  hMLD::start()[d::Y], hMLD::start()[d::Z]}), iQ);

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
    void sendAndReceiveHaloY(T * RESTRICT haloDistributionPtr) {
      INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloY",4)

        //TODO - PACK AND UNPACK


      // for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
      //     sendToRightBeginY = hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X],
      //             L::halo()[d::Y]+lD::length()[d::Y]-1,
      //             hMLD::start()[d::Z]}), iQ);
      //     receivedFromLeftBeginY = hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X], 0, hMLD::start()[d::Z]}), iQ);

      //     sendToLeftBeginY = hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X],
      //             L::halo()[d::Y], hMLD::start()[d::Z]}), iQ);

      //     receivedFromRightBeginY = hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X], L::halo()[d::Y]+lD::length()[d::Y], hMLD::start()[d::Z]}), iQ);


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
    void sendAndReceiveHaloZ(T * RESTRICT haloDistributionPtr) {
      INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloZ",4)

        //TODO: PACK AND UNPACK

      // for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
      //     sendToRightBeginZ = hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X],
      //             hMLD::start()[d::Y],
      //             L::halo()[d::Z]+lD::length()[d::Z]-1}), iQ);

      //     receivedFromLeftBeginZ = hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X],
      //             hMLD::start()[d::Y], 0}), iQ);

      //     sendToLeftBeginZ = hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X],
      //             hMLD::start()[d::Y], L::halo()[d::Z]}), iQ);

      //     receivedFromRightBeginZ = hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X], hMLD::start()[d::Y], L::halo()[d::Z]+lD::length()[d::Z]}), iQ);


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
                      MemoryLayout::Generic,
                      PartitionningType::Generic, 0>(rankMPI_in, sizeMPI_in,
                                                     processorName_in)
      , sizeStripeX(hMLD::volume()/hMLD::length()[d::X])
      , sendToRightBeginX(0)
      , receivedFromLeftBeginX(0)
      , sendToLeftBeginX(0)
      , receivedFromRightBeginX(0)
      , sizeStripeY(hMLD::volume()/hMLD::length()[d::Y])
      , sendToRightBeginY(0)
      , receivedFromLeftBeginY(0)
      , sendToLeftBeginY(0)
      , receivedFromRightBeginY(0)
      , sizeStripeZ(hMLD::volume()/hMLD::length()[d::Z])
      , sendToRightBeginZ(0)
      , receivedFromLeftBeginZ(0)
      , sendToLeftBeginZ(0)
      , receivedFromRightBeginZ(0)

    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::reduce;
  };


  template<class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::AoS, PartitionningType::Generic, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::Generic, PartitionningType::Generic, 0> {
  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::rightXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::leftXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::statusXMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::requestXMPI;

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::rightYRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::leftYRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::statusYMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::requestYMPI;

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::rightZRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::leftZRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::statusZMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::requestZMPI;


    typedef Domain<DomainType::Halo, PartitionningType::Generic,
                   MemoryLayout::AoS, L::dimQ> hMLD;

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
                      MemoryLayout::Generic,
                      PartitionningType::Generic, 0>(rankMPI_in, sizeMPI_in,
                                                     processorName_in)
      , sizeStripeX(L::dimQ*hMLD::volume()*L::halo()[d::X]/hMLD::length()[d::X])
      , sendToRightBeginX(hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X]+lD::length()[d::X]-1,
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), 0))
      , receivedFromLeftBeginX(hMLD::getIndex(MathVector<unsigned int, 3>({0,
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), 0))
      , sendToLeftBeginX(hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X],
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), 0))
      , receivedFromRightBeginX(hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X]+lD::length()[d::X],
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), 0))
      , sizeStripeY(L::dimQ*hMLD::volume()*L::halo()[d::Y]/hMLD::length()[d::Y])

      , sendToRightBeginY(hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X], L::halo()[d::Y]+lD::length()[d::Y]-1, hMLD::start()[d::Z]}), 0))

      , receivedFromLeftBeginY(hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X],
                0, hMLD::start()[d::Z]}), 0))
      , sendToLeftBeginY(hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X],
                L::halo()[d::Y], hMLD::start()[d::Z]}), 0))
      , receivedFromRightBeginY(hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X], L::halo()[d::Y]+lD::length()[d::Y], hMLD::start()[d::Z]}), 0))


      , sizeStripeZ(L::dimQ*hMLD::volume()*L::halo()[d::Z]/hMLD::length()[d::Z])

      , sendToRightBeginZ(hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X], hMLD::start()[d::Y], L::halo()[d::Z]+lD::length()[d::Z]-1}), 0))

      , receivedFromLeftBeginZ(hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X],
                hMLD::start()[d::Y], 0}), 0))
      , sendToLeftBeginZ(hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X],
                hMLD::start()[d::Y], L::halo()[d::Z]}), 0))
      , receivedFromRightBeginZ(hMLD::getIndex(MathVector<unsigned int, 3>({hMLD::start()[d::X], hMLD::start()[d::Y], L::halo()[d::Z]+lD::length()[d::Z]}), 0))

    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::reduce;

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


  template<class T, LatticeType latticeType,
           MemoryLayout memoryLayout, unsigned int Dimension>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::OneD, Dimension>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic, 0> {
  public:
    HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::Generic, 0>(rankMPI_in, sizeMPI_in,
                                                                processorName_in)
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::reduce;

    HOST
    inline void communicateHalos(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<6>::communicateHalos",3)
      sendAndReceiveHaloX(haloDistributionPtr);
    }

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendAndReceiveHaloX;
  };

  template<class T, LatticeType latticeType,
           MemoryLayout memoryLayout, unsigned int Dimension>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::TwoD, Dimension>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic, 0> {
  public:
    HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::Generic, 0>(rankMPI_in, sizeMPI_in,
                                                                processorName_in)
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::reduce;

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendAndReceiveHaloX;
  };

  template<class T, LatticeType latticeType, MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::TwoD, 2>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic, 0> {
  public:
    HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::Generic, 0>(rankMPI_in, sizeMPI_in,
                                                                processorName_in)
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::reduce;

    HOST
    inline void communicateHalos(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<6>::communicateHalos",3)

      sendAndReceiveHaloX(haloDistributionPtr);
      sendAndReceiveHaloY(haloDistributionPtr);
    }

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendAndReceiveHaloX;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendAndReceiveHaloY;
  };

  template<class T, LatticeType latticeType, MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::TwoD, 3>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::TwoD, 2> {
  public:
    HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::TwoD, 2>(rankMPI_in, sizeMPI_in,
                                                                processorName_in)
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::TwoD, 2>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::TwoD, 2>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::TwoD, 2>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::TwoD, 2>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::TwoD, 2>::reduce;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::TwoD, 2>::communicateHalos;

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::TwoD, 2>::sendAndReceiveHaloX;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::TwoD, 2>::sendAndReceiveHaloY;
  };


  template<class T, LatticeType latticeType,
           MemoryLayout memoryLayout, unsigned int Dimension>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::ThreeD, Dimension>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic, 0> {
  public:
    HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::Generic, 0>(rankMPI_in, sizeMPI_in,
                                                                processorName_in)
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::reduce;

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendAndReceiveHaloX;
  };

  template<class T, LatticeType latticeType, MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::ThreeD, 3>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::Generic, 0> {
  public:
    HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::Generic, 0>(rankMPI_in, sizeMPI_in,
                                                                processorName_in)
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::reduce;

    HOST
    inline void communicateHalos(T * haloDistributionPtr) {
      INSTRUMENT_ON("Communication<6>::communicateHalos",3)

      sendAndReceiveHaloZ(haloDistributionPtr);
      sendAndReceiveHaloY(haloDistributionPtr);
      sendAndReceiveHaloX(haloDistributionPtr);
    }

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendAndReceiveHaloX;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendAndReceiveHaloY;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::Generic, 0>::sendAndReceiveHaloZ;
  };

}

#endif // COMMUNICATION_H
