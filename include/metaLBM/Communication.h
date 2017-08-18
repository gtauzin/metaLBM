#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <mpi.h>
#include <omp.h>
#include <string>

#include "Options.h"
#include "Commons.h"
#include "Lattice.h"
#include "Domain.h"
#include "DynamicArray.h"
#include "Boundary.h"
#include "Computation.h"

namespace lbm {
  // TODO: 2D loop defined in Computation.h ? Maybe -- not for now???

  template<class T, LatticeType latticeType, AlgorithmType algorithmType,
           MemoryLayout memoryLayout, PartitionningType partitionningType,
           unsigned int Dimension>
  class Communication
  {};

  template<class T, LatticeType latticeType, PartitionningType partitionningType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::Generic, partitionningType, 0>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull> {
  protected:
    const MathVector<int, 3> rankMPI;
    const MathVector<int, 3> sizeMPI;
    const std::string processorName;

    const unsigned int rightXRankMPI;
    const unsigned int leftXRankMPI;
    MPI_Status statusMPI;

    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>()
      , rankMPI(rankMPI_in)
      , sizeMPI(sizeMPI_in)
      , processorName(processorName_in)
      , rightXRankMPI((rankMPI_in[d::X] + sizeMPI_in[d::X] - 1) % sizeMPI_in[d::X])
      , leftXRankMPI((rankMPI_in[d::X] + 1) % sizeMPI_in[d::X])
      , statusMPI()
    {}

    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>::applyY;
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>::applyZ;

    void printInputs() {
      std::cout << "MPI #" << rankMPI[d::X] << " of " << sizeMPI[d::X]
                << " running on host " << processorName << std::endl
                << "Right MPI #" << rightXRankMPI
                << ", left MPI #" << leftXRankMPI << std::endl;
    }

    MathVector<int, 3> getRankMPI() {
      return rankMPI;
    }

    void sendGlobalToLocal(T * __restrict__ globalFieldArray,
                           T * __restrict__ localFieldArray,
                           unsigned int numberComponents) {
      MPI_Scatter(globalFieldArray, numberComponents*lD::volume(), MPI_DOUBLE,
                  localFieldArray, numberComponents*lD::volume(), MPI_DOUBLE,
                  0, MPI_COMM_WORLD);
    }

    void sendLocalToGlobal(T * __restrict__ localFieldArray,
                           T * __restrict__ globalFieldArray,
                           unsigned int numberComponents) {
      MPI_Gather(localFieldArray, numberComponents*lD::volume(), MPI_DOUBLE,
                 globalFieldArray, numberComponents*lD::volume(), MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    }

    T reduce(T * __restrict__ localFieldArray) {
      T localSum = (T) 0;
      for(unsigned int iZ = lD::start()[d::Z]; iZ < lD::end()[d::Z]; ++iZ) {
        for(unsigned int iY = lD::start()[d::Y]; iY < lD::end()[d::Y]; ++iY) {
          for(unsigned int iX = lD::start()[d::X]; iX < lD::end()[d::X]; ++iX) {
            localSum += localFieldArray[lD::getIndex({iX, iY, iZ})];
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

  template<class T, LatticeType latticeType, PartitionningType partitionningType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::Default, partitionningType, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::Generic, partitionningType, 0> {
  protected:
    const MathVector<int, 3> rankMPI;
    const MathVector<int, 3> sizeMPI;
    const std::string processorName;

    void sendAndReceiveHaloX(T * __restrict__ haloFieldArray) {
      packToSendX(haloFieldArray);

      MPI_Sendrecv(sendRightBufferX.data(), L::dimQ*sendRightBufferX.getVolume(),
                   MPI_DOUBLE, rightXRankMPI, 17,
                   receivedLeftBufferX.data(), L::dimQ*receivedLeftBufferX.getVolume(),
                   MPI_DOUBLE, leftXRankMPI, 17,
                   MPI_COMM_WORLD, &statusMPI);
      MPI_Sendrecv(sendLeftBufferX.data(), L::dimQ*sendLeftBufferX.getVolume(),
                   MPI_DOUBLE, leftXRankMPI, 23,
                   receivedRightBufferX.data(), L::dimQ*receivedRightBufferX.getVolume(),
                   MPI_DOUBLE, rightXRankMPI, 23,
                   MPI_COMM_WORLD, &statusMPI);

      unpackReceivedX(haloFieldArray);
    }

    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>::applyY;
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>::applyZ;


  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::Generic, partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                                                   processorName_in)
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::reduce;

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::rightXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::leftXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::statusMPI;

    LocalizedField<T, L::dimQ> sendRightBufferX;
    LocalizedField<T, L::dimQ> sendLeftBufferX;
    LocalizedField<T, L::dimQ> receivedRightBufferX;
    LocalizedField<T, L::dimQ> receivedLeftBufferX;

    void packToSendX(T * __restrict__ haloFieldArray) {
      MathVector<unsigned int, 3> iP_Origin;
      MathVector<unsigned int, 3> iP_Destination;

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for (unsigned int iZ = hD::start()[d::Z]; iZ < hD::end()[d::Z]; ++iZ) {
#pragma omp simd
        for (unsigned int iY = hD::start()[d::Y]; iY < hD::end()[d::Y]; ++iY) {
          iP_Destination = {0, iY, iZ};

          iP_Origin = {L::halo()[d::X], iY, iZ};
          UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
              sendLeftBufferX[bXD::getIndex(iP_Destination, iQ)]
                = haloFieldArray[hD::getIndex(iP_Origin, iQ)];
            });

          iP_Origin = {L::halo()[d::X]+lD::length()[d::X]-1, iY, iZ};
          UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
              sendRightBufferX[bXD::getIndex(iP_Destination, iQ)]
                = haloFieldArray[hD::getIndex(iP_Origin, iQ)];
            });
        }
      }
    }

    void unpackReceivedX(T * __restrict__ haloFieldArray) {
      MathVector<unsigned int, 3> iP_Origin;
      MathVector<unsigned int, 3> iP_Destination;

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for (unsigned int iZ = hD::start()[d::Z]; iZ < hD::end()[d::Z]; ++iZ) {
#pragma omp simd
        for (unsigned int iY = hD::start()[d::Y]; iY < hD::end()[d::Y]; ++iY) {
          iP_Origin = {0, iY, iZ};

          iP_Destination = {L::halo()[d::X]+lD::length()[d::X], iY, iZ};

          UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
              haloFieldArray[hD::getIndex(iP_Destination, iQ)]
                = receivedRightBufferX[bXD::getIndex(iP_Origin, iQ)];
            });

          iP_Destination = {0, iY, iZ};

          UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
              haloFieldArray[hD::getIndex(iP_Destination, iQ)]
                = receivedLeftBufferX[bXD::getIndex(iP_Origin, iQ)];
            });
        }
      }
    }
  };


  template<class T, LatticeType latticeType, PartitionningType partitionningType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::AoS, partitionningType, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::Generic, partitionningType, 0> {
  protected:
    void sendAndReceiveHaloX(T * __restrict__ haloFieldArray) {
      MPI_Sendrecv(haloFieldArray+sendToRightBeginX, sizeStripeX,
                   MPI_DOUBLE, rightXRankMPI, 17,
                   haloFieldArray+receivedFromLeftBeginX, sizeStripeX,
                   MPI_DOUBLE, leftXRankMPI, 17,
                   MPI_COMM_WORLD, &statusMPI);

      MPI_Sendrecv(haloFieldArray+sendToLeftBeginX, sizeStripeX,
                   MPI_DOUBLE, leftXRankMPI, 23,
                   haloFieldArray+receivedFromRightBeginX, sizeStripeX,
                   MPI_DOUBLE, rightXRankMPI, 23,
                   MPI_COMM_WORLD, &statusMPI);
    }

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::applyY;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::applyZ;


  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::Generic, partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                                                   processorName_in)
      , sizeStripeX(L::dimQ*hMLD::volume()*L::halo()[d::X]/hMLD::length()[d::X])
      , sendToRightBeginX(hMLD::getIndex({L::halo()[d::X]+lD::length()[d::X]-1,
              hMLD::start()[d::Y], hMLD::start()[d::Z]}, 0))
      , receivedFromLeftBeginX(hMLD::getIndex({0,
              hMLD::start()[d::Y], hMLD::start()[d::Z]}, 0))
      , sendToLeftBeginX(hMLD::getIndex({L::halo()[d::X],
              hMLD::start()[d::Y], hMLD::start()[d::Z]}, 0))
      , receivedFromRightBeginX(hMLD::getIndex({L::halo()[d::X]+lD::length()[d::X],
              hMLD::start()[d::Y], hMLD::start()[d::Z]}, 0))
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::reduce;
  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::rightXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::leftXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::statusMPI;

    typedef Domain<DomainType::Halo, PartitionningType::Generic,
                   MemoryLayout::AoS, L::dimQ> hMLD;

    unsigned int sizeStripeX;
    unsigned int sendToRightBeginX;
    unsigned int receivedFromLeftBeginX;
    unsigned int sendToLeftBeginX;
    unsigned int receivedFromRightBeginX;
  };


  template<class T, LatticeType latticeType, PartitionningType partitionningType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::SoA, partitionningType, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::Generic, partitionningType, 0> {
  protected:
    void sendAndReceiveHaloX(T * __restrict__ haloFieldArray) {
      UnrolledFor<0, L::dimQ>::Do([&] (unsigned int iQ) {
          sendToRightBeginX = hMLD::getIndex({L::halo()[d::X]+lD::length()[d::X]-1,
            hMLD::start()[d::Y], hMLD::start()[d::Z]}, iQ);
          receivedFromLeftBeginX = hMLD::getIndex({0,
            hMLD::start()[d::Y], hMLD::start()[d::Z]}, iQ);

          MPI_Sendrecv(haloFieldArray+sendToRightBeginX, sizeStripeX,
                       MPI_DOUBLE, rightXRankMPI, 17,
                       haloFieldArray+receivedFromLeftBeginX, sizeStripeX,
                       MPI_DOUBLE, leftXRankMPI, 17,
                       MPI_COMM_WORLD, &statusMPI);

          sendToLeftBeginX = hMLD::getIndex({L::halo()[d::X],
            hMLD::start()[d::Y], hMLD::start()[d::Z]}, iQ);

          receivedFromRightBeginX = hMLD::getIndex({L::halo()[d::X]+lD::length()[d::X],
            hMLD::start()[d::Y], hMLD::start()[d::Z]}, iQ);
          MPI_Sendrecv(haloFieldArray+sendToLeftBeginX, sizeStripeX,
                       MPI_DOUBLE, leftXRankMPI, 23,
                       haloFieldArray+receivedFromRightBeginX, sizeStripeX,
                       MPI_DOUBLE, rightXRankMPI, 23,
                       MPI_COMM_WORLD, &statusMPI);
        });
    }

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::applyY;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::applyZ;

  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::Generic, partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                                                   processorName_in)
      , sizeStripeX(hMLD::volume()/hMLD::length()[d::X])
      , sendToRightBeginX(0)
      , receivedFromLeftBeginX(0)
      , sendToLeftBeginX(0)
      , receivedFromRightBeginX(0)

    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::reduce;

    inline void periodic(T * __restrict__ f) {
      sendAndReceiveHaloX(f);
    }

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::rightXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::leftXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, partitionningType, 0>::statusMPI;

    typedef Domain<DomainType::Halo, PartitionningType::Generic,
                   MemoryLayout::SoA, L::dimQ> hMLD;

    unsigned int sizeStripeX;
    unsigned int sendToRightBeginX;
    unsigned int receivedFromLeftBeginX;
    unsigned int sendToLeftBeginX;
    unsigned int receivedFromRightBeginX;
  };


  template<class T, LatticeType latticeType,
           MemoryLayout memoryLayout, PartitionningType partitionningType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, partitionningType, 1>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, partitionningType, 0> {
  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                                          processorName_in)
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::reduce;

    inline void periodic(T * __restrict__ f) {
      sendAndReceiveHaloX(f);
    }


  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::sendAndReceiveHaloX;


  };

  template<class T, LatticeType latticeType,
           MemoryLayout memoryLayout, PartitionningType partitionningType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, partitionningType, 2>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, partitionningType, 0> {
  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                                          processorName_in)
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::reduce;

    inline void periodic(T * __restrict__ f) {
      sendAndReceiveHaloX(f);

      MathVector<unsigned int, 3> iP;

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for (unsigned int iZ = hD::start()[d::Z]; iZ < hD::end()[d::Z]; ++iZ) {
#pragma omp simd
        for (unsigned int iX = hD::start()[d::X]; iX < hD::end()[d::X]; ++iX) {
          iP = MathVector<unsigned int, 3>({iX, 0, iZ});
          applyY(f, iP);
        }
      }

    }

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::sendAndReceiveHaloX;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::applyY;
  };

  template<class T, LatticeType latticeType,
           MemoryLayout memoryLayout, PartitionningType partitionningType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, partitionningType, 3>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, partitionningType, 0> {
  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                                          processorName_in)
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::getRankMP;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::reduce;

    inline void periodic(T * __restrict__ f) {
      sendAndReceiveHaloX(f);

      MathVector<unsigned int, 3> iP;

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for (unsigned int iZ = hD::start()[d::Z]; iZ < hD::end()[d::Z]; ++iZ) {
#pragma omp simd
        for (unsigned int iX = hD::start()[d::X]; iX < hD::end()[d::X]; ++iX) {
          iP = MathVector<unsigned int, 3>({iX, 0, iZ});
          applyY(f, iP);
        }
      }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for (unsigned int iY = hD::start()[d::Y]; iY < hD::end()[d::Y]; ++iY) {
#pragma omp simd
        for (unsigned int iX = hD::start()[d::X]; iX < hD::end()[d::X]; ++iX) {
          iP = MathVector<unsigned int, 3>({iX, iY, 0});
          applyZ(f, iP);
        }
      }

    }


  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::sendAndReceiveHaloX;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::applyY;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, partitionningType, 0>::applyZ;
  };

  typedef Communication<dataT, latticeT, algorithmT,
                        memoryL, partitionningT, L::dimD> Communication_;
}

#endif // COMMUNICATION_H
