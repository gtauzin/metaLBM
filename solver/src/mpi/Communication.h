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
#include "Computation.h"
#include "Boundary.h"

namespace lbm {
  // TODO: 2D loop defined in Computation.h ? Maybe -- not for now???

  template<class T, AlgorithmType algorithmType,
           PartitionningType partitionningType, unsigned int Dimension>
  class Communication
  {};

  template<class T, PartitionningType partitionningType>
  class Communication<T, AlgorithmType::Pull, partitionningType, 0>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull> {
  protected:
    const MathVector<int, 3> rankMPI;
    const MathVector<int, 3> sizeMPI;
    const std::string processorName;

    const unsigned int rightXRankMPI;
    const unsigned int leftXRankMPI;
    MPI_Status statusMPI;

    LocalizedField<T, L::dimQ> sendRightBufferX;
    LocalizedField<T, L::dimQ> sendLeftBufferX;
    LocalizedField<T, L::dimQ> receivedRightBufferX;
    LocalizedField<T, L::dimQ> receivedLeftBufferX;

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
      , sendRightBufferX("sendRightBufferX", bXD::volume())
      , sendLeftBufferX("sendLeftBufferX", bXD::volume())
      , receivedRightBufferX("receivedRightBufferX", bXD::volume())
      , receivedLeftBufferX("receivedLeftBufferX", bXD::volume())
    {}

    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>::applyY;
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>::applyZ;

  public:
    void printInputs() {
      std::cout << "MPI #" << rankMPI[d::X] << " of " << sizeMPI[d::X]
                << " running on host " << processorName << std::endl
                << "Right MPI #" << rightXRankMPI
                << ", left MPI #" << leftXRankMPI << std::endl;
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


  template<class T, PartitionningType partitionningType>
  class Communication<T, AlgorithmType::Pull, partitionningType, 1>
    : public Communication<T, AlgorithmType::Pull, partitionningType, 0> {
  private:
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::sendAndReceiveHaloX;

  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, AlgorithmType::Pull,
                      partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                            processorName_in)
    {}

    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::printInputs;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::reduce;

    inline void periodic(T * __restrict__ f) {
      sendAndReceiveHaloX(f);
    }

  };

  template<class T, PartitionningType partitionningType>
  class Communication<T, AlgorithmType::Pull, partitionningType, 2>
    : public Communication<T, AlgorithmType::Pull, partitionningType, 0> {
  private:
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::sendAndReceiveHaloX;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::applyY;

  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, AlgorithmType::Pull,
                      partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                            processorName_in)
    {}

    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::printInputs;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::reduce;

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

  };

  template<class T, PartitionningType partitionningType>
  class Communication<T, AlgorithmType::Pull, partitionningType, 3>
    : public Communication<T, AlgorithmType::Pull, partitionningType, 0> {
  private:
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::sendAndReceiveHaloX;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::applyY;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::applyZ;

  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, AlgorithmType::Pull,
                      partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                            processorName_in)
    {}

    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::printInputs;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, AlgorithmType::Pull, partitionningType, 0>::reduce;

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

  };

  typedef Communication<dataT, algorithmT, partitionningT, L::dimD> Communication_;
}

#endif // COMMUNICATION_H
