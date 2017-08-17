#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <omp.h>
#include <string>

#include "Options.h"
#include "Lattice.h"
#include "Domain.h"
#include "Computation.h"
#include "Boundary.h"

namespace lbm {
  // TODO: 2D loop defined in Computation.h ? Maybe -- not for now???

  template<class T, LatticeType latticeType, AlgorithmType algorithmType,
           PartitionningType partitionningType, unsigned int Dimension>
  class Communication {};

  template<class T, LatticeType latticeType, PartitionningType partitionningType>
  class Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 0>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull> {
  protected:
    const MathVector<int, 3> rankMPI;
    const MathVector<int, 3> sizeMPI;
    const std::string processorName;

    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>()
      , rankMPI(rankMPI_in)
      , sizeMPI(sizeMPI_in)
      , processorName(processorName_in)
    {}

    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>::applyX;
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>::applyY;
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>::applyZ;

  public:
    void printInputs() {
      std::cout << "MPI #" << rankMPI << " of " << sizeMPI
                << " processes running on host " << processorName << std::endl;;
    }

    MathVector<int, 3> getRankMPI() {
      return rankMPI;
    }

    void sendGlobalToLocal(T * __restrict__ globalFieldArray,
                           T * __restrict__ localFieldArray,
                           unsigned int numberComponents) {
      memcpy(localFieldArray, globalFieldArray, numberComponents*lD::volume()*sizeof(T));
    }

    void sendLocalToGlobal(T * __restrict__ localFieldArray,
                           T * __restrict__ globalFieldArray,
                           unsigned int numberComponents) {
      memcpy(globalFieldArray, localFieldArray, numberComponents*gD::volume()*sizeof(T));
    }

    T reduce(T * __restrict__ localFieldArray) {
      T sum = (T) 0;
      for(unsigned int iZ = lD::start()[d::Z]; iZ < lD::end()[d::Z]; ++iZ) {
        for(unsigned int iY = lD::start()[d::Y]; iY < lD::end()[d::Y]; ++iY) {
          for(unsigned int iX = lD::start()[d::X]; iX < lD::end()[d::X]; ++iX) {
            sum += localFieldArray[lD::getIndex({iX, iY, iZ})];
          }
        }
      }
      return sum;
    }

  };

  template<class T, LatticeType latticeType, PartitionningType partitionningType>
  class Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 1>
    : public Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 0> {
  private:
    using Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 0>::applyX;

  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                            processorName_in)
    {}

    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::printInputs;
    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::getRankMPI;
    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::reduce;

    inline void periodic(T * __restrict__ f) {
      MathVector<unsigned int, 3> iP;

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for (unsigned int iZ = hD::start()[d::Z]; iZ < hD::end()[d::Z]; ++iZ) {
#pragma omp simd
        for (unsigned int iY = hD::start()[d::Y]; iY < hD::end()[d::Y]; ++iY) {
          iP = MathVector<unsigned int, 3>({0, iY, iZ});
          applyX(f, iP);
        }
      }

    }

  };

  template<class T, LatticeType latticeType, PartitionningType partitionningType>
  class Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 2>
    : public Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 0> {
  private:
    using Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 0>::applyX;
    using Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 0>::applyY;

  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                            processorName_in)
    {}

    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::printInputs;
    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::getRankMPI;
    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::reduce;

    inline void periodic(T * __restrict__ f) {
      MathVector<unsigned int, 3> iP;

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for (unsigned int iZ = hD::start()[d::Z]; iZ < hD::end()[d::Z]; ++iZ) {
#pragma omp simd
        for (unsigned int iY = hD::start()[d::Y]; iY < hD::end()[d::Y]; ++iY) {
          iP = MathVector<unsigned int, 3>({0, iY, iZ});
          applyX(f, iP);
        }
      }

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

  template<class T, LatticeType latticeType, PartitionningType partitionningType>
  class Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 3>
    : public Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 0> {
  private:
    using Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 0>::applyX;
    using Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 0>::applyY;
    using Communication<T, latticeType, AlgorithmType::Pull, partitionningType, 0>::applyZ;

  public:
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      partitionningType, 0>(rankMPI_in, sizeMPI_in,
                                            processorName_in)
    {}

    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::printInputs;
    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::getRankMPI;

    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, latticeType,
                        AlgorithmType::Pull, partitionningType, 0>::reduce;

    inline void periodic(T * __restrict__ f) {
      MathVector<unsigned int, 3> iP;

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for (unsigned int iZ = hD::start()[d::Z]; iZ < hD::end()[d::Z]; ++iZ) {
#pragma omp simd
        for (unsigned int iY = hD::start()[d::Y]; iY < hD::end()[d::Y]; ++iY) {
          iP = MathVector<unsigned int, 3>({0, iY, iZ});
          applyX(f, iP);
        }
      }

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


  typedef Communication<dataT, latticeT,
                        algorithmT, partitionningT, L::dimD> Communication_;
}

#endif // COMMUNICATION_H
