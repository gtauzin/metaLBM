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

  template<class T, PartitionningType partitionningType, unsigned int dimension>
  class Communication : public Boundary<T, BoundaryType::Periodic>
  {
  protected:
    const MathVector<unsigned int, 3> rankMPI;
    const MathVector<unsigned int, 3> sizeMPI;
    const std::string processorName;

    using Boundary<T, BoundaryType::Periodic>::applyX;
    using Boundary<T, BoundaryType::Periodic>::applyY;
    using Boundary<T, BoundaryType::Periodic>::applyZ;

  public:
    Communication(const MathVector<unsigned int, 3>& rankMPI_in,
                  const MathVector<unsigned int, 3>& sizeMPI_in,
                  const std::string& processorName_in)
      : rankMPI(rankMPI_in)
      , sizeMPI(sizeMPI_in)
      , processorName(processorName_in)
    {}

    void printInputs() {
      std::cout << "MPI #" << rankMPI << " of " << sizeMPI
                << " processes running on host " << processorName << "\n";
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

    T reduceLocal(T * __restrict__ localFieldArray) {
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

  template<class T, PartitionningType partitionningType>
  class Communication<T, partitionningType, 1>
    : public Communication<T, partitionningType, 0> {
  private:
    using Communication<T, partitionningType, 0>::rankMPI;
    using Communication<T, partitionningType, 0>::sizeMPI;
    using Communication<T, partitionningType, 0>::processorName;

    using Communication<T, partitionningType, 0>::applyX;
    using Communication<T, partitionningType, 0>::applyY;
    using Communication<T, partitionningType, 0>::applyZ;

  public:
    using Communication<T, partitionningType, 0>::Communication;

    using Communication<T, partitionningType, 0>::printInputs;
    using Communication<T, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, partitionningType, 0>::reduceLocal;

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

  template<class T, PartitionningType partitionningType>
  class Communication<T, partitionningType, 2>
    : public Communication<T, partitionningType, 0> {
  private:
    using Communication<T, partitionningType, 0>::rankMPI;
    using Communication<T, partitionningType, 0>::sizeMPI;
    using Communication<T, partitionningType, 0>::processorName;

    using Communication<T, partitionningType, 0>::applyX;
    using Communication<T, partitionningType, 0>::applyY;
    using Communication<T, partitionningType, 0>::applyZ;

  public:
    using Communication<T, partitionningType, 0>::Communication;

    using Communication<T, partitionningType, 0>::printInputs;
    using Communication<T, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, partitionningType, 0>::reduceLocal;

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

  template<class T, PartitionningType partitionningType>
  class Communication<T, partitionningType, 3>
    : public Communication<T, partitionningType, 0> {
  private:
    using Communication<T, partitionningType, 0>::rankMPI;
    using Communication<T, partitionningType, 0>::sizeMPI;
    using Communication<T, partitionningType, 0>::processorName;

    using Communication<T, partitionningType, 0>::applyX;
    using Communication<T, partitionningType, 0>::applyY;
    using Communication<T, partitionningType, 0>::applyZ;

  public:
    using Communication<T, partitionningType, 0>::Communication;

    using Communication<T, partitionningType, 0>::printInputs;
    using Communication<T, partitionningType, 0>::sendGlobalToLocal;
    using Communication<T, partitionningType, 0>::sendLocalToGlobal;
    using Communication<T, partitionningType, 0>::reduceLocal;

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
          std::cout << "iP in periodic Z: " << iP << std::endl;
          applyZ(f, iP);
        }
      }

    }

  };


  typedef Communication<dataT, partitionningT, L::dimD> Communication_;
}

#endif // COMMUNICATION_H
