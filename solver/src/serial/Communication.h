#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <omp.h>

#include "Options.h"
#include "Computation.h"
#include "Boundary.h"
#include "Lattice.h"

namespace lbm {
  // call PBC in Boundary.h
  // loop defined in Computation.h ? Maybe not for now???

  template<class T, AlgorithmType algorithmType,
           PartitionningType partitionningType>
  class Communication {} : public Boundary<T, BoundaryType::Periodic>;

  template<class T, AlgorithmType algorithmType,
           PartitionningType partitionningType>
  class Communication {

    inline void pull_handleHalos(T * __restrict__ f) {
#pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for (int iZ = 0; iZ < 2*L::hZ + L::lZ_l; ++iZ) {
#pragma omp simd
        for (int iX = 0; iX < 2*L::hX + L::lX_l; ++iX) {
          MathVector<unsigned int, 3> iP = {iX, 0, iZ};

          apply<d::Y>(f, iP);

        }
      }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for (int iZ = 0; iZ < 2*L::hZ+L::lZ_l; ++iZ) {
#pragma omp simd
        for (int iY = 0; iY < 2*L::hY+L::lY_l; ++iY) {
          MathVector<unsigned int, 3> iP = {0, iY, iZ};

          apply<d::X>(f, iP);
        }
      }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for (int iY = 0; iY < 2*L::hY + L::lY_l; ++iY) {
#pragma omp simd
        for (int iX = 0; iX < 2*L::hX + L::lX_l; ++iX) {
          MathVector<unsigned int, 3> iP = {iZ, iY, 0};

          apply<d::Z>(f, iP);

        }
      }

    }
  };
}

#endif // COMMUNICATION_H
