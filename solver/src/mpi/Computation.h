#ifndef COMPUTATION_H
#define COMPUTATION_H

#include <omp.h>

#include "Options.h"
#include "Domain.h"
#include "MathVector.h"

namespace lbm {

  struct Computation {
  template<typename Callback>
    static void Do(Callback f) {
    MathVector<unsigned int, 3> iP;
      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for(unsigned int iZ = lD::start()[d::Z]+L::halo()[d::Z];
          iZ < lD::end()[d::Z]+L::halo()[d::Z]; ++iZ) {
        for(unsigned int iY = lD::start()[d::Y]+L::halo()[d::Y];
            iY < lD::end()[d::Y]+L::halo()[d::Y]; ++iY) {
          #pragma omp simd
          for(unsigned int iX = lD::start()[d::X]+L::halo()[d::X];
              iX < lD::end()[d::X]+L::halo()[d::X]; ++iX) {
            iP = {iX, iY, iZ};
            f(iP);
          }
        }
      }
    }
  };

}


#endif // COMPUTATION_H
