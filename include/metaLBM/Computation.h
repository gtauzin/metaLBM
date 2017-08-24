#ifndef COMPUTATION_H
#define COMPUTATION_H

#include <omp.h>

#include "Options.h"
#include "Domain.h"

namespace lbm {

  struct Computation {
  template<typename Callback>
  static void Do(MathVector<unsigned int, 3> start,
                 MathVector<unsigned int, 3> end,
                 Callback f) {
    SCOREP_INSTRUMENT_ON("Computation<Callback>::Do")

    MathVector<unsigned int, 3> iP;
      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(unsigned int iZ = start[d::Z]; iZ < end[d::Z]; ++iZ) {
      for(unsigned int iY = start[d::Y]; iY < end[d::Y]; ++iY) {
          #pragma omp simd
        for(unsigned int iX = start[d::X]; iX < end[d::X]; ++iX) {
            iP = {iX, iY, iZ};
            f(iP);
          }
        }
      }
    }
  };

}


#endif // COMPUTATION_H
