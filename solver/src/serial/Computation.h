#ifndef COMPUTATION_H
#define COMPUTATION_H

#include <omp.h>

#include "Options.h"
#include "Domain.h"
#include "MathVector.h"

namespace lbm {

  // for 1 loop 2 loops and 3 loops possibilities

  template<DomainType domainType>
  struct Computation {
  private:
    MathVector<int, 3> iP;

    template<typename F>
    static void Do (F f, const int iteration) {
     #pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for(int iZ = D::start[d::Z]; iZ < D::end[d::Z]-1; ++iZ) {
        for(int iY = D::start[d::Y]; iY < D::end[d::Y]-1; ++iY) {
         #pragma omp simd
          for(int iX = D::start[d::X]; iX < D::end[d::X]-1; ++iX) {
            BOOST_LOG_TRIVIAL(debug) << " - (x, y, z) = "
                                     << "(" << iX
                                     << ", " << iY
                                     << ", " << iZ << ")";
            iP = {iX, iY, iZ};
            F(iP, iteration);
          }
        }
      }
    }
  };


}


#endif // COMPUTATION_H
