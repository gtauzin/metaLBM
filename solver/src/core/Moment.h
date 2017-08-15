#ifndef MOMENT_H
#define MOMENT_H

#include <omp.h>

#include "Options.h"
#include "Commons.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Helpers.h"

namespace lbm {

  template <class T>
  class Moment {
  private:
    T density;
    MathVector<T, L::dimD> velocity;
    T entropy;

  public:
#pragma omp declare simd
    inline void computeMoments(const T * RESTRICT f,
                               const MathVector<unsigned int, 3>& iP) {
      computeDensity(f, iP);
      computeVelocity(f, iP, density);
    }

    #pragma omp declare simd
    inline const T& getDensity() {
      return density;
    }

#pragma omp declare simd
    inline const MathVector<T, L::dimD>& getVelocity() {
      return velocity;
    }

#pragma omp declare simd
    inline void computeDensity(const T * RESTRICT f,
                               const MathVector<unsigned int, 3>& iP) {
      density = f[hD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)];

      UnrolledFor<1, L::dimQ>::Do([&] (unsigned int iQ) {
          density += f[hD::getIndex(iP-uiL::celerity()[iQ], iQ)];
      });
    }

    #pragma omp declare simd
     inline void computeVelocity(const T * RESTRICT f,
                                const MathVector<unsigned int, 3>& iP,
                                const T density_in) {
      velocity = L::celerity()[0]
        * f[hD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)];

      UnrolledFor<1, L::dimQ>::Do([&] (unsigned int iQ) {
          velocity += L::celerity()[iQ]
            * f[hD::getIndex(iP-uiL::celerity()[iQ], iQ)];
        });

      velocity /= density_in;
    }

    #pragma omp declare simd
    inline T computeEntropy(const T * RESTRICT f,
                            const MathVector<unsigned int, 3>& iP) {
      entropy = f[hD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)]
        * log(f[hD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)]
              /L::weight()[0]);

      UnrolledFor<1, L::dimQ>::Do([&] (unsigned int iQ) {
          int indexPop_iQ = hD::getIndex(iP-uiL::celerity()[iQ], iQ);
          entropy += f[indexPop_iQ]
            * log(f[indexPop_iQ]/L::weight()[iQ]);
        });
    }

  };

}

#endif // MOMENT_H
