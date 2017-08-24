#ifndef MOMENT_H
#define MOMENT_H

#include <omp.h>

#include "Commons.h"
#include "Options.h"
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
    inline void calculateMoments(const T * RESTRICT f,
                               const MathVector<unsigned int, 3>& iP) {
      SCOREP_INSTRUMENT_OFF("Moment<T>::calculateMoments")

      calculateDensity(f, iP);
      calculateVelocity(f, iP, density);
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
    inline void calculateDensity(const T * RESTRICT f,
                               const MathVector<unsigned int, 3>& iP) {
      SCOREP_INSTRUMENT_OFF("Moment<T>::calculateDensity")

      density = f[hD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)];

      UnrolledFor<1, L::dimQ>::Do([&] (unsigned int iQ) {
          density += f[hD::getIndex(iP-uiL::celerity()[iQ], iQ)];
      });
    }

    #pragma omp declare simd
     inline void calculateVelocity(const T * RESTRICT f,
                                const MathVector<unsigned int, 3>& iP,
                                const T density_in) {
       SCOREP_INSTRUMENT_OFF("Moment<T>::calculateVelocity")

       velocity = L::celerity()[0]
        * f[hD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)];

      UnrolledFor<1, L::dimQ>::Do([&] (unsigned int iQ) {
          velocity += L::celerity()[iQ]
            * f[hD::getIndex(iP-uiL::celerity()[iQ], iQ)];
        });

      velocity /= density_in;
    }

    #pragma omp declare simd
    inline T calculateEntropy(const T * RESTRICT f,
                            const MathVector<unsigned int, 3>& iP) {
      SCOREP_INSTRUMENT_OFF("Moment<T>::calculateEntropy")

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
