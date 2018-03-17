#ifndef MOMENT_H
#define MOMENT_H

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
    DEVICE HOST
    inline void calculateMoments(const T * f,
                                 const MathVector<unsigned int, 3>& iP) {
      INSTRUMENT_OFF("Moment<T>::calculateMoments",4)

      calculateDensity(f, iP);
      calculateVelocity(f, iP, density);
    }

    #pragma omp declare simd
    DEVICE HOST
    inline const T& getDensity() {
      return density;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline const MathVector<T, L::dimD>& getVelocity() {
      return velocity;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline void calculateDensity(const T * f,
				 const MathVector<unsigned int, 3>& iP) {
      INSTRUMENT_OFF("Moment<T>::calculateDensity",5)

      density = f[hSD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)];

      for(unsigned int iQ = 1; iQ < L::dimQ; ++iQ) {
        density += f[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)];
      }
    }

    #pragma omp declare simd
    DEVICE HOST
    inline void calculateVelocity(const T * f,
                                  const MathVector<unsigned int, 3>& iP,
                                  const T density_in) {
      INSTRUMENT_OFF("Moment<T>::calculateVelocity",5)

       velocity = L::celerity()[0]
        * f[hSD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)];

       for(unsigned int iQ = 1; iQ < L::dimQ; ++iQ) {
         velocity += L::celerity()[iQ]
           * f[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)];
       }

      velocity /= density_in;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline T calculateEntropy(const T * f,
                            const MathVector<unsigned int, 3>& iP) {
      INSTRUMENT_OFF("Moment<T>::calculateEntropy",5)

      entropy = f[hSD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)]
        * log(f[hSD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)]
              /L::weight()[0]);

      for(unsigned int iQ = 1; iQ < L::dimQ; ++iQ) {
        int indexPop_iQ = hSD::getIndex(iP-uiL::celerity()[iQ], iQ);
        entropy += f[indexPop_iQ]
          * log(f[indexPop_iQ]/L::weight()[iQ]);
      }
    }

  };

}

#endif // MOMENT_H
