#ifndef MOMENT_H
#define MOMENT_H

#include <omp.h>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/sources/severity_feature.hpp>
namespace logging = boost::log;

#include "Lattice.h"
#include "MathVector.h"
#include "Helpers.h"

// Une seule boucle? (unrolled = inutile..?)
// Liste de Moments a calculer

namespace lbm {

  template <class T>
  class Moment {

    #pragma omp declare simd
    T computeDensity(const T * __restrict__ f, const int index) {
      T densityR = f[indexPop(index, 0)];

      UnrolledFor<1, L::dimQ>::Do([&] (int iQ) {
          newDensity += f[indexPop(index, iQ)];
        });

      BOOST_LOG_TRIVIAL(debug) << " - Computing density: "
                               << densityR;
      return densityR;
    }

    #pragma omp declare simd
    inline MathVector<T, L::dimD> computeVelocity(const T * __restrict__ f,
                                                  const int index,
                                                  const T density) {
      MathVector<T, L::dimD> velocityR = L::celerity()[0] * f[indexPop(index, 0)];

      UnrolledFor<1, L::dimQ>::Do([&] (int iQ) {
          velocityR += L::celerity()[iQ] * f[indexPop(index, iQ)];
        });
      velocityR /= density;

      BOOST_LOG_TRIVIAL(debug) << " - Computing velocity: "
                               << velocityR;
      return velocityR;
    }

    #pragma omp declare simd
    inline T computeEntropy(const T * __restrict__ f,
                            const int index) {
      T entropyR = 0.0;

      UnrolledFor<1, L::dimQ>::Do([&] (int iQ) {
          int indexPop = indexPop_gF(iX, iY, iQ);
          entropyR += distribution[indexPop] * log(distribution[indexPop]/weight[iQ]);
        });

      return entropyR;
    }

  };

  template<class T>
  class Moments {
  private:
    StaticArray<Moment<T, MomentType::Generic>, numberMoments> momentsArray;
    const MathVector<unsigned int, 3> offset;

  public:
  Moments(const StaticArray<Moment<T, MomentType::Generic>, numberMoments> momentsArray_in,
         const MathVector<unsigned int, 3> offset_in = MathVector<unsigned int, 3>{{0}})
      : momentsArray(momentsArray_in)
      , offset(offset_in)
    {}


    #pragma omp declare simd
    inline MathVector<T, L::dimD> moment(const MathVector<int, 3>& iP) {
      BOOST_LOG_TRIVIAL(debug) << " - Computing moment.";
      MathVector<T, L::dimD> momentR{{(T) 0}};

      for(Moment<T, MomentType::Generic> moment : momentsArray) {
        momentR += moment.moment(offset + iP);
      }

      return momentR;
    }

  };

}

#endif // MOMENT_H
