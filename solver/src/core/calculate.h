#ifndef CALCULATE_H
#define CALCULATE_H

#include <iostream>
#include <math.h>
#include <omp.h>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/sources/severity_feature.hpp>
namespace logging = boost::log;

#include "input.h"
#include "structure.h"
#include "commons.h"
#include "field.h"
#include "distribution.h"

namespace lbm {

#pragma omp declare simd
  template <class T>
    T computeDensity(const T * __restrict__ f, const int idx) {
    T newDensity = f[idxPop(idx, 0)];

    UnrolledFor<1, L::dimQ>::Do([&] (int iQ) {
        newDensity += f[idxPop(idx, iQ)];
      });

    return newDensity;
  }

#pragma omp declare simd
  template <class T>
    inline MathVector<T, L::dimD> computeVelocity(const T * __restrict__ f,
                                                  const int idx,
                                                  const T density) {
    MathVector<T, L::dimD> velocityR = L::celerity()[0] * f[idxPop(idx, 0)];

    UnrolledFor<1, L::dimQ>::Do([&] (int iQ) {
        velocityR += L::celerity()[iQ] * f[idxPop(idx, iQ)];
      });

    return velocityR/density;
  }

#pragma omp declare simd
  template <class T>
    inline void calculateMoments(const T * __restrict__ f, const int idx_lattice,
                                 T& density, MathVector<T, L::dimD>& velocity) {

    BOOST_LOG_TRIVIAL(debug) << " - Computing density.";
    density = computeDensity<T>(f, idx_lattice);

    BOOST_LOG_TRIVIAL(debug) << " - Computing velocity.";
    velocity = computeVelocity<T>(f, idx_lattice, density);
  }


  template <class T>
    void calculateMomentsField(Distribution<T>& l_previous,
                               LocalField<T>& field) {
#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(int iZ = L::hZ; iZ < L::hZ + L::lZ_l; ++iZ) {
      for(int iY = L::hY; iY < L::hY + L::lY_l; ++iY) {
#pragma omp simd
        for(int iX = L::hX; iX < L::hX + L::lX_l; ++iX) {
          BOOST_LOG_TRIVIAL(debug) << " - (x, y, z) = "
                                   << "(" << iX
                                   << ", " << iY
                                   << ", " << iZ << ")";
          int idx_lattice = idxL(iX, iY, iZ);
          int idx_field = idx_inF(iX, iY, iZ);

          T previousDensity;
          MathVector<T, L::dimD> previousVelocity;
          calculateMoments<T>(l_previous.f_distribution.data(), idx_lattice,
                              previousDensity, previousVelocity);

          field.nextDensity[idx_field] = previousDensity;
          field.nextVelocity[idx_field] = previousVelocity;
        }
      }
    }
  }


#pragma omp declare simd
  template <class T>
    inline T powerOfCelerity(const T arg, const int celerity) {
    if(celerity == 0) {
      return 1.0;
    }
    else if(celerity == -1) {
      return 1.0/arg;
    }
    else {
      return arg;
    }
  }

#pragma omp declare simd
  template <class T>
    inline T computeEquilibrium(const int iQ, T rho,
                                const MathVector<T, L::dimD>& u,
                                const T u2) {

    T fEq_iQ = 1.0;

    UnrolledFor<0, L::dimD>::Do([&] (int d) {
        fEq_iQ *= (2.0 - sqrt(1.0 + 3.0*u[d]*u[d]))
          * powerOfCelerity<T>((2* u[d] + sqrt(1.0 + 3.0*u[d]*u[d]))/(1.0 - u[d]),
                               (int) L::celerity()[iQ][d]);
      });

    return rho * L::weight()[iQ]*fEq_iQ;
  }

}

#endif // CALCULATE_H
