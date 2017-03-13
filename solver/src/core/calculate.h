#ifndef CALCULATE_H
#define CALCULATE_H

#include <iostream>
#include <random>
#include <array>
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

namespace lbm {

#pragma omp declare simd
  template <class T, LatticeType L>
    T computeDensity(const T * __restrict__ f, const int idx) {
    T newDensity = f[idxPop<T, L>(idx, 0)];

    UnrolledFor<1, dimQ<T, L>()>::Do([&] (int iQ) {
        newDensity += f[idxPop<T, L>(idx, iQ)];
      });
    return newDensity;
  }


#pragma omp declare simd
  template <class T, LatticeType L>
    inline MathVector<T, dimD<T, L>()> computeVelocity(const T * __restrict__ f,
                                                       const int idx,
                                                       const T density) {
    MathVector<T, dimD<T, L>()> velocityR(celerity<T, L>(0)* f[idxPop<T, L>(idx, 0)]);

    UnrolledFor<1, dimQ<T, L>()>::Do([&] (int iQ) {
        velocityR += celerity<T, L>(iQ)* f[idxPop<T, L>(idx, iQ)];
      });

    return velocityR/density;
  }

  #pragma omp declare simd
  template<class T, LatticeType L>
    inline void calculateMoments(const T * __restrict__ f, const int idx_lattice,
                                 T& density, MathVector<T, dimD<T, L>()>& velocity) {

    BOOST_LOG_TRIVIAL(debug) << " - Computing density.";
    density = computeDensity<T, L>(f, idx_lattice);

    BOOST_LOG_TRIVIAL(debug) << " - Computing velocity.";
    velocity = computeVelocity<T, L>(f, idx_lattice, density);
  }


  template<class T, LatticeType L>
    void calculateMomentsField(Lattice<T, L>& l_previous,
                               LocalField<T, L>& field) {
#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(int iZ = hZ<T, L>(); iZ < hZ<T, L>()+lZ_l<T, L>(); ++iZ) {
      for(int iY = hY<T, L>(); iY < hY<T, L>()+lY_l<T, L>(); ++iY) {
#pragma omp simd
        for(int iX = hX<T, L>(); iX < hX<T, L>()+lX_l<T, L>(); ++iX) {
          BOOST_LOG_TRIVIAL(debug) << " - (x, y, z) = "
                                   << "(" << iX
                                   << ", " << iY
                                   << ", " << iZ << ")";
          int idx_lattice = idxL<T, L>(iX, iY, iZ);
          int idx_field = idx_inF<T, L>(iX, iY, iZ);

          T previousDensity;
          MathVector<T, dimD<T, L>()> previousVelocity;
          calculateMoments<T, L>(l_previous.f_distribution.data(), idx_lattice,
                                 previousDensity, previousVelocity);

          field.nextDensity[idx_field] = previousDensity;
          field.nextVelocity[idx_field] = previousVelocity;
        }
      }
    }
  }


#pragma omp declare simd
  template <class T, LatticeType L>
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
  template <class T, LatticeType L>
    inline T computeEquilibrium(const int iQ, T rho,
                                const MathVector<T, dimD<T, L>()>& u,
                                const T u2) {

    T fEq_iQ = 1.0;

    UnrolledFor<0, dimD<T, L>()>::Do([&] (int d) {
        fEq_iQ *= (2.0 - sqrt(1.0 + 3.0*u[d]*u[d]))
          * powerOfCelerity<T, L>((2* u[d] + sqrt(1.0 + 3.0*u[d]*u[d]))/(1.0 - u[d]),
                                  (int) celerity<T, L>(d, iQ));
      });

    return rho * weight<T, L>(iQ)*fEq_iQ;
  }

}

#endif // CALCULATE_H
