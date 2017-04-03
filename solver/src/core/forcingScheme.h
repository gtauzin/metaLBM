#ifndef FORCINGSCHEME_H
#define FORCINGSCHEME_H

#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include <math.h>
#include <omp.h>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include "input.h"
#include "structure.h"
#include "commons.h"
#include "calculate.h"

namespace lbm {

  template <class T>
    class ForcingScheme {
  public:
    MathVector<T, L::dimD> force;

    ForcingScheme<T>()
      : force()
      {}

    virtual MathVector<T, L::dimD> getEqVelocityForcing(const T density) const = 0;

    virtual T getCollisionForcing(const int i, const T density,
                                  const MathVector<T, L::dimD>& velocity,
                                  const T velocity2) const = 0;

    virtual MathVector<T, L::dimD> getHydroVelocityForcing(const T density) const = 0;

  };

  template <class T>
    class Guo : public ForcingScheme<T> {
  public:
    using ForcingScheme<T>::force;

    Guo<T>()
      : ForcingScheme<T>()
      {}

#pragma omp declare simd
    inline MathVector<T, L::dimD> getEqVelocityForcing(const T density) const {
      return 0.5/density * force;
    }


#pragma omp declare simd
    inline T getCollisionForcing(const int iQ, const T density,
                                 const MathVector<T, L::dimD>& velocity,
                                 const T velocity2) const {
      T celerity_iQDotVelocity = L::celerity()[iQ].dot(velocity);

      T collisionForcingR = (L::celerity()[iQ] - velocity
                             + L::inv_cs2 * celerity_iQDotVelocity
                             * L::celerity()[iQ]).dot(force);

      return (1.0 - beta) * L::weight()[iQ] * L::inv_cs2 * collisionForcingR;
    }

#pragma omp declare simd
    inline MathVector<T, L::dimD> getHydroVelocityForcing(const T density) const {
      return 0.5/density * force;
    }

  };

  template <class T>
    class ShanChen : public ForcingScheme<T> {

  public:
    using ForcingScheme<T>::force;

    ShanChen<T>()
      : ForcingScheme<T>()
      {}

#pragma omp declare simd
    inline MathVector<T, L::dimD> getEqVelocityForcing(const T density) const {
      return tau/density * force;
    }


#pragma omp declare simd
    inline T getCollisionForcing(const int iQ, const T density,
                                 const MathVector<T, L::dimD>& velocity,
                                 const T velocity2) const {

      return 0.0;
    }

#pragma omp declare simd
    inline MathVector<T, L::dimD> getHydroVelocityForcing(const T density) const {
      return 0.5/density * force;
    }

  };

  template <class T>
    class ExactDifferenceMethod : public ForcingScheme<T> {
  public:
    using ForcingScheme<T>::force;

    ExactDifferenceMethod<T>()
      : ForcingScheme<T>()
      {}

#pragma omp declare simd
    inline MathVector<T, L::dimD> getEqVelocityForcing(const T density) const {
      return MathVector<T, L::dimD>();
    }

#pragma omp declare simd
    inline T getCollisionForcing(const int i, const T density,
                                 const MathVector<T, L::dimD>& velocity,
                                 const T velocity2) const {

      return computeEquilibrium<T>(i, density,
                                      velocity + 1.0/density * force,
                                      velocity2)
        - computeEquilibrium<T>(i, density,
                                   velocity,
                                   velocity2);
    }

#pragma omp declare simd
    inline MathVector<T, L::dimD> getHydroVelocityForcing(const T density) const {
      return 0.5/density * force;
    }

  };


  template<class T>
    std::shared_ptr<ForcingScheme<T>> Create(const ForcingSchemeMethod& forcingSchemeMethod) {
    switch(forcingSchemeMethod){
    case ForcingSchemeMethod::guo : {
      return std::shared_ptr<ForcingScheme<T>>(new Guo<T>());
    }
    case ForcingSchemeMethod::SC : {
      return std::shared_ptr<ForcingScheme<T>>(new ShanChen<T>());
    }
    case ForcingSchemeMethod::EDM : {
      return std::shared_ptr<ForcingScheme<T>>(new ExactDifferenceMethod<T>());
    }
    default:{
      BOOST_LOG_TRIVIAL(error) << "Unknown type of forcing.";
      return nullptr;
    }
    }
  }

}

#endif // FORCINGSCHEME_H
