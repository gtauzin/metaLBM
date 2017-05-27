#ifndef FORCINGSCHEME_H
#define FORCINGSCHEME_H

#include <cmath>
#include <omp.h>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include "Options.h"
#include "MathVector.h"
#include "Helpers.h"

namespace lbm {
    template <class T, ForcingSchemeType forcingSchemeType>
    class ForcingScheme {};

  template <class T>
  class ForcingScheme<T, ForcingSchemeType::Generic> {
  protected:
    MathVector<T, L::dimD> force;
    T density;
    MathVector<T, L::dimD> velocity;
    T velocity2;

  public:
    ForcingScheme()
    : force()
    , density( (T) 0)
    , velocity(MathVector<T, L::dimD>{{ (T) 0}})
    , velocity2(0)
    {}

    inline void setForce(const MathVector<T, L::dimD>& force_in) {
      force = force_in;
    }

    inline void setDensity(const T density_in) {
      density = density_in;
    }

    inline void setVelocity(const MathVector<T, L::dimD>& velocity_in) {
      velocity = velocity_in;
      velocity2 = velocity_in.norm2();
    }

    inline void setVariables(const MathVector<T, L::dimD> force_in,
                             const T density_in, const MathVector<T, L::dimD> velocity_in) {
      setForce(force_in);
      setDensity(density_in);
      setVelocity(velocity_in);
  }

    #pragma omp declare simd
    inline MathVector<T, L::dimD> getVelocityHydroForcing
    (const MathVector<T, L::dimD>& velocity, const T density) const {
      return velocity + 0.5/density * force;
    }

  };

  template <class T>
  class ForcingScheme <T, ForcingSchemeType::Guo>
    : public ForingScheme<T, ForcingScheme::Generic> {
  private:
    using ForcingScheme<T, ForcingSchemeType::Generic>::force;
    using ForcingScheme<T, ForcingSchemeType::Generic>::density;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity2;

  public:
    using ForcingScheme<T, ForcingSchemeType::Generic>::ForcingScheme;

    using ForcingScheme<T, ForcingSchemeType::Generic>::setForce;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setDensity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setVelocity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setVariables;

    using ForcingScheme<T, ForcingSchemeType::Generic>::getVelocityHydroForcing;

    #pragma omp declare simd
    inline MathVector<T, L::dimD> getVelocityEqForcing(const MathVector<T, L::dimD>& velocity,
                                                       const T density) const {
      return velocity + 0.5/density * force;
    }

    #pragma omp declare simd
    template <unsigned int iQ>
    inline T getCollisionForcing() const {
      T celerity_iQDotVelocity = L::celerity()[iQ].dot(velocity);

      T collisionForcingR = (L::celerity()[iQ] - velocity
                             + L::inv_cs2 * celerity_iQDotVelocity
                             * L::celerity()[iQ]).dot(force);

      return (1.0 - beta) * L::weight()[iQ] * L::inv_cs2 * collisionForcingR;
    }

  };

  template <class T>
  class ForcingScheme<T, ForcingSchemeType::ShanChen>
    : public ForcingScheme<T, ForcingSchemeType::Generic> {
  private:
    using ForcingScheme<T, ForcingSchemeType::Generic>::force;
    using ForcingScheme<T, ForcingSchemeType::Generic>::density;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity2;

  public:
    using ForcingScheme<T, ForcingSchemeType::Generic>::ForcingScheme;

    using ForcingScheme<T, ForcingSchemeType::Generic>::setForce;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setDensity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setVelocity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setVariables;

    using ForcingScheme<T, ForcingSchemeType::Generic>::ForcingScheme;
    using ForcingScheme<T, ForcingSchemeType::Generic>::getVelocityHydroForcing;

    #pragma omp declare simd
    inline MathVector<T, L::dimD> getVelocityEqForcing
    (const MathVector<T, L::dimD>& velocity, const T density) const {
      return velocity + tau/density * force;
    }


    #pragma omp declare simd
    template <unsigned int iQ>
    inline T getCollisionForcing() const {

      return 0.0;
    }
  };

  template <class T>
  class ForcingScheme<T, ForcingSchemeType::ExactDifferenceMethod>
    : public ForcingScheme<T, ForcingSchemeType::Generic> {
  private:
    using ForcingScheme<T, ForcingSchemeType::Generic>::force;
    using ForcingScheme<T, ForcingSchemeType::Generic>::density;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity2;

    Equilibrium<T, L::Type, P::equilibrium::Type> deltaEquilibrium;
    MathVector<T, L::dimD> deltaVelocity;

  public:
    using ForcingScheme<T, ForcingSchemeType::Generic>::ForcingScheme;

    using ForcingScheme<T, ForcingSchemeType::Generic>::setForce;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setDensity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setVelocity;

    inline void setVariables(const MathVector<T, L::dimD>& force_in,
                             const T density_in, const MathVector<T, L::dimD>& velocity_in) {
      setForce(force_in);
      setDensity(density_in);
      setVelocity(velocity_in);

      deltaVelocity = velocity + 1.0/density * force;
      deltaEquilibrium.setVariables(density, deltaVelocity);
  }

    #pragma omp declare simd
    inline MathVector<T, L::dimD> getVelocityEqForcing(const MathVector<T, L::dimD>& velocity,
                                                       const T density) const {
      return velocity;
    }

    #pragma omp declare simd
    template <unsigned int iQ>
    inline T getCollisionForcing() const {
      return deltaEquilibium.computeEquilibrium<iQ>()
        - P::equilibrium.computeEquilibrium<iQ>();
    }

  };

}

#endif // FORCINGSCHEME_H
