#ifndef FORCINGSCHEME_H
#define FORCINGSCHEME_H

#include <cmath>
#include <omp.h>

#include "Options.h"
#include "Commons.h"
#include "MathVector.h"
#include "Equilibrium.h"
#include "Helpers.h"

namespace lbm {
  template <class T, ForcingSchemeType forcingSchemeType>
    class ForcingScheme {};

  template <class T>
    class ForcingScheme<T, ForcingSchemeType::Generic> {
  protected:
    T tau;
    T density;
    MathVector<T, L::dimD> velocity;
    T velocity2;

  ForcingScheme(const T tau_in)
    : tau(tau_in)
    , density( (T) 0)
    , velocity(MathVector<T, L::dimD>{{ (T) 0}})
    , velocity2(0)
    {}

  public:
    inline void setDensity(const T density_in) {
      density = density_in;
    }

    inline void setVelocity(const MathVector<T, L::dimD>& velocity_in) {
      velocity = velocity_in;
      velocity2 = velocity_in.norm2();
    }

    inline void setVariables(const MathVector<T, L::dimD>& force,
                             const T density_in, const MathVector<T, L::dimD>& velocity_in) {
      setDensity(density_in);
      setVelocity(velocity_in);
    }

#pragma omp declare simd
    inline MathVector<T, L::dimD> getHydrodynamicVelocity(const MathVector<T, L::dimD>& force) const {
      return velocity + 0.5/density * force;
    }

  };

  template <class T>
    class ForcingScheme<T, ForcingSchemeType::Guo>
    : public ForcingScheme<T, ForcingSchemeType::Generic> {
  private:
    using ForcingScheme<T, ForcingSchemeType::Generic>::tau;
    using ForcingScheme<T, ForcingSchemeType::Generic>::density;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity2;

  public:
    using ForcingScheme<T, ForcingSchemeType::Generic>::ForcingScheme;

    using ForcingScheme<T, ForcingSchemeType::Generic>::setDensity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setVelocity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setVariables;

    using ForcingScheme<T, ForcingSchemeType::Generic>::getVelocityHydroForcing;

#pragma omp declare simd
    inline MathVector<T, L::dimD> getEquilibriumVelocity(const MathVector<T, L::dimD>& force) const {
      return velocity + 0.5/density * force;
    }

#pragma omp declare simd
    inline T getCollisionSource(const MathVector<T, L::dimD>& force,
                                const unsigned int iQ) const {
      T celerity_iQDotVelocity = L::celerity()[iQ].dot(velocity);

      T collisionForcingR = (L::celerity()[iQ] - velocity
                             + L::inv_cs2 * celerity_iQDotVelocity
                             * L::celerity()[iQ]).dot(force);

      return (1.0 - 1.0/(2.0 * tau) ) * L::weight()[iQ] * L::inv_cs2 * collisionForcingR;
    }

  };

  template <class T>
    class ForcingScheme<T, ForcingSchemeType::ShanChen>
    : public ForcingScheme<T, ForcingSchemeType::Generic> {
  private:
    using ForcingScheme<T, ForcingSchemeType::Generic>::tau;
    using ForcingScheme<T, ForcingSchemeType::Generic>::density;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity2;

  public:
    using ForcingScheme<T, ForcingSchemeType::Generic>::ForcingScheme;

    using ForcingScheme<T, ForcingSchemeType::Generic>::setDensity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setVelocity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setVariables;

    using ForcingScheme<T, ForcingSchemeType::Generic>::getVelocityHydroForcing;

#pragma omp declare simd
    inline MathVector<T, L::dimD> getEquilibriumVelocity(const MathVector<T, L::dimD>& force) const {
      return velocity + tau/density * force;
    }


#pragma omp declare simd
    inline T getCollisionSource(const MathVector<T, L::dimD>& force,
                                const unsigned int iQ) const {
      return 0.0;
    }
  };

  template <class T>
    class ForcingScheme<T, ForcingSchemeType::ExactDifferenceMethod>
    : public ForcingScheme<T, ForcingSchemeType::Generic> {
  private:
    using ForcingScheme<T, ForcingSchemeType::Generic>::density;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::velocity2;

    Equilibrium_ equilibrium;

    Equilibrium_ deltaEquilibrium;
    MathVector<T, L::dimD> deltaVelocity;

  public:
    ForcingScheme(const T& tau_in)
      : ForcingScheme<T, ForcingSchemeType::Generic>(tau_in)
      , equilibrium()
      , deltaEquilibrium()
      , deltaVelocity()
    {}

    using ForcingScheme<T, ForcingSchemeType::Generic>::setDensity;
    using ForcingScheme<T, ForcingSchemeType::Generic>::setVelocity;

    inline void setVariables(const MathVector<T, L::dimD>& force,
                             const T density_in, const MathVector<T, L::dimD>& velocity_in) {
      setDensity(density_in);
      setVelocity(velocity_in);
      equilibrium.setVariables(density, velocity);

      deltaVelocity = velocity + 1.0/density * force;
      deltaEquilibrium.setVariables(density, deltaVelocity);
    }

#pragma omp declare simd
    inline MathVector<T, L::dimD> getEquilibriumVelocity(const MathVector<T, L::dimD>& force) const {
      return velocity;
    }

#pragma omp declare simd
    inline T getCollisionSource(const MathVector<T, L::dimD>& force,
                                const unsigned int iQ) const {

      return deltaEquilibrium.compute(iQ)
          - equilibrium.compute(iQ);
    }

  };

  typedef ForcingScheme<dataT, forcingSchemeT> ForcingScheme_;

}

#endif // FORCINGSCHEME_H
