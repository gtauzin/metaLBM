#ifndef FORCINGSCHEME_H
#define FORCINGSCHEME_H

#include <cmath>

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "Lattice.h"
#include "Equilibrium.h"
#include "Helpers.h"

namespace lbm {

  template <class T, ForcingSchemeType forcingSchemeType>
  class ForcingScheme {};

  template <class T>
  class ForcingScheme<T, ForcingSchemeType::Generic> {
  protected:
    T tau;

  public:
    ForcingScheme(const T tau_in)
      : tau(tau_in)
    {}

  public:
    #pragma omp declare simd
    DEVICE HOST
      inline MathVector<T, L::dimD> calculateHydrodynamicVelocity(const MathVector<T, L::dimD>& force, const T& density, const MathVector<T, L::dimD>& velocity) const {
      INSTRUMENT_OFF("ForcingScheme<T, ForcingSchemeType::Generic>::calculateHydrodynamicVelocity",5)

      return velocity + 0.5/density * force;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline T setVariables(const MathVector<T, L::dimD>& force,
                          const T& density,
                          const MathVector<T, L::dimD>& velocity) {
    }
  };


  template <class T>
  class ForcingScheme<T, ForcingSchemeType::Guo>
    : public ForcingScheme<T, ForcingSchemeType::Generic> {
  private:
    using Base = ForcingScheme<T, ForcingSchemeType::Generic>;

  public:
    using Base::ForcingScheme;

    #pragma omp declare simd
    DEVICE HOST
    inline MathVector<T, L::dimD>
    calculateEquilibriumVelocity(const MathVector<T, L::dimD>& force,
                                 const T& density,
                                 const MathVector<T, L::dimD>& velocity) const {
      INSTRUMENT_OFF("ForcingScheme<T, ForcingSchemeType::Guo>::calculateEquilibriumVelocity",5)

      return velocity + 0.5/density * force;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline T calculateCollisionSource(const MathVector<T, L::dimD>& force,
                                      const T& density,
                                      const MathVector<T, L::dimD>& velocity,
                                      const T& velocity2, const T equilibrium_iQ,
                                      const unsigned int iQ) const {
      INSTRUMENT_OFF("ForcingScheme<T, ForcingSchemeType::Guo>::calculateCollisionSource",5)

      T celerity_iQDotVelocity = L::celerity()[iQ].dot(velocity);

      T collisionForcingR = (L::celerity()[iQ] - velocity
                             + L::inv_cs2 * celerity_iQDotVelocity
                             * L::celerity()[iQ]).dot(force);

      return (1.0 - 1.0/(2.0 * Base::tau)) * L::weight()[iQ] * L::inv_cs2 * collisionForcingR;
    }

  };


  template <class T>
  class ForcingScheme<T, ForcingSchemeType::ShanChen>
    : public ForcingScheme<T, ForcingSchemeType::Generic> {
  private:
    using Base = ForcingScheme<T, ForcingSchemeType::Generic>;

    using Base::tau;

  public:
    using Base::ForcingScheme;
    using Base::calculateVelocityHydroForcing;
    using Base::setVariables;

    #pragma omp declare simd
    DEVICE HOST
    inline MathVector<T, L::dimD> calculateEquilibriumVelocity(const MathVector<T, L::dimD>& force, const T& density, const MathVector<T, L::dimD>& velocity) const {
      { INSTRUMENT_OFF("ForcingScheme<T, ForcingSchemeType::ShanChen>::calculateEquilibriumVelocity",5) }

      return velocity + tau/density * force;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline T calculateCollisionSource(const MathVector<T, L::dimD>& force,
                                      const T& density,
                                      const MathVector<T, L::dimD>& velocity,
                                      const T& velocity2, const T equilibrium_iQ,
                                      const unsigned int iQ) const {
      INSTRUMENT_OFF("ForcingScheme<T, ForcingSchemeType::ShanChen>::calculateCollisionSource",5)

      return 0.0;
    }

  };


  template <class T>
  class ForcingScheme<T, ForcingSchemeType::ExactDifferenceMethod>
    : public ForcingScheme<T, ForcingSchemeType::Generic> {
  private:
    using Base = ForcingScheme<T, ForcingSchemeType::Generic>;

    MathVector<T, L::dimD> deltaVelocity;
    T deltaVelocity2;

  public:
    ForcingScheme(const T& tau_in)
      : Base(tau_in)
      , deltaVelocity()
      , deltaVelocity2()
    {}

    #pragma omp declare simd
    DEVICE HOST
    inline MathVector<T, L::dimD> calculateEquilibriumVelocity(const MathVector<T, L::dimD>& force, const T& density, const MathVector<T, L::dimD>& velocity) const {
      { INSTRUMENT_OFF("ForcingScheme<T, ForcingSchemeType::ExactDifferenceMethod>::calculateEquilibriumVelocity",5) }

      return velocity;
    }

    #pragma omp declare simd
    DEVICE HOST
      inline void setVariables(const MathVector<T, L::dimD>& force,
                          const T& density,
                          const MathVector<T, L::dimD>& velocity) {
      deltaVelocity = velocity + 1.0/density * force;
      deltaVelocity2 = deltaVelocity.norm2();
    }

    #pragma omp declare simd
    DEVICE HOST
    inline T calculateCollisionSource(const MathVector<T, L::dimD>& force,
                                      const T& density,
                                      const MathVector<T, L::dimD>& velocity,
                                      const T& velocity2, const T equilibrium_iQ,
                                      const unsigned int iQ) const {
      { INSTRUMENT_OFF("ForcingScheme<T, ForcingSchemeType::ExactDifferenceMethod>::calculateCollisionSource",5) }

        return Equilibrium_::calculate(density, velocity + 1.0/density * force,
                                       (velocity + 1.0/density * force).norm2(), iQ)
          - equilibrium_iQ;
    }

  };

  typedef ForcingScheme<dataT, forcingSchemeT> ForcingScheme_;

}

#endif // FORCINGSCHEME_H
