#ifndef FORCE_H
#define FORCE_H

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

  template <class T, LatticeType L>
    class Forcing {
  public:
    MathVector<T, P::dimD> force;

    Forcing<T, L>()
      : force()
      {}

    virtual MathVector<T, P::dimD> getEqVelocityForcing(const T density) const = 0;

    virtual T getCollisionForcing(const int i, const T density,
                                  const MathVector<T, P::dimD>& velocity,
                                  const T velocity2) const = 0;

    virtual MathVector<T, P::dimD> getHydroVelocityForcing(const T density) const = 0;

  };

  template <class T, LatticeType L>
    class ForcingGuo : public Forcing<T, L> {
  public:
    using Forcing<T, L>::force;

    ForcingGuo<T, L>()
      : Forcing<T, L>()
      {}

#pragma omp declare simd
    inline MathVector<T, P::dimD> getEqVelocityForcing(const T density) const {
      return 0.5/density * force;
    }


#pragma omp declare simd
    inline T getCollisionForcing(const int iQ, const T density,
                                 const MathVector<T, P::dimD>& velocity,
                                 const T velocity2) const {
      T celerity_iQDotVelocity = P::celerity()[iQ].dot(velocity);

      T collisionForcingR = (P::celerity()[iQ] - velocity
                             + P::inv_cs2 * celerity_iQDotVelocity
                             * P::celerity()[iQ]).dot(force);

      return (1.0 - beta) * P::weight()[iQ] * P::inv_cs2 * collisionForcingR;
    }

#pragma omp declare simd
    inline MathVector<T, P::dimD> getHydroVelocityForcing(const T density) const {
      return 0.5/density * force;
    }

  };

  template <class T, LatticeType L>
    class ForcingSC : public Forcing<T, L> {

  public:
    using Forcing<T, L>::force;

    ForcingSC<T, L>()
      : Forcing<T, L>()
      {}

#pragma omp declare simd
    inline MathVector<T, P::dimD> getEqVelocityForcing(const T density) const {
      return tau/density * force;
    }


#pragma omp declare simd
    inline T getCollisionForcing(const int iQ, const T density,
                                 const MathVector<T, P::dimD>& velocity,
                                 const T velocity2) const {

      return 0.0;
    }

#pragma omp declare simd
    inline MathVector<T, P::dimD> getHydroVelocityForcing(const T density) const {
      return 0.5/density * force;
    }

  };

  template <class T, LatticeType L>
    class ForcingEDM : public Forcing<T, L> {
  public:
    using Forcing<T, L>::force;

    ForcingEDM<T, L>()
      : Forcing<T, L>()
      {}

#pragma omp declare simd
    inline MathVector<T, P::dimD> getEqVelocityForcing(const T density) const {
      return MathVector<T, P::dimD>();
    }

#pragma omp declare simd
    inline T getCollisionForcing(const int i, const T density,
                                 const MathVector<T, P::dimD>& velocity,
                                 const T velocity2) const {

      return computeEquilibrium<T, L>(i, density,
                                      velocity + 1.0/density * force,
                                      velocity2)
        - computeEquilibrium<T, L>(i, density,
                                   velocity,
                                   velocity2);
    }

#pragma omp declare simd
    inline MathVector<T, P::dimD> getHydroVelocityForcing(const T density) const {
      return 0.5/density * force;
    }

  };


  template<class T, LatticeType L>
    std::shared_ptr<Forcing<T, L>> Create(const ForcingMethod& forcingMethod) {
    switch(forcingMethod){
    case ForcingMethod::guo : {
      return std::shared_ptr<Forcing<T, L>>(new ForcingGuo<T, L>());
    }
    case ForcingMethod::SC : {
      return std::shared_ptr<Forcing<T, L>>(new ForcingSC<T, L>());
    }
    case ForcingMethod::EDM : {
      return std::shared_ptr<Forcing<T, L>>(new ForcingEDM<T, L>());
    }
    default:{
      BOOST_LOG_TRIVIAL(error) << "Unknown type of forcing.";
      return nullptr;
    }
    }
  }

  template<class T, LatticeType L>
    class Force {

  public:
    Force()
      {}

#pragma omp declare simd
    virtual MathVector<T, P::dimD> force(const int iX, const int iY, const int iZ) = 0;

  };


  template<class T, LatticeType L>
    class ConstantForce : public Force<T, L> {
  private:
    const MathVector<T, P::dimD> amplitude;

  public:

  ConstantForce(const MathVector<T, P::dimD> amplitude_in)
    : Force<T, L>()
      , amplitude(amplitude_in)
    {}

#pragma omp declare simd
    inline MathVector<T, P::dimD> force(const int iX, const int iY, const int iZ) {
      return amplitude;
    }
  };


  template<class T, LatticeType L>
    class SinusoidalForce : public Force<T, L> {
  private:
    const MathVector<T, P::dimD> amplitude;
    const MathVector<T, P::dimD> waveLength;

  public:

  SinusoidalForce(const MathVector<T, P::dimD>& amplitude_in,
                  const MathVector<T, P::dimD>& waveLength_in)
    : Force<T, L>()
      , amplitude(amplitude_in)
      , waveLength(waveLength_in)
    {}

#pragma omp declare simd
    inline MathVector<T, P::dimD> force(const int iX, const int iY, const int iZ){
      MathVector<T, P::dimD> forceR;

      forceR[0] = amplitude[0] * sin(iX*2*M_PI/waveLength[0]);
      forceR[1] = amplitude[1] * sin(iY*2*M_PI/waveLength[1]);
      forceR[2] = amplitude[2] * sin(iZ*2*M_PI/waveLength[2]);

      return forceR;
    }

  };

  template<class T, LatticeType L>
    std::shared_ptr<Force<T, L>> Create(const ForceType& forceType,
                                        const MathVector<T, P::dimD>& amplitude,
                                        const MathVector<T, P::dimD>& waveLength) {

    switch(forceType){
    case ForceType::constant : {
      return std::shared_ptr<Force<T, L>>(new ConstantForce<T, L>(amplitude));
    }
    case ForceType::sinusoidal : {
      return std::shared_ptr<Force<T, L>>(new SinusoidalForce<T, L>(amplitude,
                                                                    waveLength));
    }

    default:{
      BOOST_LOG_TRIVIAL(error) << "Unknown type of force.";
      return nullptr;
    }
    }
  }


  template<class T, LatticeType L>
    class Forces {
  private:
    std::array<std::shared_ptr<Force<T, L>>, numberForces> forcesArray;

  public:
  Forces(std::array<std::shared_ptr<Force<T, L>>, numberForces> forcesArray_in)
    : forcesArray(forcesArray_in)
    {}


#pragma omp declare simd
    inline MathVector<T, P::dimD> force(const int iX, const int iY, const int iZ) {
      MathVector<T, P::dimD> forceR = MathVector<T, P::dimD>{{0.0}};

      for(std::shared_ptr<Force<T, L>> bodyForce : forcesArray) {
        forceR += bodyForce->force(iX, iY, iZ);
      }
      return forceR;
    }

  };

}

#endif // FORCE_H
