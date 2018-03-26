#ifndef COLLISION_H
#define COLLISION_H

#include "Commons.h"
#include "Options.h"
#include "Domain.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Helpers.h"
#include "Moment.h"
#include "Equilibrium.h"
#include "Force.h"
#include "ForcingScheme.h"
#include "EntropicStep.h"

namespace lbm {

  template <class T, CollisionType collisionType>
  class Collision {};

  template <class T>
  class Collision<T, CollisionType::GenericSRT> {
  protected:
    const T tau;

    Force_ forcing;
    ForcingScheme_ forcingScheme;

    T density;
    MathVector<T, L::dimD> velocity;
    T velocity2;
    MathVector<T, L::dimD> force;
    T entropy;


    Collision(const T tau_in,
              const MathVector<T, 3>& amplitude_in,
              const MathVector<T, 3>& waveLength_in,
              const unsigned int kMin_in, const unsigned int kMax_in)
      : tau(tau_in)
      , forcing(amplitude_in, waveLength_in, kMin_in, kMax_in)
      , forcingScheme(tau_in)
      , density()
      , velocity{{0}}
      , velocity2()
      , force{{0}}
      , entropy()
    {}

  public:
    DEVICE HOST
    inline const T& getDensity() {
      return density;
    }

    DEVICE HOST
    inline const MathVector<T, L::dimD>& getVelocity() {
      return velocity;
    }

    DEVICE HOST
    inline const MathVector<T, L::dimD>& getForce() {
      return force;
    }

    DEVICE HOST
    inline void calculateMoments(const T * haloDistributionPtr,
                                 const Position& iP) {
      { INSTRUMENT_OFF("Moment<T>::calculateMoments",4) }

      Moment_::calculateDensity(haloDistributionPtr, iP, density);
      Moment_::calculateVelocity(haloDistributionPtr, iP, density, velocity);

      velocity2 = velocity.norm2();
    }

    DEVICE HOST
    inline void setForce(T * localForceArray,
                         const Position& iP,
                         const Position& offset) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::setForce",4) }

      forcing.setForce(localForceArray, iP-L::halo(), force);
      forcingScheme.setVariables(force, density, velocity);
    }

    DEVICE HOST
    inline const MathVector<T, L::dimD> getHydrodynamicVelocity() {
      return forcingScheme.calculateHydrodynamicVelocity(force, density, velocity);
    }

    DEVICE HOST
    inline void calculateRelaxationTime(const T * haloDistributionNext_Ptr,
                                        const T * haloDistributionPrevious_Ptr,
                                        const Position& iP) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::calculate",4) }
    }

    DEVICE HOST
    inline T calculate(const T * haloDistributionPtr,
                       const Position& iP,
                       const unsigned int iQ) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::calculate",4) }
      T equilibrium_iQ = Equilibrium_::calculate(density, velocity, velocity2, iQ);

        return ( (T) 1.0 - (T) 1.0/tau) * haloDistributionPtr[hSD::getIndex(iP, iQ)]
        + forcingScheme.calculateCollisionSource(force, density, velocity, velocity2,
                                                 equilibrium_iQ, iQ)
        + (T) 1.0/tau * equilibrium_iQ;
    }

    DEVICE HOST
    inline void update(const unsigned int iteration) {
      forcing.update(iteration);
    }

  };


  template <class T>
  class Collision<T, CollisionType::BGK>
    : public Collision<T, CollisionType::GenericSRT> {
  private:
    using Base = Collision<T, CollisionType::GenericSRT>;
  protected:
    using Base::tau;
    using Base::forcingScheme;

    T alpha;

  public:
    Collision(const T tau_in,
              const MathVector<T, 3>& amplitude_in,
              const MathVector<T, 3>& waveLength_in,
              const unsigned int kMin_in, const unsigned int kMax_in)
      : Base(tau_in, amplitude_in, waveLength_in, kMin_in, kMax_in)
      , alpha( (T) 2)
    {}

    using Base::update;
    using Base::setForce;

    using Base::calculate;

    using Base::getHydrodynamicVelocity;
    using Base::getDensity;
    using Base::getVelocity;
    using Base::getForce;

    DEVICE HOST
    inline T getAlpha() {
      return alpha;
    }
  };


  template <class T>
  class Collision<T, CollisionType::ELBM>
    : public Collision<T, CollisionType::BGK> {
  private:
    using Base = Collision<T, CollisionType::BGK>;

  public:
    Collision(const T tau_in,
              const MathVector<T, 3>& amplitude_in,
              const MathVector<T, 3>& waveLength_in,
              const unsigned int kMin_in, const unsigned int kMax_in)
      : Base(tau_in, amplitude_in, waveLength_in,
                                         kMin_in, kMax_in)
      , beta( (T) 1.0/(2.0 * tau_in))
    {}

    using Base::setForce;

    DEVICE HOST
    inline void setVariables(const T * haloDistributionPtr,
                             const Position& iP,
                             const T density,
                             const MathVector<T, L::dimD>& velocity) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::setVariables",4) }

        Base::setVariables(haloDistributionPtr, iP, density, velocity);

      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        T equilibrium_iQ = Equilibrium_::calculate(Base::density, Base::velocity,
                                                   Base::velocity2, iQ);

        f_Forced[iQ] = haloDistributionPtr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
          + forcingScheme.calculateCollisionSource(Base::force, Base::density,
                                                   Base::velocity, Base::velocity2,
                                                   equilibrium_iQ, iQ);
        f_NonEq[iQ] = haloDistributionPtr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
          - equilibrium_iQ;
      }

      calculateAlpha();
      tau = (T) 1.0/(Base::alpha*beta);
    }

    using Base::update;
    using Base::calculate;

    using Base::getHydrodynamicVelocity;
    using Base::getDensity;
    using Base::getVelocity;
    using Base::getForce;
    using Base::getAlpha;

  protected:
    using Base::forcingScheme;
    using Base::equilibrium;

    using Base::tau;

    const T beta;

    MathVector<T, L::dimQ> f_Forced;
    MathVector<T, L::dimQ> f_NonEq;

    DEVICE HOST
    inline bool isDeviationSmall(const T error) {
      INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::isDeviationSmall",6)

        bool isDeviationSmallR = true;
      T deviation;

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        deviation = fabs(f_NonEq[iQ]/f_Forced[iQ]);

        if(deviation > error) {
          isDeviationSmallR = false;
        }
      }

      return isDeviationSmallR;
    }

    DEVICE HOST
    T calculateAlphaMax() {
      INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::calculateAlphaMax",6)

        T alphaMaxR = 2.5;
      T alphaMaxTemp;

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        if(f_NonEq[iQ] > 0) {
          alphaMaxTemp = fabs(f_Forced[iQ]/f_NonEq[iQ]);

          if(alphaMaxTemp < alphaMaxR) {
            alphaMaxR = alphaMaxTemp;
          }
        }
      }

      return alphaMaxR;
    }

    DEVICE HOST
    inline T solveAlpha(const T alphaMin, const T alphaMax) {
      INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::solveAlpha",6)

        std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T>(f_Forced, f_NonEq));
      const T tolerance = 1e-5;
      const int iterationMax = 20;
      T alphaR = Base::alpha;

      bool hasConverged = NewtonRaphsonSolver(entropicStepFunctor,
                                              tolerance, iterationMax,
                                              alphaR, alphaMin, alphaMax);

      if(!hasConverged) {
        return 2.0;
      }

      return alphaR;
    }

    DEVICE HOST
    inline void calculateAlpha() {
      INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::calculateAlpha",5)

        if(isDeviationSmall( (T) 1.0e-3)) {
          Base::alpha = 2.0;
        }

        else {
          T alphaMax = calculateAlphaMax();

          if(alphaMax < 2.) {
            Base::alpha = 0.95 * alphaMax;
          }

          else {
            T alphaMin = 1.;
            Base::alpha = solveAlpha(alphaMin, alphaMax);

          }
        }
    }
  };


  template <class T>
  class Collision<T, CollisionType::Approached_ELBM>
    : public Collision<T, CollisionType::ELBM> {
  private:
    using Base = Collision<T, CollisionType::ELBM>;

  public:
    Collision(const T tau_in,
              const MathVector<T, 3>& amplitude_in,
              const MathVector<T, 3>& waveLength_in,
              const unsigned int kMin_in, const unsigned int kMax_in)
      : Base(tau_in, amplitude_in, waveLength_in,
                                          kMin_in, kMax_in)
    {}

    using Base::update;
    using Base::setForce;
    using Base::setVariables;

    using Base::calculate;

    using Base::getHydrodynamicVelocity;
    using Base::getDensity;
    using Base::getVelocity;
    using Base::getForce;
    using Base::getAlpha;

  private:
    using Base::tau;
    using Base::alpha;
    using Base::forcingScheme;
    using Base::equilibrium;
    using Base::f_Forced;
    using Base::f_NonEq;

    using Base::isDeviationSmall;
    using Base::calculateAlphaMax;
    using Base::solveAlpha;

    DEVICE HOST
    inline T approximateAlpha() {
      { INSTRUMENT_OFF("Collision<T, CollisionType::Approached_ELBM>::approximateAlpha",6) }

      T a1 = 0.0;
      T a2 = 0.0;
      T a3 = 0.0;
      T a4 = 0.0;

      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        T temp = f_NonEq[iQ]/f_Forced[iQ];
        a1 += f_NonEq[iQ]*temp;
        a2 += f_NonEq[iQ]*temp*temp;
        a3 += f_NonEq[iQ]*temp*temp*temp;
        a4 += f_NonEq[iQ]*temp*temp*temp*temp;
      }

      a1 *= 1.0/2.0;
      a2 *= 1.0/6.0;
      a3 *= 1.0/12.0;
      a4 *= 1.0/20.0;

      T alphaR = 2 - 1/a1 * (4.0*a2 + 16.0*a2*a2/a1 - 8.0*a3
                             + 80.0*a2*a3/a1 - 80.0*a2*a2*a2/(a1*a1) - 16.0*a4);

      return alphaR;
    }

    DEVICE HOST
    inline void calculateAlpha() {
      { INSTRUMENT_OFF("Collision<T, CollisionType::Approached_ELBM>::calculateAlpha",5) }

        if(isRelativeDeviationSmall( (T) 1.0e-3)) {
          T alphaApproximated = approximateAlpha();
          alpha = alphaApproximated;
        }
        else {
          T alphaMax = calculateAlphaMax();

          if(alphaMax < 2.) {
            alpha = 0.95 * alphaMax;
          }

          else {
            T alphaMin = 1.;
            alpha = solveAlpha(alphaMin, alphaMax);

          }
        }
    }
  };


  template <class T>
  class Collision<T, CollisionType::ForcedNR_ELBM>
    : public Collision<T, CollisionType::ELBM> {
  private:
    using Base = Collision<T, CollisionType::ELBM>;

  public:
    using Base::update;
    using Base::setForce;
    using Base::setVariables;

    using Base::calculate;

    using Base::getHydrodynamicVelocity;
    using Base::getDensity;
    using Base::getVelocity;
    using Base::getForce;
    using Base::getAlpha;

  private:
    using Base::alpha;
    using Base::forcingScheme;
    using Base::equilibrium;

    using Base::f_Forced;
    using Base::f_NonEq;

    using Base::isDeviationSmall;
    using Base::calculateAlphaMax;

    DEVICE HOST
    inline T solveAlpha(const T alphaMin, const T alphaMax) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::ForcedNR_ELBM>::solveAlpha",6) }

      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T>(f_Forced, f_NonEq));
      const T tolerance = 1e-5;
      const int iterationMax = 20;
      T alphaR = alpha;

      bool hasConverged = NewtonRaphsonSolver(entropicStepFunctor, tolerance, iterationMax,
                                              alphaR, alphaMin, alphaMax);

      if(!hasConverged) {
        return 2.0;
      }

      return alphaR;
    }

    DEVICE HOST
    inline T calculateAlpha() {
      INSTRUMENT_OFF("Collision<T, CollisionType::ForcedNR_ELBM>::calculateAlpha",5)

        T alphaMax = calculateAlphaMax();

      if(alphaMax < 2.) {
        alpha = 0.95 * alphaMax;
      }

      else {
        T alphaMin = 1.;
        alpha = solveAlpha(alphaMin, alphaMax);
      }
    }
  };

  template <class T>
  class Collision<T, CollisionType::ForcedBNR_ELBM>
    : public Collision<T, CollisionType::ForcedNR_ELBM> {
  private:
    using Base = Collision<T, CollisionType::ForcedNR_ELBM>;

  public:
    using Base::update;
    using Base::setForce;
    using Base::setVariables;

    using Base::calculate;

    using Base::getHydrodynamicVelocity;
    using Base::getDensity;
    using Base::getVelocity;
    using Base::getForce;
    using Base::getAlpha;

  private:
    using Base::alpha;
    using Base::f_Forced;
    using Base::f_NonEq;

    using Base::forcingScheme;
    using Base::equilibrium;

    using Base::isDeviationSmall;
    using Base::calculateAlphaMax;
    using Base::calculateAlpha;

    DEVICE HOST
    inline T solveAlpha(const T alphaMin, const T alphaMax) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::ForcedBNR_ELBM>::solveAlpha",6) }

      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T>(f_Forced, f_NonEq));
      const T tolerance = 1e-5;
      const int iterationMax = 20;
      T alphaR = alpha;

      bool hasConverged = Bisection_NewtonRaphsonSolver(entropicStepFunctor,
                                                        tolerance, iterationMax,
                                                        alphaR, alphaMin, alphaMax);

      if(!hasConverged) {
        return 2.0;
      }

      return alphaR;
    }
  };


  typedef Collision<dataT, collisionT> Collision_;

}

#endif // COLLISION_H
