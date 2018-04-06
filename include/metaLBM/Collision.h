#pragma once

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
    T tau;

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
    DEVICE HOST INLINE
    const T& getDensity() {
      return density;
    }

    DEVICE HOST INLINE
    const MathVector<T, L::dimD>& getVelocity() {
      return velocity;
    }

    DEVICE HOST INLINE
    const MathVector<T, L::dimD>& getForce() {
      return force;
    }

    DEVICE HOST
    void calculateMoments(const T * haloDistributionPtr,
                                 const Position& iP) {
      { INSTRUMENT_OFF("Moment<T>::calculateMoments",4) }

      Moment_::calculateDensity(haloDistributionPtr, iP, density);
      Moment_::calculateVelocity(haloDistributionPtr, iP, density, velocity);

      velocity2 = velocity.norm2();
    }

    DEVICE HOST INLINE
    void setForce(T * localForceArray,
                  const Position& iP,
                  const unsigned int numberElements,
                  const Position& offset) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::setForce",4) }

      forcing.setForce(localForceArray, iP-L::halo(), numberElements, force);
      forcingScheme.setVariables(force, density, velocity);
    }

    DEVICE HOST
    inline const MathVector<T, L::dimD> getHydrodynamicVelocity() {
      return forcingScheme.calculateHydrodynamicVelocity(force, density, velocity);
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
    T alpha;
    const T beta;

  public:
    Collision(const T tau_in,
              const MathVector<T, 3>& amplitude_in,
              const MathVector<T, 3>& waveLength_in,
              const unsigned int kMin_in, const unsigned int kMax_in)
      : Base(tau_in, amplitude_in, waveLength_in, kMin_in, kMax_in)
      , alpha( (T) 2)
      , beta( (T) 1.0/(2.0 * tau_in))
  {}

    using Base::update;
    using Base::setForce;

    DEVICE HOST
    inline void calculateRelaxationTime(const T * haloDistributionNext_Ptr,
                                        const T * haloDistributionPrevious_Ptr,
                                        const Position& iP) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::calculate",4) }
    }


    DEVICE HOST
    inline void collideAndStream(T * haloDistributionNext_Ptr,
                                 const T * haloDistributionPrevious_Ptr,
                                 const Position& iP, const unsigned int iQ) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::calculate",4) }
      T equilibrium_iQ = Equilibrium_::calculate(Base::density, Base::velocity,
                                                 Base::velocity2, iQ);

      haloDistributionNext_Ptr[hSD::getIndex(iP, iQ)]
        = haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
        + Base::forcingScheme.calculateCollisionSource(Base::force, Base::density,
                                                       Base::velocity, Base::velocity2,
                                                       equilibrium_iQ, iQ)
        - 2*beta * (haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ],
                                                               iQ)] - equilibrium_iQ);
    }

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
    using Base::Collision;

    using Base::setForce;

    DEVICE HOST
    inline void calculateRelaxationTime(T * haloDistributionNext_Ptr,
                                        T * haloDistributionPrevious_Ptr,
                                        const Position& iP) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::calculate",4) }

      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        T equilibrium_iQ = Equilibrium_::calculate(Base::density, Base::velocity,
                                                   Base::velocity2, iQ);

        haloDistributionNext_Ptr[hSD::getIndex(iP, iQ)] =
          haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
          + Base::forcingScheme.calculateCollisionSource(Base::force, Base::density,
                                                         Base::velocity, Base::velocity2,
                                                         equilibrium_iQ, iQ);
        haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
          -= equilibrium_iQ;
      }

      calculateAlpha(haloDistributionNext_Ptr, haloDistributionPrevious_Ptr, iP);
      Base::tau = (T) 1.0/(Base::alpha*Base::beta);
    }

    DEVICE HOST
    inline void collideAndStream(T * haloDistributionNext_Ptr,
                                 const T * haloDistributionPrevious_Ptr,
                                 const Position& iP, const unsigned int iQ) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::calculate",4) }

      haloDistributionNext_Ptr[hSD::getIndex(iP, iQ)] -= (T) 1.0/Base::tau
        * haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)];
    }

    using Base::update;

    using Base::getHydrodynamicVelocity;
    using Base::getDensity;
    using Base::getVelocity;
    using Base::getForce;
    using Base::getAlpha;

  protected:
    using Base::forcingScheme;

    DEVICE HOST
    inline bool isDeviationSmall(const T * haloDistributionNext_Ptr,
                                 const T * haloDistributionPrevious_Ptr,
                                 const Position& iP, const T error) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::isDeviationSmall",6) }

      bool isDeviationSmallR = true;
      T deviation;

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        deviation = fabs(haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
                         /haloDistributionNext_Ptr[hSD::getIndex(iP, iQ)]);

        if(deviation > error) {
          isDeviationSmallR = false;
        }
      }

      return isDeviationSmallR;
    }

    HOST DEVICE
    T calculateAlphaMax(const T * haloDistributionNext_Ptr,
                        const T * haloDistributionPrevious_Ptr,
                        const Position& iP) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::calculateAlphaMax",6) }

      T alphaMaxR = 2.5;
      T alphaMaxTemp;

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        if(haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)] > 0) {
          alphaMaxTemp = fabs(haloDistributionNext_Ptr[hSD::getIndex(iP, iQ)]
             / haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]);

          if(alphaMaxTemp < alphaMaxR) {
            alphaMaxR = alphaMaxTemp;
          }
        }
      }

      return alphaMaxR;
    }

    HOST DEVICE
    inline T solveAlpha(const T * haloDistributionNext_Ptr,
                        const T * haloDistributionPrevious_Ptr,
                        const Position& iP, const T alphaMin, const T alphaMax) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::solveAlpha",6) }

      EntropicStepFunctor<T> entropicStepFunctor(haloDistributionNext_Ptr,
                                                 haloDistributionPrevious_Ptr, iP);
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

    HOST DEVICE
    inline void calculateAlpha(const T * haloDistributionNext_Ptr,
                               const T * haloDistributionPrevious_Ptr,
                               const Position& iP) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::calculateAlpha",5) }

      if(isDeviationSmall(haloDistributionNext_Ptr,
                          haloDistributionPrevious_Ptr,
                          iP, (T) 1.0e-3)) {
        Base::alpha = 2.0;
      }

      else {
        T alphaMax = calculateAlphaMax(haloDistributionNext_Ptr,
                                       haloDistributionPrevious_Ptr, iP);

        if(alphaMax < 2.) {
          Base::alpha = 0.95 * alphaMax;
        }

        else {
          T alphaMin = 1.;
          Base::alpha = solveAlpha(haloDistributionNext_Ptr,
                                   haloDistributionPrevious_Ptr, iP,
                                   alphaMin, alphaMax);
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
    using Base::Collision;

    using Base::update;
    using Base::setForce;

    using Base::collideAndStream;

    using Base::getHydrodynamicVelocity;
    using Base::getDensity;
    using Base::getVelocity;
    using Base::getForce;
    using Base::getAlpha;

  private:
    using Base::tau;
    using Base::alpha;

    using Base::isDeviationSmall;
    using Base::calculateAlphaMax;
    using Base::solveAlpha;

    HOST DEVICE
    inline T approximateAlpha(const T * haloDistributionNext_Ptr,
                              const T * haloDistributionPrevious_Ptr,
                              const Position& iP) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::Approached_ELBM>::approximateAlpha",6) }

      T a1 = (T) 0;
      T a2 = (T) 0;
      T a3 = (T) 0;
      T a4 = (T) 0;

      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        T temp = haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
          / haloDistributionNext_Ptr[hSD::getIndex(iP, iQ)];
        a1 += haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
          *temp;
        a2 += haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
          *temp*temp;
        a3 += haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
          *temp*temp*temp;
        a4 += haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
          *temp*temp*temp*temp;
      }

      a1 *= 1.0/2.0;
      a2 *= 1.0/6.0;
      a3 *= 1.0/12.0;
      a4 *= 1.0/20.0;

      T alphaR = 2 - 1/a1 * (4.0*a2 + 16.0*a2*a2/a1 - 8.0*a3
                             + 80.0*a2*a3/a1 - 80.0*a2*a2*a2/(a1*a1) - 16.0*a4);

      return alphaR;
    }

    HOST DEVICE
    inline void calculateAlpha(const T * haloDistributionNext_Ptr,
                               const T * haloDistributionPrevious_Ptr,
                               const Position& iP) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::Approached_ELBM>::calculateAlpha",5) }

        if(isRelativeDeviationSmall( (T) 1.0e-3)) {
          T alphaApproximated = approximateAlpha(haloDistributionNext_Ptr,
                                                 haloDistributionPrevious_Ptr, iP);
          Base::alpha = alphaApproximated;
        }
        else {
          T alphaMax = Base::calculateAlphaMax(haloDistributionNext_Ptr,
                                               haloDistributionPrevious_Ptr, iP);

          if(alphaMax < 2.) {
          Base::alpha = 0.95 * alphaMax;
          }

          else {
            T alphaMin = 1.;
            Base::alpha = Base::solveAlpha(haloDistributionNext_Ptr,
                                           haloDistributionPrevious_Ptr, iP,
                                           alphaMin, alphaMax);
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
    using Base::Collision;

    using Base::update;
    using Base::setForce;

    using Base::collideAndStream;

    using Base::getHydrodynamicVelocity;
    using Base::getDensity;
    using Base::getVelocity;
    using Base::getForce;
    using Base::getAlpha;

  private:
    DEVICE HOST
    inline void calculateAlpha(const T * haloDistributionNext_Ptr,
                               const T * haloDistributionPrevious_Ptr,
                               const Position& iP) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::ForcedNR_ELBM>::calculateAlpha",5) }

      T alphaMax = Base::calculateAlphaMax(haloDistributionNext_Ptr,
                                           haloDistributionPrevious_Ptr, iP);

      if(alphaMax < 2.) {
        Base::alpha = 0.95 * alphaMax;
      }

      else {
        T alphaMin = 1.;
        Base::alpha = Base::solveAlpha(haloDistributionNext_Ptr,
                                       haloDistributionPrevious_Ptr, iP,
                                       alphaMin, alphaMax);
      }
    }
  };

  template <class T>
  class Collision<T, CollisionType::ForcedBNR_ELBM>
    : public Collision<T, CollisionType::ForcedNR_ELBM> {
  private:
    using Base = Collision<T, CollisionType::ForcedNR_ELBM>;

  public:
    using Base::Collision;

    using Base::update;
    using Base::setForce;

    using Base::collideAndStream;

    using Base::getHydrodynamicVelocity;
    using Base::getDensity;
    using Base::getVelocity;
    using Base::getForce;
    using Base::getAlpha;

  private:
    DEVICE HOST
    inline T solveAlpha(const T * haloDistributionNext_Ptr,
                        const T * haloDistributionPrevious_Ptr,
                        const Position& iP, const T alphaMin, const T alphaMax) {
      { INSTRUMENT_OFF("Collision<T, CollisionType::ForcedBNR_ELBM>::solveAlpha",6) }

      EntropicStepFunctor<T> entropicStepFunctor(haloDistributionNext_Ptr,
                                                 haloDistributionPrevious_Ptr, iP);
      const T tolerance = 1e-5;
      const int iterationMax = 20;
      T alphaR = Base::alpha;

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
