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
  class Collision<T, CollisionType::BGK> {
  protected:
    T alpha;
    const T beta;

    Force_ forcing;
    ForcingScheme_ forcingScheme;

    T density;
    MathVector<T, L::dimD> velocity;
    T velocity2;
    MathVector<T, L::dimD> force;
    T entropy;

  public:
    Collision(const T tau_in,
              const MathVector<T, 3>& amplitude_in,
              const MathVector<T, 3>& waveLength_in,
              const unsigned int kMin_in, const unsigned int kMax_in)
      : alpha( (T) 2)
      , beta( (T) 1.0/(2.0 * tau_in))
      , forcing(amplitude_in, waveLength_in, kMin_in, kMax_in)
      , forcingScheme(tau_in)
      , density()
      , velocity{{0}}
      , velocity2()
      , force{{0}}
      , entropy()
    {}

    LBM_DEVICE LBM_HOST LBM_INLINE
    const T& getDensity() {
      return density;
    }

    LBM_DEVICE LBM_HOST LBM_INLINE
    const MathVector<T, L::dimD>& getVelocity() {
      return velocity;
    }

    LBM_DEVICE LBM_HOST LBM_INLINE
    const MathVector<T, L::dimD>& getForce() {
      return force;
    }

    LBM_DEVICE LBM_HOST
    inline T getAlpha() {
      return alpha;
    }


    LBM_DEVICE LBM_HOST
    void calculateMoments(const T * haloDistributionPtr,
                                 const Position& iP) {
      { LBM_INSTRUMENT_OFF("Moment<T>::calculateMoments",4) }

      Moment_::calculateDensity(haloDistributionPtr, iP, density);
      Moment_::calculateVelocity(haloDistributionPtr, iP, density, velocity);

      velocity2 = velocity.norm2();
    }

    LBM_DEVICE LBM_HOST LBM_INLINE
    void setForce(T * localForceArray,
                  const Position& iP,
                  const unsigned int numberElements,
                  const Position& offset) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::setForce",4) }

      forcing.setForce(localForceArray, iP-L::halo(), numberElements, force);
      forcingScheme.setVariables(force, density, velocity);
    }

    LBM_DEVICE LBM_HOST
    inline const MathVector<T, L::dimD> getHydrodynamicVelocity() {
      return forcingScheme.calculateHydrodynamicVelocity(force, density, velocity);
    }

    LBM_DEVICE LBM_HOST
    inline void update(const unsigned int iteration) {
      forcing.update(iteration);
    }

    LBM_DEVICE LBM_HOST
    inline void calculateRelaxationTime(T * haloDistributionPrevious_Ptr,
                                        const Position& iP,
                                        MathVector<T, L::dimQ>& f_Forced,
                                        MathVector<T, L::dimQ>& f_NonEq) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::calculate",4) }


      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        T equilibrium_iQ = Equilibrium_::calculate(density, velocity, velocity2, iQ);

        f_NonEq[iQ] =
          haloDistributionPrevious_Ptr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)]
          - equilibrium_iQ;

        f_Forced[iQ] =
          f_NonEq[iQ] + equilibrium_iQ
          + forcingScheme.calculateCollisionSource(force, density, velocity, velocity2,
                                                   equilibrium_iQ, iQ);
      }
    }

    LBM_DEVICE LBM_HOST
    inline void collideAndStream(T * haloDistributionNext_Ptr,
                                 const Position& iP, const unsigned int iQ,
                                 const MathVector<T, L::dimQ>& f_Forced,
                                 const MathVector<T, L::dimQ>& f_NonEq) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::calculate",4) }

      haloDistributionNext_Ptr[hSD::getIndex(iP, iQ)] = f_Forced[iQ]
        - alpha*beta * f_NonEq[iQ];
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

    LBM_DEVICE LBM_HOST
    inline void calculateRelaxationTime(T * haloDistributionPrevious_Ptr,
                                        const Position& iP,
                                        MathVector<T, L::dimQ>& f_Forced,
                                        MathVector<T, L::dimQ>& f_NonEq) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::calculate",4) }

      Base::calculateRelaxationTime(haloDistributionPrevious_Ptr, iP, f_Forced, f_NonEq);

      calculateAlpha(f_Forced, f_NonEq);
    }

    using Base::collideAndStream;
    using Base::update;

    using Base::getHydrodynamicVelocity;
    using Base::getDensity;
    using Base::getVelocity;
    using Base::getForce;
    using Base::getAlpha;

  protected:
    LBM_DEVICE LBM_HOST
    inline bool isDeviationSmall(const MathVector<T, L::dimQ>& f_Forced,
                                 const MathVector<T, L::dimQ>& f_NonEq,
                                 const T error) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::isDeviationSmall",6) }

      bool isDeviationSmallR = true;
      T deviation;

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        deviation = fabs(f_NonEq[iQ] / f_Forced[iQ]);

        if(deviation > error) {
          isDeviationSmallR = false;
        }
      }

      return isDeviationSmallR;
    }

    LBM_HOST LBM_DEVICE
    T calculateAlphaMax(const MathVector<T, L::dimQ>& f_Forced,
                        const MathVector<T, L::dimQ>& f_NonEq) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::calculateAlphaMax",6) }

      T alphaMaxR = 2.5;
      T alphaMaxTemp;

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        if(f_NonEq[iQ] > 0) {
          alphaMaxTemp = fabs(f_Forced[iQ] / f_NonEq[iQ]);

          if(alphaMaxTemp < alphaMaxR) {
            alphaMaxR = alphaMaxTemp;
          }
        }
      }

      return alphaMaxR;
    }

    LBM_HOST LBM_DEVICE
    inline T solveAlpha(const MathVector<T, L::dimQ>& f_Forced,
                        const MathVector<T, L::dimQ>& f_NonEq,
                        const T alphaMin, const T alphaMax) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::solveAlpha",6) }

      EntropicStepFunctor<T> entropicStepFunctor(f_Forced, f_NonEq);

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

    LBM_HOST LBM_DEVICE
    inline void calculateAlpha(const MathVector<T, L::dimQ>& f_Forced,
                               const MathVector<T, L::dimQ>& f_NonEq) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::calculateAlpha",5) }

      if(isDeviationSmall(f_Forced, f_NonEq, (T) 1.0e-3)) {
        Base::alpha = 2.0;
      }

      else {
        T alphaMax = calculateAlphaMax(f_Forced, f_NonEq);

        if(alphaMax < 2.) {
          Base::alpha = 0.95 * alphaMax;
        }

        else {
          T alphaMin = 1.;
          Base::alpha = solveAlpha(f_Forced, f_NonEq, alphaMin, alphaMax);

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
    LBM_HOST LBM_DEVICE
    inline T approximateAlpha(const MathVector<T, L::dimQ>& f_Forced,
                              const MathVector<T, L::dimQ>& f_NonEq) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::Approached_ELBM>::approximateAlpha",6) }

      T a1 = (T) 0;
      T a2 = (T) 0;
      T a3 = (T) 0;
      T a4 = (T) 0;

      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        T temp = Base:: f_NonEq[iQ] / f_Forced[iQ];
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

    LBM_HOST LBM_DEVICE
    inline void calculateAlpha(const MathVector<T, L::dimQ>& f_Forced,
                               const MathVector<T, L::dimQ>& f_NonEq) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::Approached_ELBM>::calculateAlpha",5) }

      if(Base::isDeviationSmall(f_Forced, f_NonEq, (T) 1.0e-5)) {
        T alphaApproximated = approximateAlpha(f_Forced, f_NonEq);
        Base::alpha = alphaApproximated;
        }

      else {
        T alphaMax = Base::calculateAlphaMax(f_Forced, f_NonEq);

        if(alphaMax < 2.) {
          Base::alpha = 0.95 * alphaMax;
        }

        else {
          T alphaMin = 1.;
          Base::alpha = Base::solveAlpha(f_Forced, f_NonEq, alphaMin, alphaMax);
        }
      }
    }
  };


  template <class T>
  class Collision<T, CollisionType::Essentially1_ELBM>
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
    LBM_HOST LBM_DEVICE
    inline void calculateAlpha(const MathVector<T, L::dimQ>& f_Forced,
                               const MathVector<T, L::dimQ>& f_NonEq) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::Essentially_ELBM>::calculateAlpha",5) }

      T term1 = (T) 0;
      T term2 = (T) 0;
      T term3 = (T) 0;

      T x_iQ;
      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        x_iQ = - f_Forced[iQ] / f_NonEq[iQ];

        term1 += f_NonEq[iQ]*x_iQ*x_iQ;

        if(x_iQ < 0) term2 += f_NonEq[iQ]*x_iQ*x_iQ*x_iQ;

        term3 += f_NonEq[iQ]*2.0*x_iQ*x_iQ/(2.0+x_iQ);
      }

      Base::alpha = (term1 - sqrt(term1*term1 - 8.0*term2*term3))/(2.0*term3);
    }
  };


  template <class T>
  class Collision<T, CollisionType::Essentially2_ELBM>
    : public Collision<T, CollisionType::Essentially1_ELBM> {
  private:
    using Base = Collision<T, CollisionType::Essentially1_ELBM>;

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
    LBM_HOST LBM_DEVICE
    inline void calculateAlpha(const MathVector<T, L::dimQ>& f_Forced,
                               const MathVector<T, L::dimQ>& f_NonEq) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::Essentially_ELBM>::calculateAlpha",5) }

      Base::calculateAlpha();

      T alpha1 = Base::alpha;

      T A = (T) 0;
      T B = (T) 0;
      T C = (T) 0;

      T x_iQ;
      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        x_iQ = - f_Forced[iQ] / f_NonEq[iQ];

        if(x_iQ < 0) {
          A -= f_NonEq[iQ]*x_iQ*x_iQ*x_iQ/6.0;
        }
        else {
          B -= Base::beta*Base::beta*f_NonEq[iQ]
            *2.0*alpha1*x_iQ*x_iQ*x_iQ/15.0
            *(2.0/(4.0+alpha1*x_iQ)+1.0/(4.0+2.0*alpha1*x_iQ)+2.0/(4.0+3.0*alpha1*x_iQ));
        }

        B += f_NonEq[iQ]*x_iQ*x_iQ/2.0;

        C += f_NonEq[iQ]
          *(x_iQ*x_iQ*(60.0*(1+x_iQ)+11.0*x_iQ*x_iQ))/(60.0+x_iQ*(90.0+x_iQ*(36.0+3.0*x_iQ)));
      }

      A *= Base::beta*Base::beta;

      T alpha2p = (-B + sqrt(B*B-4.0*A*C))/(2.0*A);
      T alpha2m = (-B - sqrt(B*B-4.0*A*C))/(2.0*A);

      T alpha2 = alpha2p;

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        x_iQ = - f_Forced[iQ] / f_NonEq[iQ];

        if(x_iQ < 0) {
          A += f_NonEq[iQ]
            *(alpha2*Base::beta*x_iQ*x_iQ*x_iQ*x_iQ*(-1.0/12.0+alpha2*Base::beta*x_iQ(1.0/20.0-alpha2*Base::beta*x_iQ*1.0/5.0)));
        }
      }

      T alphap = (-B + sqrt(B*B-4.0*A*C))/(2.0*A);
      T alpham = (-B - sqrt(B*B-4.0*A*C))/(2.0*A);

      Base::alpha = (T) 0;
      if(alphap > alpha1 && alphap < alpha2) Base::alpha = alphap;
      else if(alpham > alpha1 && alpham < alpha2) Base::alpha = alpham;
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
    LBM_DEVICE LBM_HOST
    inline void calculateAlpha(const MathVector<T, L::dimQ>& f_Forced,
                               const MathVector<T, L::dimQ>& f_NonEq) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ForcedNR_ELBM>::calculateAlpha",5) }

      T alphaMax = Base::calculateAlphaMax(f_Forced, f_NonEq);

      if(alphaMax < 2.) {
        Base::alpha = 0.95 * alphaMax;
      }

      else {
        T alphaMin = 1.;
        Base::alpha = Base::solveAlpha(f_Forced, f_NonEq, alphaMin, alphaMax);
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
    LBM_DEVICE LBM_HOST
    inline T solveAlpha(const MathVector<T, L::dimQ>& f_Forced,
                        const MathVector<T, L::dimQ>& f_NonEq,
                        const T alphaMin, const T alphaMax) {
      { LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ForcedBNR_ELBM>::solveAlpha",6) }

      EntropicStepFunctor<T> entropicStepFunctor(f_Forced, f_NonEq);

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
