#pragma once

#include "Commons.h"
#include "Domain.h"
#include "EntropicStep.h"
#include "Equilibrium.h"
#include "Force.h"
#include "ForcingScheme.h"
#include "Helpers.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Moment.h"
#include "Options.h"

namespace lbm {

  template <class T, CollisionType collisionType, Architecture architecture>
class Collision {};

template <class T, Architecture architecture>
  class Collision<T, CollisionType::GenericSRT, architecture> {
 protected:
  T tau;

  Force_<architecture> forcing;
  ForcingScheme_ forcingScheme;
  FieldList<T, architecture>& fieldList;

  T density;
  MathVector<T, L::dimD> velocity;
  T velocity2;
  MathVector<T, L::dimD> force;
  T entropy;

  Collision(const T tau_in,
            FieldList<T, architecture>& fieldList_in,
            const MathVector<T, 3>& amplitude_in,
            const MathVector<T, 3>& waveLength_in,
            const unsigned int kMin_in,
            const unsigned int kMax_in)
    : tau(tau_in)
    , fieldList(fieldList_in)
    , forcing(gFD::offset(MPIInit::rank), amplitude_in, waveLength_in, kMin_in, kMax_in)
    , forcingScheme(tau_in)
    , density()
    , velocity{{0}}
    , velocity2()
    , force{{0}}
    , entropy()
  {
    if (writeForce) {
      forcing.setForceArray(fieldList.force.getData(FFTWInit::numberElements),
                            fieldList);
    }

  }

 public:
  LBM_DEVICE LBM_HOST LBM_INLINE
  const T& getDensity() { return density; }

  LBM_DEVICE LBM_HOST LBM_INLINE
  const MathVector<T, L::dimD>& getVelocity() {
    return velocity;
  }

  LBM_DEVICE LBM_HOST LBM_INLINE
  const MathVector<T, L::dimD>& getForce() {
    return force;
  }

  LBM_DEVICE LBM_HOST void
  calculateMoments(const T* haloDistributionPreviousPtr, const Position& iP) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateMoments", 4)

    Moment_::calculateDensity(haloDistributionPreviousPtr, iP, density);
    Moment_::calculateVelocity(haloDistributionPreviousPtr, iP, density, velocity);
    velocity2 = velocity.norm2();
  }

  LBM_DEVICE LBM_HOST LBM_INLINE
  void setForce(T * forcePtr, const Position& iP,
                const Position& offset, const unsigned int numberElements) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::setForce", 4)

    forcing.setForce(forcePtr, iP - L::halo(), force, numberElements);
    forcingScheme.setVariables(force, density, velocity);
  }

  LBM_DEVICE LBM_HOST inline
  const MathVector<T, L::dimD> getHydrodynamicVelocity() {
    return forcingScheme.calculateHydrodynamicVelocity(force, density,
                                                       velocity);
  }

  LBM_HOST inline
  void update(const unsigned int iteration, const unsigned int numberElements) {
    forcing.update(fieldList.force.getData(numberElements),
                   fieldList, iteration);
  }
};

template <class T, Architecture architecture>
  class Collision<T, CollisionType::BGK, architecture>
  : public Collision<T, CollisionType::GenericSRT, architecture> {
 private:
  using Base = Collision<T, CollisionType::GenericSRT, architecture>;

 protected:
  T alpha;
  const T beta;

 public:
  Collision(const T tau_in,
            FieldList<T, architecture>& fieldList_in,
            const MathVector<T, 3>& amplitude_in,
            const MathVector<T, 3>& waveLength_in,
            const unsigned int kMin_in,
            const unsigned int kMax_in)
    : Base(tau_in, fieldList_in, amplitude_in, waveLength_in, kMin_in, kMax_in)
    , alpha((T)2)
    , beta((T)1.0 / (2.0 * tau_in))
  {}

  using Base::setForce;
  using Base::update;

  LBM_DEVICE LBM_HOST LBM_INLINE
  void calculateRelaxationTime(const T* haloDistributionNextPtr, const T* haloDistributionPreviousPtr,
                               const Position& iP, const T alphaGuess) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::calculate", 4)
  }

  LBM_DEVICE LBM_HOST inline
    void collideAndStream(T * haloDistributionNextPtr, const T* haloDistributionPreviousPtr,
                                                   const Position& iP,
                                                   const unsigned int iQ) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::GenericSRT>::calculate", 4)

    T equilibrium_iQ =
      Equilibrium_::calculate(Base::density, Base::velocity, Base::velocity2,
                                iQ);

    haloDistributionNextPtr[hSD::getIndex(iP, iQ)] =
      (1.-2.*beta)
      * haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
      + 2.*beta * equilibrium_iQ
      + Base::forcingScheme.calculateCollisionSource(Base::force, Base::density,
                                                     Base::velocity, Base::velocity2,
                                                     equilibrium_iQ, iQ);
  }

  LBM_DEVICE LBM_HOST
  void calculateTs(const T* haloDistributionPreviousPtr, const T* haloDistributionNextPtr,
                   const Position& iP) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateMoments", 4)

  }

  LBM_DEVICE LBM_HOST inline T getT2() { return 0; }
  LBM_DEVICE LBM_HOST inline T getT3() { return 0; }
  LBM_DEVICE LBM_HOST inline T getT4() { return 0; }

  using Base::getDensity;
  using Base::getForce;
  using Base::getHydrodynamicVelocity;
  using Base::getVelocity;

  LBM_DEVICE LBM_HOST inline T getAlpha() { return alpha; }
};

template <class T, Architecture architecture>
  class Collision<T, CollisionType::ELBM, architecture>
  : public Collision<T, CollisionType::BGK, architecture> {
 private:
  using Base = Collision<T, CollisionType::BGK, architecture>;

protected:
  T T2, T3, T4;

 public:
  Collision(const T tau_in,
            FieldList<T, architecture>& fieldList_in,
            const MathVector<T, 3>& amplitude_in,
            const MathVector<T, 3>& waveLength_in,
            const unsigned int kMin_in,
            const unsigned int kMax_in)
    : Base(tau_in, fieldList_in, amplitude_in, waveLength_in, kMin_in, kMax_in)
    , T2((T)0), T3((T)0), T4((T)0)
  {}

  using Base::setForce;

  LBM_DEVICE LBM_HOST void
  calculateTs(const T* haloDistributionPreviousPtr, const T* haloDistributionNextPtr,
              const Position& iP) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateMoments", 4)

    if(writeT) {
      Moment_::calculateT(haloDistributionPreviousPtr, haloDistributionNextPtr, iP,
                          T2, T3, T4);
    }
  }


  LBM_DEVICE LBM_HOST inline
  void calculateRelaxationTime(T * haloDistributionNextPtr, T * haloDistributionPreviousPtr,
                               const Position& iP, const T alphaGuess) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::calculate", 4)

    for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
      haloDistributionNextPtr[hSD::getIndex(iP, iQ)] =
        haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
        - Equilibrium_::calculate(Base::density, Base::velocity, Base::velocity2, iQ);
    }

    Base::alpha = alphaGuess;
    calculateAlpha(haloDistributionNextPtr, haloDistributionPreviousPtr, iP);
    Base::tau = (T)1.0 / (Base::alpha * Base::beta);
  }

  LBM_DEVICE LBM_HOST inline
  void collideAndStream(T* haloDistributionNextPtr, const T* haloDistributionPreviousPtr,
                        const Position& iP, const unsigned int iQ) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::calculate", 4)

    T equilibrium_iQ =
      haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
      - haloDistributionNextPtr[hSD::getIndex(iP, iQ)];

    haloDistributionNextPtr[hSD::getIndex(iP, iQ)] =
      haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
      - (T)1.0 / Base::tau * haloDistributionNextPtr[hSD::getIndex(iP, iQ)]
      + Base::forcingScheme.calculateCollisionSource(Base::force, Base::density,
                                                     Base::velocity, Base::velocity2,
                                                     equilibrium_iQ, iQ);
  }

  using Base::update;

  using Base::getAlpha;
  using Base::getDensity;
  using Base::getForce;
  using Base::getHydrodynamicVelocity;
  using Base::getVelocity;

  LBM_DEVICE LBM_HOST inline T getT2() { return T2; }
  LBM_DEVICE LBM_HOST inline T getT3() { return T3; }
  LBM_DEVICE LBM_HOST inline T getT4() { return T4; }

 protected:
  using Base::forcingScheme;

  LBM_DEVICE LBM_HOST inline
  bool isDeviationSmall(const T* haloDistributionNextPtr, const T* haloDistributionPreviousPtr,
                        const Position& iP, const T error) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::isDeviationSmall", 6)

    bool isDeviationSmallR = true;
    T deviation;

    for (unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
      deviation =
        fabs(haloDistributionNextPtr[hSD::getIndex(iP, iQ)])
             / haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)];

      if (deviation > error) {
        isDeviationSmallR = false;
      }
    }

    return isDeviationSmallR;
  }

  LBM_HOST LBM_DEVICE
  T calculateAlphaMax(const T* haloDistributionNextPtr, const T* haloDistributionPreviousPtr,
                      const Position& iP) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::calculateAlphaMax", 6)

    T alphaMaxR = 2.5;
    T alphaMaxTemp;

    for (unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
      if (haloDistributionNextPtr[hSD::getIndex(iP, iQ)] > 0) {
        alphaMaxTemp =
          fabs(haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)])
               / haloDistributionNextPtr[hSD::getIndex(iP, iQ)];

        if (alphaMaxTemp < alphaMaxR) {
          alphaMaxR = alphaMaxTemp;
        }
      }
    }

    return alphaMaxR;
  }

  LBM_HOST LBM_DEVICE inline
  T solveAlpha(const T* haloDistributionNextPtr, const T* haloDistributionPreviousPtr,
               const Position& iP, const T alphaMin, const T alphaMax) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::solveAlpha",6)

    EntropicStepFunctor<T>
        entropicStepFunctor(haloDistributionNextPtr,
                            haloDistributionPreviousPtr, iP);
    const T tolerance = 1e-8;
    const int iterationMax = 50;
    T alphaR = Base::alpha;

    bool hasConverged =
        NewtonRaphsonSolver(entropicStepFunctor, tolerance, iterationMax,
                            alphaR, alphaMin, alphaMax);

    if(!hasConverged) {
      return 2.0;
    }

    return alphaR;
  }

  LBM_HOST LBM_DEVICE inline
  void calculateAlpha(const T* haloDistributionNextPtr, const T* haloDistributionPreviousPtr,
                      const Position& iP) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ELBM>::calculateAlpha", 5)

    if(isDeviationSmall(haloDistributionNextPtr, haloDistributionPreviousPtr,
                         iP, (T)1.0e-3)) {
      Base::alpha = 2.0;
    }

    else {
      T alphaMax = calculateAlphaMax(haloDistributionNextPtr,
                                     haloDistributionPreviousPtr, iP);

      if(alphaMax < 2.) {
        Base::alpha = 0.95 * alphaMax;
      }

      else {
        T alphaMin = 1.;
        Base::alpha = solveAlpha(haloDistributionNextPtr, haloDistributionPreviousPtr,
                                 iP, alphaMin, alphaMax);
      }
    }
  }
};

template <class T, Architecture architecture>
  class Collision<T, CollisionType::Approached_ELBM, architecture>
  : public Collision<T, CollisionType::ELBM, architecture> {
 private:
  using Base = Collision<T, CollisionType::ELBM, architecture>;

 public:
  using Base::Collision;

  using Base::setForce;
  using Base::update;

  using Base::collideAndStream;

  using Base::getAlpha;
  using Base::getDensity;
  using Base::getForce;
  using Base::getHydrodynamicVelocity;
  using Base::getVelocity;

 private:
  LBM_HOST LBM_DEVICE inline T approximateAlpha(
      const T* haloDistributionNextPtr,
      const T* haloDistributionPreviousPtr,
      const Position& iP) {
    LBM_INSTRUMENT_OFF( "Collision<T, CollisionType::Approached_ELBM>::approximateAlpha", 6)

    T a1 = (T)0;
    T a2 = (T)0;
    T a3 = (T)0;
    T a4 = (T)0;

    for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
      T temp = haloDistributionNextPtr[hSD::getIndex(iP, iQ)]
        / haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)];

      a1 += haloDistributionNextPtr[hSD::getIndex(iP, iQ)] * temp;
      a2 += haloDistributionNextPtr[hSD::getIndex(iP, iQ)] * temp * temp;
      a3 += haloDistributionNextPtr[hSD::getIndex(iP, iQ)] * temp * temp * temp;
      a4 += haloDistributionNextPtr[hSD::getIndex(iP, iQ)] * temp * temp * temp * temp;
    }

    a1 *= 1.0 / 2.0;
    a2 *= 1.0 / 6.0;
    a3 *= 1.0 / 12.0;
    a4 *= 1.0 / 20.0;

    T alphaR = 2 - 1 / a1 *
                       (4.0 * a2 + 16.0 * a2 * a2 / a1 - 8.0 * a3 +
                        80.0 * a2 * a3 / a1 - 80.0 * a2 * a2 * a2 / (a1 * a1) -
                        16.0 * a4);

    return alphaR;
  }

  LBM_HOST LBM_DEVICE inline void calculateAlpha(
      const T* haloDistributionNextPtr,
      const T* haloDistributionPreviousPtr,
      const Position& iP) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::Approached_ELBM>::calculateAlpha", 5)

    if (isRelativeDeviationSmall((T)1.0e-3)) {
      T alphaApproximated = approximateAlpha(haloDistributionNextPtr,
                                             haloDistributionPreviousPtr, iP);
      Base::alpha = alphaApproximated;
    } else {
      T alphaMax = Base::calculateAlphaMax(haloDistributionNextPtr,
                                           haloDistributionPreviousPtr, iP);

      if (alphaMax < 2.) {
        Base::alpha = 0.95 * alphaMax;
      }

      else {
        T alphaMin = 1.;
        Base::alpha = Base::solveAlpha(haloDistributionNextPtr,
                                       haloDistributionPreviousPtr, iP,
                                       alphaMin, alphaMax);
      }
    }
  }
};

template <class T, Architecture architecture>
  class Collision<T, CollisionType::Malaspinas_ELBM, architecture>
  : public Collision<T, CollisionType::ELBM, architecture> {
 private:
  using Base = Collision<T, CollisionType::ELBM, architecture>;

 public:
  using Base::Collision;

  using Base::setForce;
  using Base::update;

  using Base::collideAndStream;

  using Base::getAlpha;
  using Base::getDensity;
  using Base::getForce;
  using Base::getHydrodynamicVelocity;
  using Base::getVelocity;

 private:

  LBM_HOST LBM_DEVICE inline void calculateAlpha(
      const T* haloDistributionNextPtr,
      const T* haloDistributionPreviousPtr,
      const Position& iP) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::Malaspinas_ELBM>::calculateAlpha", 5)

    MathVector<T, 3> dissipativeTensor_diag{{0}};
    MathVector<T, 3> dissipativeTensor_sym{{0}};

    for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
      dissipativeTensor_diag[d::X] += L::celerity()[iQ][d::X] * L::celerity()[iQ][d::X]
        * haloDistributionNextPtr[hSD::getIndex(iP, iQ)];

      dissipativeTensor_diag[d::Y] += L::celerity()[iQ][d::Y] * L::celerity()[iQ][d::Y]
        * haloDistributionNextPtr[hSD::getIndex(iP, iQ)];

      dissipativeTensor_diag[d::Z] += L::celerity()[iQ][d::Z] * L::celerity()[iQ][d::Z]
        * haloDistributionNextPtr[hSD::getIndex(iP, iQ)];

      dissipativeTensor_sym[d::X] += L::celerity()[iQ][d::X] * L::celerity()[iQ][d::Y]
        * haloDistributionNextPtr[hSD::getIndex(iP, iQ)];

      dissipativeTensor_sym[d::Y] += L::celerity()[iQ][d::X] * L::celerity()[iQ][d::Z]
        * haloDistributionNextPtr[hSD::getIndex(iP, iQ)];

      dissipativeTensor_sym[d::Z] += L::celerity()[iQ][d::Y] * L::celerity()[iQ][d::Z]
        * haloDistributionNextPtr[hSD::getIndex(iP, iQ)];
    }

    T traceDissipativeTensor2 =
      dissipativeTensor_diag[d::X] * dissipativeTensor_diag[d::X]
      + dissipativeTensor_diag[d::Y] * dissipativeTensor_diag[d::Y]
      + dissipativeTensor_diag[d::Z] * dissipativeTensor_diag[d::Z]
      + 2. * dissipativeTensor_sym[d::X] * dissipativeTensor_sym[d::X]
      + 2. * dissipativeTensor_sym[d::Y] * dissipativeTensor_sym[d::Y]
      + 2. * dissipativeTensor_sym[d::Z] * dissipativeTensor_sym[d::Z];

    T traceDissipativeTensor3 =
      dissipativeTensor_diag[d::X]*dissipativeTensor_diag[d::X]*dissipativeTensor_diag[d::X]
      + dissipativeTensor_diag[d::Y]*dissipativeTensor_diag[d::Y]*dissipativeTensor_diag[d::Y]
      + dissipativeTensor_diag[d::Z]*dissipativeTensor_diag[d::Z]*dissipativeTensor_diag[d::Z]
      + dissipativeTensor_diag[d::X]
      *(3.*(dissipativeTensor_sym[d::X]*dissipativeTensor_sym[d::X]
            + dissipativeTensor_sym[d::Y]*dissipativeTensor_sym[d::Y]))
      + dissipativeTensor_diag[d::Y]
      *(3.*(dissipativeTensor_sym[d::X]*dissipativeTensor_sym[d::X]
            + dissipativeTensor_sym[d::Z]*dissipativeTensor_sym[d::Z]))
      + dissipativeTensor_diag[d::Z]
      *(3.*(dissipativeTensor_sym[d::Y]*dissipativeTensor_sym[d::Y]
            + dissipativeTensor_sym[d::Z]*dissipativeTensor_sym[d::Z]))
      + 6.*(dissipativeTensor_sym[d::X]*dissipativeTensor_sym[d::Y]
            * dissipativeTensor_sym[d::Z]);

    Base::alpha = 2. - 2./(3.*Base::density*L::cs2)
      * traceDissipativeTensor3/traceDissipativeTensor2;
  }
};



template <class T, Architecture architecture>
  class Collision<T, CollisionType::Essentially1_ELBM, architecture>
  : public Collision<T, CollisionType::ELBM, architecture> {
 private:
  using Base = Collision<T, CollisionType::ELBM, architecture>;

 public:
  using Base::Collision;

  using Base::setForce;
  using Base::update;

  using Base::collideAndStream;

  using Base::getAlpha;
  using Base::getDensity;
  using Base::getForce;
  using Base::getHydrodynamicVelocity;
  using Base::getVelocity;

 private:
  LBM_HOST LBM_DEVICE inline void calculateAlpha(
      const T* haloDistributionNextPtr,
      const T* haloDistributionPreviousPtr,
      const Position& iP) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::Essentially_ELBM>::calculateAlpha", 5)

    T term1 = (T)0;
    T term2 = (T)0;
    T term3 = (T)0;

    T x_iQ;
    for (unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
      x_iQ = - haloDistributionNextPtr[hSD::getIndex(iP, iQ)]
        / haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)];

      term1 += haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
        * x_iQ * x_iQ;

      if (x_iQ < 0)
        term2 += haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
          * x_iQ * x_iQ * x_iQ;

      term3 += haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
        * 2.0 * x_iQ * x_iQ / (2.0 + x_iQ);
    }

    Base::alpha =
        (term1 - sqrt(term1 * term1 - 8.0 * term2 * term3)) / (2.0 * term3);
  }
};

template <class T, Architecture architecture>
  class Collision<T, CollisionType::Essentially2_ELBM, architecture>
  : public Collision<T, CollisionType::Essentially1_ELBM, architecture> {
 private:
  using Base = Collision<T, CollisionType::Essentially1_ELBM, architecture>;

 public:
  using Base::Collision;

  using Base::setForce;
  using Base::update;

  using Base::collideAndStream;

  using Base::getAlpha;
  using Base::getDensity;
  using Base::getForce;
  using Base::getHydrodynamicVelocity;
  using Base::getVelocity;

 private:
  LBM_HOST LBM_DEVICE inline void calculateAlpha(
      const T* haloDistributionNextPtr,
      const T* haloDistributionPreviousPtr,
      const Position& iP) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::Essentially_ELBM>::calculateAlpha", 5)

    Base::calculateAlpha(haloDistributionNextPtr, haloDistributionPreviousPtr, iP);

    T alpha1 = Base::alpha;

    T A = (T)0;
    T B = (T)0;
    T C = (T)0;

    T x_iQ;
    for (unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
      x_iQ = - haloDistributionNextPtr[hSD::getIndex(iP, iQ)]
        / haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)];

      if (x_iQ < 0) {
        A -= haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
          * x_iQ * x_iQ * x_iQ / 6.0;
      } else {
        B -= Base::beta * Base::beta
          * haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
          * 2.0 * alpha1 * x_iQ * x_iQ * x_iQ / 15.0 *
          (2.0 / (4.0 + alpha1 * x_iQ) + 1.0 / (4.0 + 2.0 * alpha1 * x_iQ) +
           2.0 / (4.0 + 3.0 * alpha1 * x_iQ));
      }

      B += haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
        * x_iQ * x_iQ / 2.0;

      C += haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
        * (x_iQ * x_iQ * (60.0 * (1 + x_iQ) + 11.0 * x_iQ * x_iQ))
        / (60.0 + x_iQ * (90.0 + x_iQ * (36.0 + 3.0 * x_iQ)));
    }

    A *= Base::beta * Base::beta;

    T alpha2p = (-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
    T alpha2m = (-B - sqrt(B * B - 4.0 * A * C)) / (2.0 * A);

    T alpha2 = alpha2p;

    for (unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
      x_iQ = - haloDistributionNextPtr[hSD::getIndex(iP, iQ)]
        / haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)];

      if (x_iQ < 0) {
        A += haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
          * (alpha2 * Base::beta * x_iQ * x_iQ * x_iQ * x_iQ *
             (-1.0 / 12.0 + alpha2 * Base::beta
              * x_iQ(1.0 / 20.0 - alpha2 * Base::beta * x_iQ * 1.0 / 5.0)));
      }
    }

    T alphap = (-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
    T alpham = (-B - sqrt(B * B - 4.0 * A * C)) / (2.0 * A);

    Base::alpha = (T)0;
    if (alphap > alpha1 && alphap < alpha2)
      Base::alpha = alphap;
    else if (alpham > alpha1 && alpham < alpha2)
      Base::alpha = alpham;
  }
};

template <class T, Architecture architecture>
  class Collision<T, CollisionType::ForcedNR_ELBM, architecture>
  : public Collision<T, CollisionType::ELBM, architecture> {
 private:
  using Base = Collision<T, CollisionType::ELBM, architecture>;

 public:
  using Base::Collision;

  using Base::setForce;
  using Base::update;

  using Base::collideAndStream;

  using Base::getAlpha;
  using Base::getDensity;
  using Base::getForce;
  using Base::getHydrodynamicVelocity;
  using Base::getVelocity;

 protected:
  LBM_DEVICE LBM_HOST inline
  void calculateAlpha(const T* haloDistributionNextPtr, const T* haloDistributionPreviousPtr,
                      const Position& iP) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ForcedNR_ELBM>::calculateAlpha", 5)

    T alphaMax = Base::calculateAlphaMax(haloDistributionNextPtr,
                                         haloDistributionPreviousPtr, iP);

    if (alphaMax < 2.) {
      Base::alpha = 0.95 * alphaMax;
    }

    else {
      T alphaMin = 1.;
      Base::alpha = Base::solveAlpha(haloDistributionNextPtr,
                                     haloDistributionPreviousPtr, iP, alphaMin,
                                     alphaMax);
    }
  }
};


template <class T, Architecture architecture>
  class Collision<T, CollisionType::ForcedNR_ELBM_Forcing, architecture>
  : public Collision<T, CollisionType::ForcedNR_ELBM, architecture> {
 private:
  using Base = Collision<T, CollisionType::ForcedNR_ELBM, architecture>;

 public:
  using Base::Collision;

  using Base::setForce;
  using Base::update;

  using Base::getAlpha;
  using Base::getDensity;
  using Base::getForce;
  using Base::getHydrodynamicVelocity;
  using Base::getVelocity;

  LBM_DEVICE LBM_HOST void
  calculateTs(const T* haloDistributionPreviousPtr, const T* haloDistributionNextPtr,
              const Position& iP) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateMoments", 4)

    if(writeT) {
      Moment_::calculateT_forcing(haloDistributionPreviousPtr, haloDistributionNextPtr, iP,
                                  Base::T2, Base::T3, Base::T4);
    }
  }

  LBM_DEVICE LBM_HOST inline
  void calculateRelaxationTime(T * haloDistributionNextPtr, T * haloDistributionPreviousPtr,
                               const Position& iP, const T alphaGuess) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ForcedNR_ELBM_Forcing>::calculate", 4)

    for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
      T equilibrium_iQ = Equilibrium_::calculate(Base::density, Base::velocity,
                                                 Base::velocity2, iQ);

      haloDistributionNextPtr[hSD::getIndex(iP, iQ)] =
        haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
        + Base::forcingScheme.calculateCollisionSource(Base::force, Base::density, Base::velocity,
                                                       Base::velocity2, equilibrium_iQ, iQ);
      haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)] -= equilibrium_iQ;
    }

    Base::alpha = alphaGuess;

    calculateAlpha(haloDistributionNextPtr, haloDistributionPreviousPtr, iP);

    Base::tau = (T)1.0 / (Base::alpha * Base::beta);
  }

  LBM_DEVICE LBM_HOST inline
    void collideAndStream(T* haloDistributionNext_Ptr, const T* haloDistributionPrevious_Ptr,
                          const Position& iP, const unsigned int iQ) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ForceNR_ELBM_Forcing>::calculate", 4)

    haloDistributionNext_Ptr[hSD::getIndex(iP, iQ)] -=
        (T)1.0 / Base::tau *
        haloDistributionPrevious_Ptr[hSD::getIndex(iP - uiL::celerity()[iQ],
                                                   iQ)];
  }


 private:
  LBM_DEVICE LBM_HOST inline
  void calculateAlpha(const T* haloDistributionNextPtr, const T* haloDistributionPreviousPtr,
                      const Position& iP) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ForcedNR_ELBM_Forcing>::calculateAlpha", 5)

    T alphaMax = calculateAlphaMax(haloDistributionNextPtr,
                                   haloDistributionPreviousPtr, iP);

    if (alphaMax < 2.) {
      Base::alpha = 0.95 * alphaMax;
    }

    else {
      T alphaMin = 1.;
      Base::alpha = solveAlpha(haloDistributionNextPtr, haloDistributionPreviousPtr,
                               iP, alphaMin, alphaMax);
    }
  }

  LBM_HOST LBM_DEVICE T calculateAlphaMax(const T* haloDistributionNext_Ptr,
                                          const T* haloDistributionPrevious_Ptr,
                                          const Position& iP) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ForcedNR_ELBM_Forcing>::calculateAlphaMax", 6)

    T alphaMaxR = 2.5;
    T alphaMaxTemp;

    for (unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
      if (haloDistributionPrevious_Ptr[hSD::getIndex(iP - uiL::celerity()[iQ],
                                                     iQ)] > 0) {
        alphaMaxTemp = fabs(haloDistributionNext_Ptr[hSD::getIndex(iP, iQ)] /
                            haloDistributionPrevious_Ptr[hSD::getIndex(
                                iP - uiL::celerity()[iQ], iQ)]);

        if (alphaMaxTemp < alphaMaxR) {
          alphaMaxR = alphaMaxTemp;
        }
      }
    }

    return alphaMaxR;
  }


  LBM_HOST LBM_DEVICE inline
  T solveAlpha(const T* haloDistributionNextPtr, const T* haloDistributionPreviousPtr,
               const Position& iP, const T alphaMin, const T alphaMax) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ForcedNR_ELBM_Forcing>::solveAlpha",6)

      EntropicStepFunctor<T, true> entropicStepFunctor(haloDistributionNextPtr,
                                                       haloDistributionPreviousPtr, iP);
    const T tolerance = 1e-8;
    const int iterationMax = 50;
    T alphaR = Base::alpha;

    bool hasConverged = NewtonRaphsonSolver(entropicStepFunctor, tolerance, iterationMax,
                                            alphaR, alphaMin, alphaMax);

    if(!hasConverged) {
      return 2.0;
    }

    return alphaR;
  }

};






template <class T, Architecture architecture>
  class Collision<T, CollisionType::ForcedBNR_ELBM, architecture>
  : public Collision<T, CollisionType::ForcedNR_ELBM, architecture> {
 private:
  using Base = Collision<T, CollisionType::ForcedNR_ELBM, architecture>;

 public:
  using Base::Collision;

  using Base::setForce;
  using Base::update;

  using Base::collideAndStream;

  using Base::getAlpha;
  using Base::getDensity;
  using Base::getForce;
  using Base::getHydrodynamicVelocity;
  using Base::getVelocity;

 private:
  LBM_DEVICE LBM_HOST inline T solveAlpha(const T* haloDistributionNextPtr,
                                          const T* haloDistributionPreviousPtr,
                                          const Position& iP,
                                          const T alphaMin,
                                          const T alphaMax) {
    LBM_INSTRUMENT_OFF("Collision<T, CollisionType::ForcedBNR_ELBM>::solveAlpha", 6)

    EntropicStepFunctor<T>
        entropicStepFunctor(haloDistributionNextPtr,
                            haloDistributionPreviousPtr, iP);
    const T tolerance = 1e-5;
    const int iterationMax = 20;
    T alphaR = Base::alpha;

    bool hasConverged =
        Bisection_NewtonRaphsonSolver(entropicStepFunctor, tolerance,
                                      iterationMax, alphaR, alphaMin, alphaMax);

    if (!hasConverged) {
      return 2.0;
    }

    return alphaR;
  }
};



template<Architecture architecture>
  using Collision_ = Collision<dataT, collisionT, architecture>;

}  // namespace lbm
