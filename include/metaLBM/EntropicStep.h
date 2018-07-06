#pragma once

#include <cmath>

#include "Commons.h"
#include "Helpers.h"
#include "Lattice.h"
#include "MathVector.h"

namespace lbm {

  template <class T>
  struct EntropicStepFunctor {
  private:
    const T* haloDistributionNext_Ptr;
    const T* haloDistributionPrevious_Ptr;
    const Position iP;

  public:
    LBM_HOST LBM_DEVICE
    EntropicStepFunctor(const T* haloDistributionNext_Ptr_in,
                        const T* haloDistributionPrevious_Ptr_in,
                        const Position& iP_in)
      : haloDistributionNext_Ptr(haloDistributionNext_Ptr_in),
        haloDistributionPrevious_Ptr(haloDistributionPrevious_Ptr_in),
        iP(iP_in) {}

    LBM_HOST LBM_DEVICE inline
    T evaluateFunction(T const& alpha) {
      T entropicStepFunction = (T)0;

      for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
        T f_iQ = haloDistributionPrevious_Ptr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)];
        T fmAlphafNeq_iQ =
          f_iQ - alpha * haloDistributionNext_Ptr[hSD::getIndex(iP, iQ)];;
        entropicStepFunction +=
          f_iQ * log(f_iQ / L::weight()[iQ]) -
          fmAlphafNeq_iQ * log(fmAlphafNeq_iQ / L::weight()[iQ]);
      }

      return entropicStepFunction;
    }

    LBM_HOST LBM_DEVICE inline
    T evaluateDerivative(T const& alpha) {
      T entropicStepFunctionDerivative = (T)0;

      for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
        T fNeq_iQ = haloDistributionNext_Ptr[hSD::getIndex(iP, iQ)];
        T fmAlphafNeq_iQ =
          haloDistributionPrevious_Ptr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
          - alpha * fNeq_iQ;

        entropicStepFunctionDerivative +=
          fNeq_iQ * (1 + log(fmAlphafNeq_iQ / L::weight()[iQ]));
      }

      return entropicStepFunctionDerivative;
    }
  };

  template <class T>
  LBM_HOST LBM_DEVICE inline
  bool NewtonRaphsonSolver(EntropicStepFunctor<T> functor,
                           const T tolerance, const int iterationMax,
                           T& xR, const T xMin, const T xMax) {
    T error = 1 + tolerance;
    T xStep = 0.0;

    for (int iteration = 1; iteration <= iterationMax; ++iteration) {
      xR = xR - xStep;

      T functionEvaluation = functor.evaluateFunction(xR);
      T derivativeEvaluation = functor.evaluateDerivative(xR);
      xStep = functionEvaluation / derivativeEvaluation;

      error = fabs(xStep);

      if (error <= tolerance) {
        if (xR > xMin && xR < xMax) {
          return true;
        }

        else {
          return false;
        }
      }
    }

    return false;
  }

  template <class T>
  LBM_HOST LBM_DEVICE inline
  bool Bisection_NewtonRaphsonSolver(EntropicStepFunctor<T> functor,
                                     const T tolerance, const int iterationMax,
                                     T& xR, const T xMin, const T xMax) {
    T xLow, xHigh;
    T function_xLow = functor.evaluateFunction(xMin);
    T function_xHigh = functor.evaluateFunction(xMax);

    if ((function_xLow > 0.0 && function_xHigh > 0.0) ||
        (function_xLow < 0.0 && function_xHigh < 0.0)) {
      return false;
    }

    if (function_xLow == 0.0) {
      xR = xMin;
      return true;
    }

    if (function_xHigh == 0.0) {
      xR = xMax;
      return true;
    }

    if (function_xLow < 0.0) {
      xLow = xMin;
      xHigh = xMax;
    }

    else {
      xLow = xMax;
      xHigh = xMin;
    }

    xR = 0.5 * (xMin + xMax);
    // T xR = 2.0;
    T xStepPrevious = fabs(xMax - xMin);
    T xStep = xStepPrevious;
    T functionEvaluation = functor.evaluateFunction(xR);
    T derivativeEvaluation = functor.evaluateDerivative(xR);

    for (int iteration = 1; iteration <= iterationMax; ++iteration) {
      if ((((xR - xHigh) * derivativeEvaluation - functionEvaluation) *
           ((xR - xLow) * derivativeEvaluation - functionEvaluation) > 0.0)
          || (fabs(2.0 * functionEvaluation) > fabs(xStepPrevious * derivativeEvaluation))) {
        xStepPrevious = xStep;
        xStep = 0.5 * (xHigh - xLow);
        xR = xLow + xStep;

        if (xLow == xR) {
          return true;
        }
      } else {
        xStepPrevious = xStep;
        xStep = functionEvaluation / derivativeEvaluation;
        T xTemp = xR;
        xR -= xStep;

        if (xTemp == xR) {
          return true;
        }
      }

      if (fabs(xStep) <= tolerance) {
        return true;
      }

      functionEvaluation = functor.evaluateFunction(xR);
      derivativeEvaluation = functor.evaluateDerivative(xR);

      if (functionEvaluation < 0.0) {
        xLow = xR;
      } else {
        xHigh = xR;
      }
    }

    return false;
  }

}  // namespace lbm
