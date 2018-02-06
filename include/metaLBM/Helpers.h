#ifndef HELPERS_H
#define HELPERS_H

#include <memory>
#include <cmath>

#include "Commons.h"

namespace lbm {

  template<int Begin, int End, int Step = 1>
  struct UnrolledFor {
    template<typename F>
    DEVICE HOST
    static void Do (F f) {
      f(Begin);
      UnrolledFor<Begin+Step, End, Step>::Do(f);
    }
  };

  template<int End>
  struct UnrolledFor<End, End> {
    template<typename F>
    DEVICE HOST
    static void Do (F f) {
    }
  };

  template <class T>
  DEVICE HOST
  inline T PowerBase(T arg, int power) {
    if(power == 1) {
      return arg;
    }
    else if(power == 0) {
      return (T) 1;
      }
    else if (power == -1) {
      return (T) 1.0/arg;
    }
    return (T) pow(arg, power);
  }


  template <class T, int power>
  class Power {
  public:
    DEVICE HOST
    static inline T Do(const T arg) {
      return (T) pow(arg, power);
    }
  };

  template <class T>
  class Power<T, 0> {
  public:
    DEVICE HOST
    static inline T Do(const T arg) {
      return (T) 1;
    }
  };

  template <class T>
  class Power<T, 1> {
  public:
    DEVICE HOST
    static inline T Do(const T arg) {
      return (T) arg;
    }
  };

  template <class T>
  class Power<T, 2> {
  public:
    DEVICE HOST
    static inline T Do(const T arg) {
      return (T) arg*arg;
    }
  };

  template <class T>
  class Power<T, 3> {
  public:
    DEVICE HOST
    static inline T Do(const T arg) {
      return (T) arg*arg*arg;
    }
  };

  template <class T>
  class Power<T, 4> {
  public:
    DEVICE HOST
    static inline T Do(const T arg) {
      return (T) arg*arg*arg*arg;
    }
  };

  template <class T>
  class Power<T, -1> {
  public:
    DEVICE HOST
    static inline T Do(const T arg) {
      return (T) 1.0/arg;
    }
  };


  template <class T>
  struct RootFinderFunctor {
  public:
    DEVICE HOST
    RootFinderFunctor(){};

    virtual T evaluateFunction(T const& x) = 0;

    virtual T evaluateDerivative(T const& x) = 0;

  };

  template <class T>
  inline bool NewtonRaphsonSolver(std::shared_ptr<RootFinderFunctor<T>> functor,
                                  const T tolerance, const int iterationMax,
                                  T& xR, const T xMin, const T xMax) {
    T error = 1 + tolerance;
    T xStep = 0.0;

    for(int iteration = 1; iteration <= iterationMax; ++iteration) {
      xR = xR - xStep;

      T functionEvaluation = functor->evaluateFunction(xR);
      T derivativeEvaluation = functor->evaluateDerivative(xR);
      xStep = functionEvaluation/derivativeEvaluation;

      error = fabs(xStep);

      if(error <= tolerance) {
        if(xR > xMin && xR < xMax) {
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
  inline bool Bisection_NewtonRaphsonSolver(std::shared_ptr<RootFinderFunctor<T>> functor,
                                            const T tolerance, const int iterationMax,
                                            T& xR, const T xMin, const T xMax) {

    T xLow, xHigh;
    T function_xLow = functor->evaluateFunction(xMin);
    T function_xHigh = functor->evaluateFunction(xMax);

    if ((function_xLow > 0.0 && function_xHigh > 0.0)
        || (function_xLow < 0.0 && function_xHigh < 0.0)) {
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
    //T xR = 2.0;
    T xStepPrevious = fabs(xMax-xMin);
    T xStep = xStepPrevious;
    T functionEvaluation = functor->evaluateFunction(xR);
    T derivativeEvaluation = functor->evaluateDerivative(xR);

    for(int iteration = 1; iteration <= iterationMax; ++iteration) {

      if ((((xR-xHigh)*derivativeEvaluation-functionEvaluation)
           *((xR-xLow)*derivativeEvaluation-functionEvaluation) > 0.0)
          || (fabs(2.0*functionEvaluation) > fabs(xStepPrevious*derivativeEvaluation))) {

        xStepPrevious = xStep;
        xStep = 0.5 * (xHigh-xLow);
        xR = xLow + xStep;

        if (xLow == xR) {
          return true;
        }
      }

      else {
        xStepPrevious = xStep;
        xStep = functionEvaluation/derivativeEvaluation;
        T xTemp = xR;
        xR -= xStep;

        if (xTemp == xR) {
          return true;
        }
      }

      if(fabs(xStep) <= tolerance) {
        return true;
      }

      functionEvaluation = functor->evaluateFunction(xR);
      derivativeEvaluation = functor->evaluateDerivative(xR);

      if (functionEvaluation < 0.0) {
        xLow = xR;
      }
      else {
        xHigh = xR;
      }
    }

    return false;

  }

}

#endif // HELPERS_H
