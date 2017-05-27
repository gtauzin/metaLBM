#ifndef HELPERS_H
#define HELPERS_H

#include <cmath>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

namespace lbm {

  template<int Begin, int End, int Step = 1>
  struct UnrolledFor {
    template<typename F>
    static void Do (F f) {
      f(Begin);
      UnrolledFor<Begin+Step, End, Step>::Do(f);
    }
  };

  template<int End>
  struct UnrolledFor<End, End> {
    template<typename F>
    static void Do (F f) {
    }
  };


#pragma omp declare simd
  template <class T, int power>
  class Power {
  public:
    inline T operator()(const T arg) {
      return (T) pow(arg, power);
    }
  };

#pragma omp declare simd
  template <class T>
  class Power<T, 0> {
  public:
    inline T operator()(const T arg) {
      return (T) 1;
    }
  };

#pragma omp declare simd
  template <class T>
  class Power<T, 1> {
  public:
    inline T operator()(const T arg) {

      return (T) arg;
    }
  };

#pragma omp declare simd
  template <class T>
  class Power<T, 2> {
  public:
    inline T operator()(const T arg) {
      return (T) arg*arg;
    }
  };

#pragma omp declare simd
  template <class T>
  class Power<T, 3> {
  public:
    inline T operator()(const T arg) {
      return (T) arg*arg*arg;
    }
  };

#pragma omp declare simd
  template <class T>
  class Power<T, 4> {
  public:
    inline T operator()(const T arg) {
      return (T) arg*arg*arg*arg;
    }
  };

#pragma omp declare simd
  template <class T>
  class Power<T, -1> {
  public:
    inline T operator()(const T arg) {
      return (T) 1.0/arg;
    }
  };



  template <class T>
  struct RootFinderFunctor {
  public:
    RootFinderFunctor(){};

#pragma omp declare simd
    virtual T evaluateFunction(T const& x) = 0;

    virtual T evaluateDerivative(T const& x) = 0;

  };

#pragma omp declare simd
  template <class T>
  inline bool NewtonRaphsonSolver(std::shared_ptr<RootFinderFunctor<T>> functor,
                                  const T tolerance, const int iterationMax,
                                  T& xR, const T xMin, const T xMax) {

    BOOST_LOG_TRIVIAL(debug) << "xR: " << xR
                             << ", xMin: " << xMin
                             << ", xMax: " << xMax
                             << ", tolerance: " << tolerance
                             << ", iterationMax: " << iterationMax;


    T error = 1 + tolerance;
    T xStep = 0.0;

    for(int iteration = 1; iteration <= iterationMax; ++iteration) {
      xR = xR - xStep;

      T functionEvaluation = functor->evaluateFunction(xR);
      T derivativeEvaluation = functor->evaluateDerivative(xR);
      xStep = functionEvaluation/derivativeEvaluation;

      error = fabs(xStep);

      BOOST_LOG_TRIVIAL(debug) << "iteration: " << iteration
                               << ", error: " << error
                               << ", xR: " << xR;

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

#pragma omp declare simd
  template <class T>
  inline bool Bisection_NewtonRaphsonSolver(std::shared_ptr<RootFinderFunctor<T>> functor,
                                            const T tolerance, const int iterationMax,
                                            T& xR, const T xMin, const T xMax) {

    BOOST_LOG_TRIVIAL(debug) << "xR: " << xR
                             << ", xMin: " << xMin
                             << ", xMax: " << xMax
                             << ", tolerance: " << tolerance
                             << ", iterationMax: " << iterationMax;

    T xLow, xHigh;
    T function_xLow = functor->evaluateFunction(xMin);
    T function_xHigh = functor->evaluateFunction(xMax);

    if ((function_xLow > 0.0 && function_xHigh > 0.0)
        || (function_xLow < 0.0 && function_xHigh < 0.0)) {
      BOOST_LOG_TRIVIAL(debug) << "Root must be in [xMin, xMax]";
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

      BOOST_LOG_TRIVIAL(debug) << "iteration: " << iteration
                               << ", error: " << fabs(xStep)
                               << ", xR: " << xR;

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

    BOOST_LOG_TRIVIAL(debug) << "Maximum number of iterations exceeded";
    return false;

  }

}

#endif // HELPERS_H
