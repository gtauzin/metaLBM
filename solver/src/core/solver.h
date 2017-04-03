#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <random>
#include <array>
#include <memory>
#include <math.h>
#include <iomanip>
#include <omp.h>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/sources/severity_feature.hpp>
namespace logging = boost::log;

#include "input.h"
#include "structure.h"
#include "commons.h"
#include "helpers.h"

namespace lbm {

#pragma omp declare simd
  template <class T>
    struct EntropicStepFunctor : public RootFinderFunctor<T> {
  private:
    const MathVector<T, L::dimQ> f;
    const MathVector<T, L::dimQ> fNeq;

  public:
  EntropicStepFunctor(const MathVector<T, L::dimQ>& f_in,
                      const MathVector<T, L::dimQ>& fNeq_in)
    : RootFinderFunctor<T>()
      , f(f_in)
      , fNeq(fNeq_in)
    {}

#pragma omp declare simd
    inline T evaluateFunction(T const& alpha) {
      BOOST_LOG_TRIVIAL(debug) << "Evaluation of entropicStepFunction for alpha: "
                               << alpha;

      T entropicStepFunction = 0.0;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          BOOST_LOG_TRIVIAL(debug) << "f[" << iQ << "] = "
                                   << std::setprecision (17) << f[iQ] << ", "
                                   << "fNeq[" << iQ << "] = "
                                   << std::setprecision (17) << fNeq[iQ];

          T fmAlphafNeq_iQ = f[iQ] - alpha*fNeq[iQ];
          entropicStepFunction += f[iQ]*log(f[iQ]/L::weight()[iQ])
            - fmAlphafNeq_iQ*log(fmAlphafNeq_iQ/L::weight()[iQ]);
        });

      BOOST_LOG_TRIVIAL(debug) << "entropicStepFunction: "
                               << entropicStepFunction;

      return entropicStepFunction;
    }

#pragma omp declare simd
    inline T evaluateDerivative(T const& alpha) {
      BOOST_LOG_TRIVIAL(debug) << "Evaluation of entropicStepFunctionDerivative for alpha: "
                               << alpha;

      T entropicStepFunctionDerivative = 0.0;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          BOOST_LOG_TRIVIAL(debug) << "f[" << iQ << "] = "
                                   << std::setprecision (17) << f[iQ] << ", "
                                   << "fNeq[" << iQ << "] = "
                                   << std::setprecision (17) << fNeq[iQ];

          T fmAlphafNeq_iQ = f[iQ] - alpha*fNeq[iQ];

          entropicStepFunctionDerivative += fNeq[iQ]*(1 + log(fmAlphafNeq_iQ/L::weight()[iQ]));
        });

      BOOST_LOG_TRIVIAL(debug) << "entropicStepFunctionDerivative: "
                               << entropicStepFunctionDerivative;

      return entropicStepFunctionDerivative;
    }
  };


  template <class T>
    class Solver {
  public:

    Solver()
      {}

#pragma omp declare simd
    virtual T computeAlpha(const MathVector<T, L::dimQ>& f,
                           const MathVector<T, L::dimQ>& fNeq,
                           const T alphaGuess) = 0;
  };


  template <class T>
    class BGKSolver : public Solver<T> {
  public:

  BGKSolver()
    : Solver<T>()
      {}

#pragma omp declare simd
    inline T computeAlpha(const MathVector<T, L::dimQ>& f,
                          const MathVector<T, L::dimQ>& fNeq,
                          const T alphaGuess) {
      return 2.0;
    }

  };


  template <class T>
    class ELBMSolver : public Solver<T> {
  public:

  ELBMSolver()
    : Solver<T>()
      {}

#pragma omp declare simd
    inline T computeAlpha(const MathVector<T, L::dimQ>& f,
                          const MathVector<T, L::dimQ>& fNeq,
                          const T alphaGuess) {
      if(isDeviationSmall(f, fNeq)) {
        BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 2.0;
        return 2.0;
      }

      else {
        T alphaMax = calculateAlphaMax(f, fNeq);

        if(alphaMax < 2.) {
          BOOST_LOG_TRIVIAL(debug) << "alphaMax < 2: " << alphaMax;
          BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 0.95 * alphaMax;
          return 0.95 * alphaMax;
        }

        else {
          BOOST_LOG_TRIVIAL(debug) << "Starting NR as alphaMax > 2: " << alphaMax;
          T alphaMin = 1.;
          return calculateAlpha(f, fNeq, alphaGuess, alphaMin, alphaMax);

        }
      }
    }

#pragma omp declare simd
    inline bool isDeviationSmall(const MathVector<T, L::dimQ>& f,
                                 const MathVector<T, L::dimQ>& fNeq) {

      bool isDeviationSmallR = true;
      T deviation;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          deviation = fabs(fNeq[iQ]/f[iQ]);
          BOOST_LOG_TRIVIAL(debug) << "Calculation of deviation " << iQ << " : "
                                   << deviation;

          if(deviation > 1.0e-3) {
            BOOST_LOG_TRIVIAL(debug) << "Deviation > 1e-3: " << deviation;
            isDeviationSmallR = false;
          }
        });

      return isDeviationSmallR;
    }

#pragma omp declare simd
    T calculateAlphaMax(const MathVector<T, L::dimQ>& f,
                        const MathVector<T, L::dimQ>& fNeq) {

      T alphaMaxR = 2.5;
      T alphaMaxTemp;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          if(fNeq[iQ] > 0) {
            alphaMaxTemp = fabs(f[iQ]/fNeq[iQ]);

            if(alphaMaxTemp < alphaMaxR) {
              alphaMaxR = alphaMaxTemp;
            }
          }
        });

      BOOST_LOG_TRIVIAL(debug) << "Calculation of alphaMax : " << alphaMaxR;
      return alphaMaxR;
    }

#pragma omp declare simd
    inline T calculateAlpha(const MathVector<T, L::dimQ>& f,
                            const MathVector<T, L::dimQ>& fNeq,
                            const T alphaGuess,
                            const T alphaMin, const T alphaMax) {

      BOOST_LOG_TRIVIAL(debug) << "Starting with alpha: " << alphaGuess
                               << ", alphaMin: " << alphaMin
                               << ", alphaMax: " << alphaMax;


      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T>(f, fNeq));
      const T tolerance = 1e-5;
      const int iterationMax = 20;
      T alphaR = alphaGuess;

      bool hasConverged = NewtonRaphsonSolver(entropicStepFunctor,
                                              tolerance, iterationMax,
                                              alphaR, alphaMin, alphaMax);

      if(!hasConverged) {
        BOOST_LOG_TRIVIAL(debug) << "Unable to locate solution";
        BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 2.0;
        return 2.0;
      }

      BOOST_LOG_TRIVIAL(debug) << "Solution located: " << alphaR;
      BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << alphaR;
      return alphaR;
    }

  };

  template <class T>
    class Approached_ELBMSolver : public Solver<T> {
  public:

  Approached_ELBMSolver()
    : Solver<T>()
      {}

#pragma omp declare simd
    inline T computeAlpha(const MathVector<T, L::dimQ>& f,
                          const MathVector<T, L::dimQ>& fNeq,
                          const T alphaGuess) {
      if(isRelativeDeviationSmall(f, fNeq)) {
        BOOST_LOG_TRIVIAL(debug) << "Approximating alpha: " << 2.0;
        return approximateAlpha(f, fNeq);
      }
      else {
        T alphaMax = calculateAlphaMax(f, fNeq);

        if(alphaMax < 2.) {
          BOOST_LOG_TRIVIAL(debug) << "alphaMax < 2: " << alphaMax;
          BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 0.95 * alphaMax;
          return 0.95 * alphaMax;
        }

        else {
          BOOST_LOG_TRIVIAL(debug) << "Starting NR as alphaMax > 2: " << alphaMax;
          T alphaMin = 1.;
          return calculateAlpha(f, fNeq, alphaGuess, alphaMin, alphaMax);

        }
      }
    }

#pragma omp declare simd
    inline bool isRelativeDeviationSmall(const MathVector<T, L::dimQ>& f,
                                         const MathVector<T, L::dimQ>& fNeq) {

      bool isDeviationSmallR = true;
      T deviation;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          deviation = fabs(fNeq[iQ]/f[iQ]);
          BOOST_LOG_TRIVIAL(debug) << "Calculation of deviation " << iQ << " : "
                                   << deviation;

          if(deviation > 1.0e-3) {
            BOOST_LOG_TRIVIAL(debug) << "Deviation > 1e-3: " << deviation;
            isDeviationSmallR = false;
          }
        });

      return isDeviationSmallR;
    }


#pragma omp declare simd
    T approximateAlpha(const MathVector<T, L::dimQ>& f,
                       const MathVector<T, L::dimQ>& fNeq) {

      T a1 = 0.0;
      T a2 = 0.0;
      T a3 = 0.0;
      T a4 = 0.0;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          T temp = fNeq[iQ]/f[iQ];
          a1 += fNeq[iQ]*temp;
          a2 += fNeq[iQ]*temp*temp;
          a3 += fNeq[iQ]*temp*temp*temp;
          a4 += fNeq[iQ]*temp*temp*temp*temp;
        });

      a1 *= 1.0/2.0;
      a2 *= 1.0/6.0;
      a3 *= 1.0/12.0;
      a4 *= 1.0/20.0;

      T alphaR = 2 - 1/a1 * (4.0*a2 + 16.0*a2*a2/a1 - 8.0*a3
                             + 80.0*a2*a3/a1 - 80.0*a2*a2*a2/(a1*a1) - 16.0*a4);

      BOOST_LOG_TRIVIAL(debug) << "Approximation of alpha : " << alphaR;
      return alphaR;
    }

#pragma omp declare simd
    T calculateAlphaMax(const MathVector<T, L::dimQ>& f,
                        const MathVector<T, L::dimQ>& fNeq) {

      T alphaMaxR = 2.5;
      T alphaMaxTemp;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          if(fNeq[iQ] > 0) {
            alphaMaxTemp = fabs(f[iQ]/fNeq[iQ]);

            if(alphaMaxTemp < alphaMaxR) {
              alphaMaxR = alphaMaxTemp;
            }
          }
        });

      BOOST_LOG_TRIVIAL(debug) << "Calculation of alphaMax : " << alphaMaxR;
      return alphaMaxR;
    }

#pragma omp declare simd
    inline T calculateAlpha(const MathVector<T, L::dimQ>& f,
                            const MathVector<T, L::dimQ>& fNeq,
                            const T alphaGuess,
                            const T alphaMin, const T alphaMax) {

      BOOST_LOG_TRIVIAL(debug) << "Starting with alpha: " << alphaGuess
                               << ", alphaMin: " << alphaMin
                               << ", alphaMax: " << alphaMax;


      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T>(f, fNeq));
      const T tolerance = 1e-5;
      const int iterationMax = 20;
      T alphaR = alphaGuess;

      bool hasConverged = NewtonRaphsonSolver(entropicStepFunctor,
                                              tolerance, iterationMax,
                                              alphaR, alphaMin, alphaMax);

      if(!hasConverged) {
        BOOST_LOG_TRIVIAL(debug) << "Unable to locate solution";
        BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 2.0;
        return 2.0;
      }

      BOOST_LOG_TRIVIAL(debug) << "Solution located: " << alphaR;
      BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << alphaR;
      return alphaR;
    }

  };


  template<class T>
    class ForcedBNR_ELBMSolver : public Solver<T> {
  public:

  ForcedBNR_ELBMSolver()
    : Solver<T>()
      {}

#pragma omp declare simd
    inline T computeAlpha(const MathVector<T, L::dimQ>& f,
                          const MathVector<T, L::dimQ>& fNeq,
                          const T alphaGuess) {
      if(isDeviationSmall(f, fNeq)) {
        BOOST_LOG_TRIVIAL(debug) << "Deviation is small, should return alpha: : " << 2.0;
      }

      T alphaMax = calculateAlphaMax(f, fNeq);

      if(alphaMax < 2.) {
        BOOST_LOG_TRIVIAL(debug) << "alphaMax < 2: " << alphaMax;
        BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 0.95 * alphaMax;
        return 0.95 * alphaMax;
      }

      else {
        BOOST_LOG_TRIVIAL(debug) << "Starting NR as alphaMax > 2: " << alphaMax;
        T alphaMin = 1.;
        return calculateAlpha(f, fNeq, alphaGuess, alphaMin, alphaMax);

      }
    }

#pragma omp declare simd
    inline bool isDeviationSmall(const MathVector<T, L::dimQ>& f,
                                 const MathVector<T, L::dimQ>& fNeq) {

      bool isDeviationSmallR = true;
      T deviation;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          deviation = fabs(fNeq[iQ]/f[iQ]);
          BOOST_LOG_TRIVIAL(debug) << "Calculation of deviation " << iQ << " : "
                                   << deviation;

          if(deviation > 1.0e-3) {
            BOOST_LOG_TRIVIAL(debug) << "Deviation > 1e-3: " << deviation;
            isDeviationSmallR = false;
          }
        });

      return isDeviationSmallR;
    }

#pragma omp declare simd
    T calculateAlphaMax(const MathVector<T, L::dimQ>& f,
                        const MathVector<T, L::dimQ>& fNeq) {

      T alphaMaxR = 2.5;
      T alphaMaxTemp;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          if(fNeq[iQ] > 0) {
            alphaMaxTemp = fabs(f[iQ]/fNeq[iQ]);

            if(alphaMaxTemp < alphaMaxR) {
              alphaMaxR = alphaMaxTemp;
            }
          }
        });

      BOOST_LOG_TRIVIAL(debug) << "Calculation of alphaMax : " << alphaMaxR;
      return alphaMaxR;
    }

    inline T calculateAlpha(const MathVector<T, L::dimQ>& f,
                            const MathVector<T, L::dimQ>& fNeq,
                            const T alphaGuess,
                            const T alphaMin, const T alphaMax) {

      BOOST_LOG_TRIVIAL(debug) << "Starting with alpha: " << alphaGuess
                               << ", alphaMin: " << alphaMin
                               << ", alphaMax: " << alphaMax;


      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T>(f, fNeq));
      const T tolerance = 1e-5;
      const int iterationMax = 20;
      T alphaR = alphaGuess;

      bool hasConverged = Bisection_NewtonRaphsonSolver(entropicStepFunctor,
                                                        tolerance, iterationMax,
                                                        alphaR, alphaMin, alphaMax);

      if(!hasConverged) {
        BOOST_LOG_TRIVIAL(debug) << "Unable to locate solution";
        BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 2.0;
        return 2.0;
      }

      BOOST_LOG_TRIVIAL(debug) << "Solution located: " << alphaR;
      BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << alphaR;
      return alphaR;
    }
  };

  template<class T>
    class ForcedNR_ELBMSolver : public Solver<T> {
  public:

  ForcedNR_ELBMSolver()
    : Solver<T>()
      {}

#pragma omp declare simd
    inline T computeAlpha(const MathVector<T, L::dimQ>& f,
                          const MathVector<T, L::dimQ>& fNeq,
                          const T alphaGuess) {
      if(isDeviationSmall(f, fNeq)) {
        BOOST_LOG_TRIVIAL(debug) << "Deviation is small, should return alpha: : " << 2.0;
      }

      T alphaMax = calculateAlphaMax(f, fNeq);

      if(alphaMax < 2.) {
        BOOST_LOG_TRIVIAL(debug) << "alphaMax < 2: " << alphaMax;
        BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 0.95 * alphaMax;
        return 0.95 * alphaMax;
      }

      else {
        BOOST_LOG_TRIVIAL(debug) << "Starting NR as alphaMax > 2: " << alphaMax;
        T alphaMin = 1.;
        return calculateAlpha(f, fNeq, alphaGuess, alphaMin, alphaMax);

      }
    }

#pragma omp declare simd
    inline bool isDeviationSmall(const MathVector<T, L::dimQ>& f,
                                 const MathVector<T, L::dimQ>& fNeq) {

      bool isDeviationSmallR = true;
      T deviation;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          deviation = fabs(fNeq[iQ]/f[iQ]);
          BOOST_LOG_TRIVIAL(debug) << "Calculation of deviation " << iQ << " : "
                                   << deviation;

          if(deviation > 1.0e-3) {
            BOOST_LOG_TRIVIAL(debug) << "Deviation > 1e-3: " << deviation;
            isDeviationSmallR = false;
          }
        });

      return isDeviationSmallR;
    }

#pragma omp declare simd
    T calculateAlphaMax(const MathVector<T, L::dimQ>& f,
                        const MathVector<T, L::dimQ>& fNeq) {

      T alphaMaxR = 2.5;
      T alphaMaxTemp;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          if(fNeq[iQ] > 0) {
            alphaMaxTemp = fabs(f[iQ]/fNeq[iQ]);

            if(alphaMaxTemp < alphaMaxR) {
              alphaMaxR = alphaMaxTemp;
            }
          }
        });

      BOOST_LOG_TRIVIAL(debug) << "Calculation of alphaMax : " << alphaMaxR;
      return alphaMaxR;
    }

    inline T calculateAlpha(const MathVector<T, L::dimQ>& f,
                            const MathVector<T, L::dimQ>& fNeq,
                            const T alphaGuess,
                            const T alphaMin, const T alphaMax) {

      BOOST_LOG_TRIVIAL(debug) << "Starting with alpha: " << alphaGuess
                               << ", alphaMin: " << alphaMin
                               << ", alphaMax: " << alphaMax;


      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T>(f, fNeq));
      const T tolerance = 1e-5;
      const int iterationMax = 20;
      T alphaR = alphaGuess;

      bool hasConverged = NewtonRaphsonSolver(entropicStepFunctor, tolerance, iterationMax,
                                              alphaR, alphaMin, alphaMax);

      if(!hasConverged) {
        BOOST_LOG_TRIVIAL(debug) << "Unable to locate solution";
        BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 2.0;
        return 2.0;
      }

      BOOST_LOG_TRIVIAL(debug) << "Solution located: " << alphaR;
      BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << alphaR;
      return alphaR;
    }
  };


  template<class T>
    std::shared_ptr<Solver<T>> Create(const SolverMethod& solverMethod) {
    switch(solverMethod){
    case SolverMethod::BGK : {
      return std::shared_ptr<Solver<T>>(new BGKSolver<T>());
    }
    case SolverMethod::ELBM : {
      return std::shared_ptr<Solver<T>>(new ELBMSolver<T>());
    }
    case SolverMethod::Approached_ELBM : {
      return std::shared_ptr<Solver<T>>(new Approached_ELBMSolver<T>());
    }
    case SolverMethod::ForcedNR_ELBM : {
      return std::shared_ptr<Solver<T>>(new ForcedNR_ELBMSolver<T>());
    }
    case SolverMethod::ForcedBNR_ELBM : {
      return std::shared_ptr<Solver<T>>(new ForcedBNR_ELBMSolver<T>());
    }

    default:{
      BOOST_LOG_TRIVIAL(error) << "Unknown type of solver.";
      return nullptr;
    }
    }
  }

}

#endif // SOLVER_H
