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
  template <class T, LatticeType L>
    T computeDensity(const T * __restrict__ f, const int idx) {
    T newDensity = f[idxPop<T, L>(idx, 0)];

    UnrolledFor<1, dimQ<T, L>()>::Do([&] (int iQ) {
        newDensity += f[idxPop<T, L>(idx, iQ)];
      });
    return newDensity;
  }


#pragma omp declare simd
  template <class T, LatticeType L>
    inline MathVector<T, dimD<T, L>()> computeVelocity(const T * __restrict__ f,
                                                        const int idx,
                                                        const T density) {
    MathVector<T, dimD<T, L>()> velocityR(celerity<T, L>(0)* f[idxPop<T, L>(idx, 0)]);

    UnrolledFor<1, dimQ<T, L>()>::Do([&] (int iQ) {
        velocityR += celerity<T, L>(iQ)* f[idxPop<T, L>(idx, iQ)];
      });

    return velocityR/density;
  }


#pragma omp declare simd
  template <class T, LatticeType L>
  inline T powerOfCelerity(const T arg, const int celerity) {
    if(celerity == 0) {
      return 1.0;
    }
    else if(celerity == -1) {
      return 1.0/arg;
    }
    else {
      return arg;
    }
  }

#pragma omp declare simd
  template <class T, LatticeType L>
    inline T computeEquilibrium(const int iQ, T rho,
                                const MathVector<T, dimD<T, L>()>& u,
                                const T u2) {

    T fEq_iQ = 1.0;

    UnrolledFor<0, dimD<T, L>()>::Do([&] (int d) {
        fEq_iQ *= (2.0 - sqrt(1.0 + 3.0*u[d]*u[d]))
          * powerOfCelerity<T, L>((2* u[d] + sqrt(1.0 + 3.0*u[d]*u[d]))/(1.0 - u[d]),
                                  (int) celerity<T, L>(d, iQ));
      });

    return rho * weight<T, L>(iQ)*fEq_iQ;
  }


#pragma omp declare simd
  template <class T, LatticeType L>
    struct EntropicStepFunctor : public RootFinderFunctor<T> {
  private:
    const MathVector<T, dimQ<T, L>()> f;
    const MathVector<T, dimQ<T, L>()> fNeq;

  public:
  EntropicStepFunctor(const MathVector<T, dimQ<T, L>()>& f_in,
                      const MathVector<T, dimQ<T, L>()>& fNeq_in)
    : RootFinderFunctor<T>()
      , f(f_in)
      , fNeq(fNeq_in)
    {}

#pragma omp declare simd
    inline T evaluateFunction(T const& alpha) {
      BOOST_LOG_TRIVIAL(debug) << "Evaluation of entropicStepFunction for alpha: "
                               << alpha;

      T entropicStepFunction = 0.0;

      UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
          BOOST_LOG_TRIVIAL(debug) << "f[" << iQ << "] = "
                                   << std::setprecision (17) << f[iQ] << ", "
                                   << "fNeq[" << iQ << "] = "
                                   << std::setprecision (17) << fNeq[iQ];

          T fmAlphafNeq_iQ = f[iQ] - alpha*fNeq[iQ];
          entropicStepFunction += f[iQ]*log(f[iQ]/weight<T, L>(iQ))
            - fmAlphafNeq_iQ*log(fmAlphafNeq_iQ/weight<T, L>(iQ));
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

      UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
          BOOST_LOG_TRIVIAL(debug) << "f[" << iQ << "] = "
                                   << std::setprecision (17) << f[iQ] << ", "
                                   << "fNeq[" << iQ << "] = "
                                   << std::setprecision (17) << fNeq[iQ];

          T fmAlphafNeq_iQ = f[iQ] - alpha*fNeq[iQ];

          entropicStepFunctionDerivative += fNeq[iQ]*(1 + log(fmAlphafNeq_iQ/weight<T, L>(iQ)));
        });

      BOOST_LOG_TRIVIAL(debug) << "entropicStepFunctionDerivative: "
                               << entropicStepFunctionDerivative;

      return entropicStepFunctionDerivative;
    }
  };


  template <class T, LatticeType L>
    class Solver {
  public:

    Solver()
      {}

#pragma omp declare simd
    virtual T computeAlpha(const MathVector<T, dimQ<T, L>()>& f,
                           const MathVector<T, dimQ<T, L>()>& fNeq,
                           const T alphaGuess) = 0;
  };


  template <class T, LatticeType L>
    class BGKSolver : public Solver<T, L> {
  public:

  BGKSolver()
    : Solver<T, L>()
      {}

#pragma omp declare simd
    inline T computeAlpha(const MathVector<T, dimQ<T, L>()>& f,
                          const MathVector<T, dimQ<T, L>()>& fNeq,
                          const T alphaGuess) {
      return 2.0;
    }

  };


  template <class T, LatticeType L>
    class ELBMSolver : public Solver<T, L> {
  public:

  ELBMSolver()
    : Solver<T, L>()
      {}

#pragma omp declare simd
    inline T computeAlpha(const MathVector<T, dimQ<T, L>()>& f,
                          const MathVector<T, dimQ<T, L>()>& fNeq,
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
    inline bool isDeviationSmall(const MathVector<T, dimQ<T, L>()>& f,
                                 const MathVector<T, dimQ<T, L>()>& fNeq) {

      T deviation;

      UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
          deviation = fabs(fNeq[iQ]/f[iQ]);
          BOOST_LOG_TRIVIAL(debug) << "Calculation of deviation " << iQ << " : " << deviation;

          if(deviation > 1.0e-3) {
            BOOST_LOG_TRIVIAL(debug) << "Deviation > 1e-3: " << deviation;
            return false;
          }
        });

      return true;
    }

#pragma omp declare simd
    T calculateAlphaMax(const MathVector<T, dimQ<T, L>()>& f,
                        const MathVector<T, dimQ<T, L>()>& fNeq) {

      T alphaMaxR = 2.5;
      T alphaMaxTemp;

      UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
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
    inline T calculateAlpha(const MathVector<T, dimQ<T, L>()>& f,
                            const MathVector<T, dimQ<T, L>()>& fNeq,
                            const T alphaGuess,
                            const T alphaMin, const T alphaMax) {

      BOOST_LOG_TRIVIAL(debug) << "Starting with alpha: " << alphaGuess
                               << ", alphaMin: " << alphaMin
                               << ", alphaMax: " << alphaMax;


      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T, L>(f, fNeq));
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

  template <class T, LatticeType L>
    class Approached_ELBMSolver : public Solver<T, L> {
  public:

  Approached_ELBMSolver()
    : Solver<T, L>()
      {}

#pragma omp declare simd
    inline T computeAlpha(const MathVector<T, dimQ<T, L>()>& f,
                          const MathVector<T, dimQ<T, L>()>& fNeq,
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
    inline bool isRelativeDeviationSmall(const MathVector<T, dimQ<T, L>()>& f,
                                         const MathVector<T, dimQ<T, L>()>& fNeq) {

      T deviation;

      UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
          deviation = fabs(fNeq[iQ]/f[iQ]);
          BOOST_LOG_TRIVIAL(debug) << "Calculation of relative deviation " << iQ << " : "
                                   << deviation;

          if(deviation > 1.0e-2) {
            BOOST_LOG_TRIVIAL(debug) << "Relative deviation > 1e-2: " << deviation;
            return false;
          }
        });

      return true;
    }


#pragma omp declare simd
    T approximateAlpha(const MathVector<T, dimQ<T, L>()>& f,
                       const MathVector<T, dimQ<T, L>()>& fNeq) {

      T a1 = 0.0;
      T a2 = 0.0;
      T a3 = 0.0;
      T a4 = 0.0;

      UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
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
    T calculateAlphaMax(const MathVector<T, dimQ<T, L>()>& f,
                        const MathVector<T, dimQ<T, L>()>& fNeq) {

      T alphaMaxR = 2.5;
      T alphaMaxTemp;

      UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
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
    inline T calculateAlpha(const MathVector<T, dimQ<T, L>()>& f,
                            const MathVector<T, dimQ<T, L>()>& fNeq,
                            const T alphaGuess,
                            const T alphaMin, const T alphaMax) {

      BOOST_LOG_TRIVIAL(debug) << "Starting with alpha: " << alphaGuess
                               << ", alphaMin: " << alphaMin
                               << ", alphaMax: " << alphaMax;


      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T, L>(f, fNeq));
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


  template<class T, LatticeType L>
    class ForcedBNR_ELBMSolver : public Solver<T, L> {
  public:

  ForcedBNR_ELBMSolver()
    : Solver<T, L>()
      {}

#pragma omp declare simd
    inline T computeAlpha(const MathVector<T, dimQ<T, L>()>& f,
                          const MathVector<T, dimQ<T, L>()>& fNeq,
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
    inline bool isDeviationSmall(const MathVector<T, dimQ<T, L>()>& f,
                                 const MathVector<T, dimQ<T, L>()>& fNeq) {

      T deviation;

      UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
          deviation = fabs(fNeq[iQ]/f[iQ]);
          BOOST_LOG_TRIVIAL(debug) << "Calculation of deviation " << iQ << " : "
                                   << deviation;

          if(deviation > 1.0e-3) {
            BOOST_LOG_TRIVIAL(debug) << "Deviation > 1e-3: " << deviation;
            return false;
          }
        });

      return true;
    }

#pragma omp declare simd
    T calculateAlphaMax(const MathVector<T, dimQ<T, L>()>& f,
                        const MathVector<T, dimQ<T, L>()>& fNeq) {

      T alphaMaxR = 2.5;
      T alphaMaxTemp;

      UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
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

    inline T calculateAlpha(const MathVector<T, dimQ<T, L>()>& f,
                            const MathVector<T, dimQ<T, L>()>& fNeq,
                            const T alphaGuess,
                            const T alphaMin, const T alphaMax) {

      BOOST_LOG_TRIVIAL(debug) << "Starting with alpha: " << alphaGuess
                               << ", alphaMin: " << alphaMin
                               << ", alphaMax: " << alphaMax;


      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T, L>(f, fNeq));
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

  template<class T, LatticeType L>
    class ForcedNR_ELBMSolver : public Solver<T, L> {
  public:

  ForcedNR_ELBMSolver()
    : Solver<T, L>()
      {}

#pragma omp declare simd
    inline T computeAlpha(const MathVector<T, dimQ<T, L>()>& f,
                          const MathVector<T, dimQ<T, L>()>& fNeq,
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
    inline bool isDeviationSmall(const MathVector<T, dimQ<T, L>()>& f,
                                 const MathVector<T, dimQ<T, L>()>& fNeq) {

      T deviation;

      UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
          deviation = fabs(fNeq[iQ]/f[iQ]);
          BOOST_LOG_TRIVIAL(debug) << "Calculation of deviation " << iQ << " : "
                                   << deviation;

          if(deviation > 1.0e-3) {
            BOOST_LOG_TRIVIAL(debug) << "Deviation > 1e-3: " << deviation;
            return false;
          }
        });

      return true;
    }

#pragma omp declare simd
    T calculateAlphaMax(const MathVector<T, dimQ<T, L>()>& f,
                        const MathVector<T, dimQ<T, L>()>& fNeq) {

      T alphaMaxR = 2.5;
      T alphaMaxTemp;

      UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
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

    inline T calculateAlpha(const MathVector<T, dimQ<T, L>()>& f,
                            const MathVector<T, dimQ<T, L>()>& fNeq,
                            const T alphaGuess,
                            const T alphaMin, const T alphaMax) {

      BOOST_LOG_TRIVIAL(debug) << "Starting with alpha: " << alphaGuess
                               << ", alphaMin: " << alphaMin
                               << ", alphaMax: " << alphaMax;


      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T, L>(f, fNeq));
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


  template<class T, LatticeType L>
    std::shared_ptr<Solver<T, L>> Create(const SolverMethod& solverMethod) {
    switch(solverMethod){
    case SolverMethod::BGK : {
      return std::shared_ptr<Solver<T, L>>(new BGKSolver<T, L>());
    }
    case SolverMethod::ELBM : {
      return std::shared_ptr<Solver<T, L>>(new ELBMSolver<T, L>());
    }
    case SolverMethod::Approached_ELBM : {
      return std::shared_ptr<Solver<T, L>>(new Approached_ELBMSolver<T, L>());
    }
    case SolverMethod::ForcedNR_ELBM : {
      return std::shared_ptr<Solver<T, L>>(new ForcedNR_ELBMSolver<T, L>());
    }
    case SolverMethod::ForcedBNR_ELBM : {
      return std::shared_ptr<Solver<T, L>>(new ForcedBNR_ELBMSolver<T, L>());
    }

    default:{
      BOOST_LOG_TRIVIAL(error) << "Unknown type of solver.";
      return nullptr;
    }
    }
  }

}

#endif // SOLVER_H
