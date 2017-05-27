#ifndef COLLISION_H
#define COLLISION_H

#include <omp.h>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/sources/severity_feature.hpp>
namespace logging = boost::log;

#include "Options.h"
#include "MathVector.h"
#include "Helpers.h"
#include "Equilibrium.h"
#include "ForcingScheme.h"
#include "EntropicStep.h"

namespace lbm {

  template <class T, CollisionType collisionType>
  class Collision {};

  template <class T>
  class Collision<T, CollisionType::GenericSRT> {
  protected:
    T tau;

  public:
    #pragma omp declare simd
    inline void collide(T* __restrict__ f) {
      UnrolledFor<0, L::dimQ>::Do([&] (unsigned int iQ) {
          f = f[Domain<DomainType::Local>:: iQ)]
        + P::forcingScheme.getCollisionForcing<iQ>()
        - (T) 1.0 / tau * (f[idxPop(idx_lattice_previous, iQ)]
                           - P::equilibrium.computeEquilibrium<iQ>())
        ;
      });
    }

  };

  // fForced[iQ] - nextAlpha*beta*fNeq[iQ];

  template <class T>
  class Collision<T, CollisionType::BGK>
    : public Collision<T, CollisionType::GenericSRT> {

  public:
    inline T getAlpha() {
      return 2.0;
    }

  };

  template <class T>
  class Collision<T, CollisionType::ELBM> {
  protected:
    MathVector<T, L::dimQ>& f_Forced;
    MathVector<T, L::dimQ>& f_NonEq;

    #pragma omp declare simd
    inline bool isDeviationSmall(const T error) {

      bool isDeviationSmallR = true;
      T deviation;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          deviation = fabs(f_NonEq[iQ]/f_Forced[iQ]);
          BOOST_LOG_TRIVIAL(debug) << "Calculation of deviation " << iQ << " : "
                                   << deviation;

          if(deviation > error) {
            BOOST_LOG_TRIVIAL(debug) << "Deviation > " << error << ": " << deviation;
            isDeviationSmallR = false;
          }
        });

      return isDeviationSmallR;
    }

    #pragma omp declare simd
    T calculateAlphaMax() {

      T alphaMaxR = 2.5;
      T alphaMaxTemp;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          if(f_NonEq[iQ] > 0) {
            alphaMaxTemp = fabs(f_Forced[iQ]/f_NonEq[iQ]);

            if(alphaMaxTemp < alphaMaxR) {
              alphaMaxR = alphaMaxTemp;
            }
          }
        });

      BOOST_LOG_TRIVIAL(debug) << "Calculation of alphaMax : " << alphaMaxR;
      return alphaMaxR;
    }

    #pragma omp declare simd
    inline T calculateAlpha(const T alphaGuess,
                            const T alphaMin, const T alphaMax) {

      BOOST_LOG_TRIVIAL(debug) << "Starting with alpha: " << alphaGuess
                               << ", alphaMin: " << alphaMin
                               << ", alphaMax: " << alphaMax;


      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T>(f_Forced, f_NonEq));
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

  public:
    #pragma omp declare simd
    inline T computeAlpha(const T alphaGuess) {
      if(isDeviationSmall( (T) 1.0e-3)) {
        BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 2.0;
        return 2.0;
      }

      else {
        T alphaMax = calculateAlphaMax();

        if(alphaMax < 2.) {
          BOOST_LOG_TRIVIAL(debug) << "alphaMax < 2: " << alphaMax;
          BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 0.95 * alphaMax;
          return 0.95 * alphaMax;
        }

        else {
          BOOST_LOG_TRIVIAL(debug) << "Starting NR as alphaMax > 2: " << alphaMax;
          T alphaMin = 1.;
          return calculateAlpha(alphaGuess, alphaMin, alphaMax);

        }
      }
    }

  };

  template <class T>
  class Collision<T, CollisionType::Approached_ELBM>
    : public Collision<T, CollisionType::ELBM> {
  private:
    using Collision<T, CollisionType::ELBM>::isRelativeDeviationSmall;
    using Collision<T, CollisionType::ELBM>::calculateAlphaMax;
    using Collision<T, CollisionType::ELBM>::calculateAlpha;

    #pragma omp declare simd
    inline T approximateAlpha() {

      T a1 = 0.0;
      T a2 = 0.0;
      T a3 = 0.0;
      T a4 = 0.0;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          T temp = f_NonEq[iQ]/f_Forced[iQ];
          a1 += f_NonEq[iQ]*temp;
          a2 += f_NonEq[iQ]*temp*temp;
          a3 += f_NonEq[iQ]*temp*temp*temp;
          a4 += f_NonEq[iQ]*temp*temp*temp*temp;
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

  public:
    #pragma omp declare simd
    inline T computeAlpha(const T alphaGuess) {
      if(isRelativeDeviationSmall( (T) 1.0e-3)) {
        T alphaApproximated = approximateAlpha();
        BOOST_LOG_TRIVIAL(debug) << "Approximating alpha: " << alphaApproximated;
        return alphaApproximated;
      }
      else {
        T alphaMax = calculateAlphaMax();

        if(alphaMax < 2.) {
          BOOST_LOG_TRIVIAL(debug) << "alphaMax < 2: " << alphaMax;
          BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 0.95 * alphaMax;
          return 0.95 * alphaMax;
        }

        else {
          BOOST_LOG_TRIVIAL(debug) << "Starting root solver as alphaMax > 2: " << alphaMax;
          T alphaMin = 1.;
          return calculateAlpha(alphaGuess, alphaMin, alphaMax);

        }
      }
    }

  };

  template <class T>
  class Collision<T, CollisionType::ForcedNR_ELBM>
  : public Collision<T, CollisionType::ELBM> {
  private:
    using Collision<T, CollisionType::ELBM>::isRelativeDeviationSmall;
    using Collision<T, CollisionType::ELBM>::calculateAlphaMax;

    inline T calculateAlpha(const T alphaGuess,
                            const T alphaMin, const T alphaMax) {

      BOOST_LOG_TRIVIAL(debug) << "Starting with alpha: " << alphaGuess
                               << ", alphaMin: " << alphaMin
                               << ", alphaMax: " << alphaMax;


      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T>(f_Forced, f_NonEq));
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

  public:
    #pragma omp declare simd
    inline T computeAlpha(const T alphaGuess) {
      if(isDeviationSmall( (T) 1.0e-3)) {
        BOOST_LOG_TRIVIAL(debug) << "Deviation is small, should return alpha: : " << 2.0;
      }

      T alphaMax = calculateAlphaMax();

      if(alphaMax < 2.) {
        BOOST_LOG_TRIVIAL(debug) << "alphaMax < 2: " << alphaMax;
        BOOST_LOG_TRIVIAL(debug) << "Returning alpha: " << 0.95 * alphaMax;
        return 0.95 * alphaMax;
      }

      else {
        BOOST_LOG_TRIVIAL(debug) << "Starting root solver as alphaMax > 2: " << alphaMax;
        T alphaMin = 1.;
        return calculateAlpha(alphaGuess, alphaMin, alphaMax);

      }
    }

  };

  template <class T>
  class Collision<T, CollisionType::ForcedBNR_ELBM>
    : public Collision<T, CollisionType::ForcedNR_ELBM> {
  private:
    using Collision<T, CollisionType::ForcedNR_ELBM>::isRelativeDeviationSmall;
    using Collision<T, CollisionType::ForcedNR_ELBM>::calculateAlphaMax;

    inline T calculateAlpha(const T alphaGuess,
                            const T alphaMin, const T alphaMax) {

      BOOST_LOG_TRIVIAL(debug) << "Starting with alpha: " << alphaGuess
                               << ", alphaMin: " << alphaMin
                               << ", alphaMax: " << alphaMax;


      std::shared_ptr<RootFinderFunctor<T>> entropicStepFunctor =
        std::shared_ptr<RootFinderFunctor<T>>(new EntropicStepFunctor<T>(f_Forced, f_NonEq));
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

  public:
    using Collision<T, CollisionType::ForcedNR_ELBM>::computeAlpha;

  };

}

#endif // COLLISION_H
