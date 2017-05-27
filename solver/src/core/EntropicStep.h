#ifndef ENTROPICSTEP_H
#define ENTROPICSTEP_H

#include <cmath>
#include <omp.h>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/sources/severity_feature.hpp>
namespace logging = boost::log;

#include "MathVector.h"
#include "Helpers.h"

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

}

#endif // ENTROPICSTEP_H
