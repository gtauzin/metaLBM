#ifndef ENTROPICSTEP_H
#define ENTROPICSTEP_H

#include <cmath>
#include <omp.h>

#include "MathVector.h"
#include "Helpers.h"

namespace lbm {

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
      T entropicStepFunction = 0.0;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          T fmAlphafNeq_iQ = f[iQ] - alpha*fNeq[iQ];
          entropicStepFunction += f[iQ]*log(f[iQ]/L::weight()[iQ])
            - fmAlphafNeq_iQ*log(fmAlphafNeq_iQ/L::weight()[iQ]);
        });

      return entropicStepFunction;
    }

#pragma omp declare simd
    inline T evaluateDerivative(T const& alpha) {
      T entropicStepFunctionDerivative = 0.0;

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          T fmAlphafNeq_iQ = f[iQ] - alpha*fNeq[iQ];

          entropicStepFunctionDerivative += fNeq[iQ]*(1 + log(fmAlphafNeq_iQ/L::weight()[iQ]));
        });

      return entropicStepFunctionDerivative;
    }
  };

}

#endif // ENTROPICSTEP_H
