#ifndef FUNCTOR_H
#define FUNCTOR_H

#include "metaLBM/MathVector.h"

namespace lbm {
  template<Architecture architecture>
    struct Functor {
      __host__ __device__
        void operator()(const MathVector<unsigned int, 3>& iP) {
      }

      void compute() {
        Computation<architecture, 1>().Do(MathVector<unsigned int, 3>{{0}},
                                         MathVector<unsigned int, 3>{{1}},
                                         MathVector<unsigned int, 3>{{1}},
                                         *this);
      }

  };
}

#endif // FUNCTOR_H
