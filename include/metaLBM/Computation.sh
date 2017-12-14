#ifndef COMPUTATION_SH
#define COMPUTATION_SH

#include "Computation.cuh"
#include "Options.h"

namespace lbm {

  template<>
  class Computation<Architecture::SHMEM, 1>
  : public Computation<Architecture::GPU, 1> {
  public:
    using Computation<Architecture::GPU, 1>::Computation;
    using Computation<Architecture::GPU, 1>::Do;
  };

  template<>
  class Computation<Architecture::SHMEM, 2>
  : public Computation<Architecture::GPU, 2> {
  public:
    using Computation<Architecture::GPU, 2>::Computation;
    using Computation<Architecture::GPU, 2>::Do;
  };

  template<>
  class Computation<Architecture::SHMEM, 3>
  : public Computation<Architecture::GPU, 3> {
  public:
    using Computation<Architecture::GPU, 3>::Computation;
    using Computation<Architecture::GPU, 3>::Do;
  };

}

#endif // COMPUTATION_SH
