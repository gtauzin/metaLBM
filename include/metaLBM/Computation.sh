#ifndef COMPUTATION_SH
#define COMPUTATION_SH

#include "Computation.cuh"
#include "Options.h"

namespace lbm {

  template<>
  struct Computation<Architecture::SHMEM, 1>
  : public Computation<Architecture::GPU, 1> {
    using Computation<Architecture::GPU, 1>::Do;
  };

  template<>
  struct Computation<Architecture::SHMEM, 2>
  : public Computation<Architecture::GPU, 2> {
    using Computation<Architecture::GPU, 2>::Do;
  };

  template<>
  struct Computation<Architecture::SHMEM, 3>
  : public Computation<Architecture::GPU, 3> {
    using Computation<Architecture::GPU, 3>::Do;
  };

}

#endif // COMPUTATION_SH
