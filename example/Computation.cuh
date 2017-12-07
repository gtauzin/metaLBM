#ifndef COMPUTATION_CUH
#define COMPUTATION_CUH

#include "Computation.h"

namespace lbm {

  template<typename Callback>
  __global__
  void kernel_1D(const MathVector<unsigned int, 3> start,
                 const MathVector<unsigned int, 3> end,
                 Callback function) {

    MathVector<unsigned int, 3> iP = {blockIdx.x*blockDim.x + threadIdx.x,
                                      start[1], start[2]};

    function(iP);
  }

  template<>
  struct Computation<Architecture::GPU, 1> {
    template<typename Callback>
    void Do(const MathVector<unsigned int, 3> start,
            const MathVector<unsigned int, 3> end,
            const MathVector<unsigned int, 3> length,
            Callback function){

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid(1, 1, 1);
      kernel_1D<Callback><<<dimBlock, dimGrid >>>(start, end, function);
    }

  };

}

#endif // COMPUTATION_CUH
