#pragma once

#include <stdio.h>

#include "Computation.h"
#include "Commons.h"
#include "Options.h"

namespace lbm {

  template<typename Callback, typename... Arguments>
  GLOBAL
  void kernel_1D(const Position start,
                 const Position end,
                 Callback function, const Arguments... arguments) {

    Position iP = {blockIdx.x*blockDim.x + threadIdx.x+ start[d::X],
                   start[d::Y], start[d::Z]};

    // delete if and make sure that this we choose blocks so that we are always
    // in the right case
    if(iP[0] < end[d::X]) {
      function(iP, arguments...);
    }
  }

  template<typename Callback, typename... Arguments>
  GLOBAL
  void kernel_2D(const Position start,
                 const Position end,
                 Callback function, const Arguments... arguments) {

    Position iP = {blockIdx.y*blockDim.y + threadIdx.y + start[d::X],
                   blockIdx.x*blockDim.x + threadIdx.x + start[d::Y],
                   start[d::Z]};

    if(iP[1] < end[d::X] && iP[0] < end[d::Y]) {
      function(iP, arguments...);
    }
  }

  template<typename Callback, typename... Arguments>
  GLOBAL
  void kernel_3D(const Position start,
                 const Position end,
                 Callback function, const Arguments... arguments) {

    Position iP = {blockIdx.z*blockDim.z + threadIdx.z + start[d::X],
                   blockIdx.y*blockDim.y + threadIdx.y + start[d::Y],
                   blockIdx.x*blockDim.x + threadIdx.x + start[d::Z]};
    if(iP[2] < end[d::X] && iP[1] < end[d::Y] && iP[0] < end[d::Z]) {
      function(iP, arguments...);
    }
  }


 template<>
  class Computation<Architecture::GPU, 1>
    : public Computation<Architecture::Generic, 1> {
  private:
   using Computation<Architecture::Generic, 1>::start;
   using Computation<Architecture::Generic, 1>::end;
   using Computation<Architecture::Generic, 1>::length;

  public:
   using Computation<Architecture::Generic, 1>::Computation;

    template<typename Callback, typename... Arguments>
    void Do(Callback function, const Arguments... arguments) {
      INSTRUMENT_OFF("Computation<Architecture::GPU, 1>::Do<Callback>",3)

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid((127+length[d::X])/128, 1, 1);
      kernel_1D<Callback><<<dimGrid, dimBlock>>>(start, end, function, arguments...);
      CUDA_CALL ( cudaGetLastError(); )
      CUDA_CALL( cudaDeviceSynchronize(); )
    }
  };

 template<>
  class Computation<Architecture::GPU, 2>
    : public Computation<Architecture::Generic, 2> {
  private:
   using Computation<Architecture::Generic, 2>::start;
   using Computation<Architecture::Generic, 2>::end;
   using Computation<Architecture::Generic, 2>::length;

  public:
   using Computation<Architecture::Generic, 2>::Computation;

    template<typename Callback, typename... Arguments>
    void Do(Callback function, const Arguments... arguments) {
      INSTRUMENT_OFF("Computation<Architecture::GPU, 2>::Do<Callback>",3)

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid((127+length[d::Y])/128, length[d::X], 1);
      kernel_2D<<<dimGrid, dimBlock>>>(start, end, function, arguments...);
      CUDA_CALL ( cudaGetLastError(); )
      CUDA_CALL( cudaDeviceSynchronize(); )
    }
  };

 template<>
  class Computation<Architecture::GPU, 3>
    : public Computation<Architecture::Generic, 3> {
  private:
   using Computation<Architecture::Generic, 3>::start;
   using Computation<Architecture::Generic, 3>::end;
   using Computation<Architecture::Generic, 3>::length;

  public:
   using Computation<Architecture::Generic, 3>::Computation;

    template<typename Callback, typename... Arguments>
    void Do(Callback function, const Arguments... arguments) {
      INSTRUMENT_OFF("Computation<Architecture::GPU, 3>::Do<Callback>",3)

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid((127+length[d::Z])/128, length[d::Y], length[d::X]);

      kernel_3D<Callback><<<dimGrid, dimBlock>>>(start, end, function, arguments...);
      CUDA_CALL ( cudaGetLastError(); )
      CUDA_CALL( cudaDeviceSynchronize(); )
    }
  };

}
