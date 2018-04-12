#pragma once

#include <stdio.h>

#include "Computation.h"
#include "Commons.h"
#include "Options.h"
#include "Stream.cuh"

namespace lbm {

  template<typename Callback, typename... Arguments>
  GLOBAL
  void kernel_1D(const Position start, const Position end, const Position dir,
                 Callback function, const Arguments... arguments) {
    Position iP = start;
    iP[dir[0]] += blockIdx.x*blockDim.x + threadIdx.x;

    // delete if and make sure that this we choose blocks so that we are always
    // in the right case
    if(iP[dir[0]] < end[dir[0]]) {
      function(iP, arguments...);
    }
  }

  template<typename Callback, typename... Arguments>
  GLOBAL
  void kernel_2D(const Position start, const Position end, const Position dir,
                 Callback function, const Arguments... arguments) {
    Position iP = start;
    iP[dir[0]] += blockIdx.y*blockDim.y + threadIdx.y;
    iP[dir[1]] += blockIdx.x*blockDim.x + threadIdx.x;

    if(iP[dir[1]] < end[dir[0]] && iP[dir[0]] < end[dir[1]]) {
      function(iP, arguments...);
    }
  }

  template<typename Callback, typename... Arguments>
  GLOBAL
  void kernel_3D(const Position start, const Position end, const Position dir,
                 Callback function, const Arguments... arguments) {
    Position iP = start;
    iP[dir[0]] += blockIdx.z*blockDim.z + threadIdx.z;
    iP[dir[1]] += blockIdx.y*blockDim.y + threadIdx.y;
    iP[dir[2]] += blockIdx.x*blockDim.x + threadIdx.x;

    if(iP[dir[2]] < end[dir[0]] && iP[dir[1]] < end[dir[1]] && iP[dir[0]] < end[dir[2]]) {
      function(iP, arguments...);
    }
  }


 template<>
  class Computation<Architecture::GPU, 1>
    : public Computation<Architecture::GPU, 0> {
  private:
   using Base = Computation<Architecture::GPU, 0>;

  public:
   using Base::Computation;

    template<typename Callback, typename... Arguments>
    void Do(const Stream<Architecture::GPU>& stream, 
	    Callback function, const Arguments... arguments) {
      { INSTRUMENT_OFF("Computation<Architecture::GPU, 1>::Do<Callback>",3) }

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid((127+Base::length[dir[0]])/128, 1, 1);

      kernel_1D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end, Base::dir,
							function, arguments...);
      CUDA_CALL ( cudaGetLastError(); );
      CUDA_CALL( cudaDeviceSynchronize() );
    }

  };

 template<>
  class Computation<Architecture::GPU, 2>
    : public Computation<Architecture::GPU, 0> {
  private:
   using Base = Computation<Architecture::GPU, 0>;

  public:
   using Base::Computation;

    template<typename Callback, typename... Arguments>
    void Do(const Stream<Architecture::GPU>& stream, 
	    Callback function, const Arguments... arguments) {
      { INSTRUMENT_OFF("Computation<Architecture::GPU, 2>::Do<Callback>",3) }

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid((127+Base::length[dir[1]])/128, Base::length[dir[0]], 1);
      kernel_2D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end, Base::dir,
							function, arguments...);
      CUDA_CALL ( cudaGetLastError() );
      CUDA_CALL( cudaDeviceSynchronize() );
    }

 };

 template<>
  class Computation<Architecture::GPU, 3>
    : public Computation<Architecture::GPU, 0> {
  private:
   using Base = Computation<Architecture::GPU, 0>;

  public:
   using Base::Computation;

   template<typename Callback, typename... Arguments>
    void Do(const Stream<Architecture::GPU>& stream, 
	    Callback function, const Arguments... arguments) {
     { INSTRUMENT_OFF("Computation<Architecture::GPU, 3>::Do<Callback>",3) }

     dim3 dimBlock(128, 1, 1);
     dim3 dimGrid((127+Base::length[dir[2]])/128, Base::length[dir[1]], Base::length[dir[0]]);

     kernel_3D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end, Base::dir,
						       function, arguments...);
     CUDA_CALL ( cudaGetLastError() );
     CUDA_CALL( cudaDeviceSynchronize() );
   }

 };

}
