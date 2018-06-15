#pragma once

#include <stdio.h>

#include "Computation.h"
#include "Commons.h"
#include "Options.h"
#include "Stream.cuh"

namespace lbm {

  template<typename Callback, typename... Arguments>
  LBM_GLOBAL
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
  LBM_GLOBAL
  void kernel_2D(const Position start, const Position end, const Position dir,
                 Callback function, const Arguments... arguments) {
    Position iP = start;
    iP[dir[0]] += blockIdx.y*blockDim.y + threadIdx.y;
    iP[dir[1]] += blockIdx.x*blockDim.x + threadIdx.x;

    if(iP[dir[0]] < end[dir[0]] && iP[dir[1]] < end[dir[1]]) {
      function(iP, arguments...);
    }
  }

  template<typename Callback, typename... Arguments>
  LBM_GLOBAL
  void kernel_3D(const Position start, const Position end, const Position dir,
                 Callback function, const Arguments... arguments) {
    Position iP = start;
    iP[dir[0]] += blockIdx.z*blockDim.z + threadIdx.z;
    iP[dir[1]] += blockIdx.y*blockDim.y + threadIdx.y;
    iP[dir[2]] += blockIdx.x*blockDim.x + threadIdx.x;

    if(iP[dir[0]] < end[dir[0]] && iP[dir[1]] < end[dir[1]] && iP[dir[2]] < end[dir[2]]) {
      function(iP, arguments...);
    }
  }

template <unsigned int Dimension>
  class Computation<Architecture::GPU, Dimension> {
 protected:
  const Position start;
  const Position end;
  const Position length;
  const Position dir;

 public:
  Computation(const Position& start_in,
              const Position& end_in,
              const Position& dir_in = {{d::X, d::Y, d::Z}})
      : start(start_in), end(end_in), length(end_in - start_in), dir(dir_in) {}

  LBM_INLINE static void synchronize() {
    LBM_CUDA_CALL(cudaDeviceSynchronize());
  }

};



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
      LBM_INSTRUMENT_OFF("Computation<Architecture::GPU, 1>::Do<Callback>",3)

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid((127+Base::length[dir[0]])/128, 1, 1);

      kernel_1D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end, Base::dir,
							function, arguments...);
      LBM_CUDA_CALL(cudaGetLastError());
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
      LBM_INSTRUMENT_OFF("Computation<Architecture::GPU, 2>::Do<Callback>",3)

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid((127+Base::length[dir[1]])/128, Base::length[dir[0]], 1);
      kernel_2D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end, Base::dir,
							function, arguments...);
      LBM_CUDA_CALL ( cudaGetLastError() );
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
     LBM_INSTRUMENT_OFF("Computation<Architecture::GPU, 3>::Do<Callback>",3)

     dim3 dimBlock(128, 1, 1);
     dim3 dimGrid((127+Base::length[dir[2]])/128, Base::length[dir[1]], Base::length[dir[0]]);

     kernel_3D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end, Base::dir,
						       function, arguments...);
     LBM_CUDA_CALL(cudaGetLastError());
   }
 };

}
