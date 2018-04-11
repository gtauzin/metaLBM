#pragma once

#include <stdio.h>

#include "Computation.h"
#include "Commons.h"
#include "Options.h"
#include "Stream.cuh"

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
    : public Computation<Architecture::GPU, 0> {
  private:
   using Base = Computation<Architecture::GPU, 0>;

  public:
   using Base::Computation;

    template<typename Callback, typename... Arguments>
    void Do(Stream<Architecture::GPU> stream,
            Callback function, const Arguments... arguments) {
      { INSTRUMENT_OFF("Computation<Architecture::GPU, 1>::Do<Callback>",3) }

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid((127+Base::length[d::X])/128, 1, 1);
      kernel_1D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end,
                                                        function, arguments...);
      CUDA_CALL ( cudaGetLastError(); );
      CUDA_CALL( cudaDeviceSynchronize() );
    }

   template<typename Callback, typename... Arguments>
   void Do(Callback function, const Arguments... arguments) {
     Do<Callback>(DefaultStream<Architecture::GPU>(), function, arguments...);
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
    void Do(Stream<Architecture::GPU> stream,
            Callback function, const Arguments... arguments) {
      { INSTRUMENT_OFF("Computation<Architecture::GPU, 2>::Do<Callback>",3) }

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid((127+Base::length[d::Y])/128, Base::length[d::X], 1);
      kernel_2D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end,
                                                        function, arguments...);
      CUDA_CALL ( cudaGetLastError() );
      CUDA_CALL( cudaDeviceSynchronize() );
    }

   template<typename Callback, typename... Arguments>
   void Do(Callback function, const Arguments... arguments) {
     Do<Callback>(DefaultStream<Architecture::GPU>(), function, arguments...);
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
   void Do(Stream<Architecture::GPU>& stream,
           Callback function, const Arguments... arguments) {
     { INSTRUMENT_OFF("Computation<Architecture::GPU, 3>::Do<Callback>",3) }

     dim3 dimBlock(128, 1, 1);
     dim3 dimGrid((127+Base::length[d::Z])/128, Base::length[d::Y], Base::length[d::X]);

     kernel_3D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end,
                                                       function, arguments...);
     CUDA_CALL ( cudaGetLastError() );
     CUDA_CALL( cudaDeviceSynchronize() );
   }

   template<typename Callback, typename... Arguments>
   void Do(Callback function, const Arguments... arguments) {
     Stream<Architecture::GPU> stream = DefaultStream<Architecture::GPU>();	
     Do<Callback>(stream,  function, arguments...);
   }
 };

}
