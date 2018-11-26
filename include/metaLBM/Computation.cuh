#pragma once

#include <stdio.h>

#include "Computation.h"
#include "Commons.h"
#include "Options.h"
#include "Stream.cuh"

#ifdef USE_NVSHMEM
  #include<shmem.h>
  #include<shmemx.h>
#endif

namespace lbm {

  template<typename Callback, typename... Arguments>
  LBM_GLOBAL
  void kernel_1D(const Position start, const Position end, const Position order,
                 Callback function, const Arguments... arguments) {
    Position iP = start;
    iP[order[0]] += blockIdx.x*blockDim.x + threadIdx.x;

    // delete if and make sure that this we choose blocks so that we are always
    // in the right case
    if(iP[order[0]] < end[order[0]]) {
      function(iP, arguments...);
    }
  }

  template<typename Callback, typename... Arguments>
  LBM_GLOBAL
  void kernel_2D(const Position start, const Position end, const Position order,
                 Callback function, const Arguments... arguments) {
    Position iP = start;
    iP[order[0]] += blockIdx.y*blockDim.y + threadIdx.y;
    iP[order[1]] += blockIdx.x*blockDim.x + threadIdx.x;

    if(iP[order[0]] < end[order[0]] && iP[order[1]] < end[order[1]]) {
      function(iP, arguments...);
    }
  }

  template<typename Callback, typename... Arguments>
  LBM_GLOBAL
  void kernel_3D(const Position start, const Position end, const Position order,
                 Callback function, const Arguments... arguments) {
    Position iP = start;
    iP[order[0]] += blockIdx.z*blockDim.z + threadIdx.z;
    iP[order[1]] += blockIdx.y*blockDim.y + threadIdx.y;
    iP[order[2]] += blockIdx.x*blockDim.x + threadIdx.x;

    if(iP[order[0]] < end[order[0]] && iP[order[1]] < end[order[1]] && iP[order[2]] < end[order[2]]) {
      function(iP, arguments...);
    }
  }

template <unsigned int Dimension>
  class Computation<Architecture::GPU, Dimension> {
 protected:
  const Position start;
  const Position end;
  const Position length;
  const Position order;

 public:
  Computation(const Position& start_in,
              const Position& end_in,
              const Position& order_in = {{d::X, d::Y, d::Z}})
      : start(start_in), end(end_in), length(end_in - start_in), order(order_in) {}

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
      dim3 dimGrid((127+Base::length[order[0]])/128, 1, 1);

      kernel_1D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end, Base::order,
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
      dim3 dimGrid((127+Base::length[order[1]])/128, Base::length[order[0]], 1);
      kernel_2D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end, Base::order,
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
     dim3 dimGrid((127+Base::length[order[2]])/128, Base::length[order[1]], Base::length[order[0]]);

     kernel_3D<<<dimGrid, dimBlock, 0, stream.get()>>>(Base::start, Base::end, Base::order,
						       function, arguments...);
     LBM_CUDA_CALL(cudaGetLastError());
   }
 };


#ifdef USE_NVSHMEM

 template<>
  class Computation<Architecture::GPU_SHMEM, 1>
    : public Computation<Architecture::GPU, 0> {
  private:
   using Base = Computation<Architecture::GPU, 0>;

  public:
   using Base::Computation;

   using Base::start;
   using Base::end;
   using Base::order;

   template<typename Callback, typename... Arguments>
   void Do(const Stream<Architecture::GPU_SHMEM>& stream,
	    Callback function, const Arguments... arguments) {
      LBM_INSTRUMENT_OFF("Computation<Architecture::GPU_SHMEM, 1>::Do<Callback>",3)

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid((127+Base::length[order[0]])/128, 1, 1);

      void *argumentArray[] = {(void *) &start, (void *) &end, (void *) &order, &function, (void *) &arguments...};
      LBM_SHMEM_CALL( shmemx_collective_launch((const void *) kernel_1D<Callback, Arguments...>, dimGrid, dimBlock,
                                               argumentArray, 0, stream.get()) );
      LBM_CUDA_CALL(cudaGetLastError());
    }
  };

 template<>
  class Computation<Architecture::GPU_SHMEM, 2>
    : public Computation<Architecture::GPU, 0> {
  private:
   using Base = Computation<Architecture::GPU, 0>;

   using Base::start;
   using Base::end;
   using Base::order;

  public:
   using Base::Computation;

    template<typename Callback, typename... Arguments>
    void Do(const Stream<Architecture::GPU_SHMEM>& stream,
            Callback function, const Arguments... arguments) {
      LBM_INSTRUMENT_OFF("Computation<Architecture::GPU_SHMEM, 2>::Do<Callback>",3)

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid((127+Base::length[order[1]])/128, Base::length[order[0]], 1);

      void *argumentArray[] = {(void *) &start, (void *) &end, (void *) &order, &function, (void *) &arguments...};
      LBM_SHMEM_CALL( shmemx_collective_launch((const void *) kernel_2D<Callback, Arguments...>, dimGrid, dimBlock,
                                               argumentArray, 0, stream.get()) );
      LBM_CUDA_CALL ( cudaGetLastError() );
    }
 };

  template<>
  class Computation<Architecture::GPU_SHMEM, 3>
    : public Computation<Architecture::GPU, 0> {
  private:
    using Base = Computation<Architecture::GPU, 0>;

    using Base::start;
    using Base::end;
    using Base::order;

  public:
    using Base::Computation;

    template<typename Callback, typename... Arguments>
    void Do(const Stream<Architecture::GPU_SHMEM>& stream,
	    Callback function, const Arguments... arguments) {
     LBM_INSTRUMENT_OFF("Computation<Architecture::GPU_SHMEM, 3>::Do<Callback>",3)

     dim3 dimBlock(128, 1, 1);
     dim3 dimGrid((127+Base::length[order[2]])/128, Base::length[order[1]], Base::length[order[0]]);

      void *argumentArray[] = {(void *) &start, (void *) &end, (void *) &order, &function, (void *) &arguments...};
      LBM_SHMEM_CALL( shmemx_collective_launch((const void *) kernel_3D<Callback, Arguments...>, dimGrid, dimBlock,
                                               argumentArray, 0, stream.get()) );
     LBM_CUDA_CALL(cudaGetLastError());
   }
 };

#endif // USE_NVSHMEM

}
