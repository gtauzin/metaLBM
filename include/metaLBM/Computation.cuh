#ifndef COMPUTATION_CUH
#define COMPUTATION_CUH

#include "Computation.h"
#include "Commons.h"
#include "Options.h"
#include "Domain.h"

namespace lbm {

#ifdef __CUDA_ARCH__

  template<typename Callback>
  GLOBAL
  void kernel_1D(const MathVector<unsigned int, 3>& start,
                 const MathVector<unsigned int, 3>& end,
                 Callback function) {

    MathVector<unsigned int, 3> iP = {blockIdx.x*blockDim.x + threadIdx.x,
                                      start[d::Y], start[d::Z]};

    if(iP[0] >= start[d::X] && iP[0] < end[d::X]) {
        function(iP);
    }
  }

  template<typename Callback>
  GLOBAL
  void kernel_2D(const MathVector<unsigned int, 3>& start,
              const MathVector<unsigned int, 3>& end,
              Callback function) {

    MathVector<unsigned int, 3> iP = {blockIdx.y*blockDim.y + threadIdx.y,
                                      blockIdx.x*blockDim.x + threadIdx.x,
                                      start[d::Z]};

    if(iP[1] >= start[d::X] && iP[1] < end[d::X]
       && iP[0] >= start[d::Y] && iP[0] < end[d::Y]) {
      function(iP);
    }
  }

  template<typename Callback>
  GLOBAL
  void kernel_3D(const MathVector<unsigned int, 3>& start,
                 const MathVector<unsigned int, 3>& end,
                 Callback function) {

    MathVector<unsigned int, 3> iP = {blockIdx.z*blockDim.z + threadIdx.z,
                                      blockIdx.y*blockDim.y + threadIdx.y,
                                      blockIdx.x*blockDim.x + threadIdx.x};

    if(iP[2] >= start[d::X] && iP[2] < end[d::X]
       && iP[1] >= start[d::Y] && iP[1] < end[d::Y]
       && iP[0] >= start[d::Z] && iP[0] < end[d::Z]) {
      function(iP);
    }
  }


  template<>
  struct Computation<Architecture::GPU, 1> {
    template<typename Callback>
    HOST
    static void Do(const MathVector<unsigned int, 3>& start,
                   const MathVector<unsigned int, 3>& end,
                   Callback function) {
      SCOREP_INSTRUMENT_OFF("Computation<Architecture::GPU, 1>::Do<Callback>")

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid(hD::length()[d::X]/128, 1, 1);
      //CUDA_CALL(
      kernel_1D<Callback><<<dimBlock, dimGrid >>>(start, end, function);
      //)

      CUDA_CALL( cudaDeviceSynchronize(); )
    }
  };

  template<>
  struct Computation<Architecture::GPU, 2> {
  public:
    template<typename Callback>
    HOST
    static void Do(const MathVector<unsigned int, 3>& start,
                   const MathVector<unsigned int, 3>& end,
                   Callback function) {
      SCOREP_INSTRUMENT_OFF("Computation<Architecture::GPU, 2>::Do<Callback>")

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid(hD::length()[d::Y]/128, hD::length()[d::X], 1);
      //      CUDA_CALL(
      kernel_2D<Callback><<<dimBlock, dimGrid>>>(start, end, function);
                //)

      CUDA_CALL( cudaDeviceSynchronize(); )
    }
  };


  template<>
  struct Computation<Architecture::GPU, 3> {
  public:
    template<typename Callback>
    HOST
    static void Do(const MathVector<unsigned int, 3>& start,
                   const MathVector<unsigned int, 3>& end,
                   Callback function) {
      SCOREP_INSTRUMENT_OFF("Computation<Architecture::GPU, 3>::Do<Callback>")

      dim3 dimBlock(128, 1, 1);
      dim3 dimGrid(hD::length()[d::Z]/128, hD::length()[d::Y], hD::length()[d::X]);
      //CUDA_CALL(
      kernel_3D<Callback><<<dimBlock, dimGrid>>>(start, end, function);
                //)

      CUDA_CALL( cudaDeviceSynchronize(); )
    }
  };

#endif

  typedef Computation<arch, L::dimD> Computation_;

}

#endif // COMPUTATION_CUH
