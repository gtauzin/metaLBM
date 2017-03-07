extern "C" {
#include "utils_gpu.h"
}

#include <cuda_runtime_api.h>
#include <cuda.h>

namespace lbm {

  extern "C"
  void Init_DeviceLattice(double *f, double *d_f, const int size) {

    cudaErrorCheck(cudaMalloc((void **) &d_f, size));
    cudaErrorCheck(cudaMemcpy(d_f, f, size, cudaMemcpyHostToDevice));
  }

  extern "C"
  void Free_DeviceLattice(double *d_f){
    cudaErrorCheck(cudaFree(d_f));
  }

  extern "C"
  void Copy_LatticeFromDevice(double *d_f, double *f, const int size) {
    cudaErrorCheck(cudaMemcpy(f, d_f, size, cudaMemcpyDeviceToHost));
  }

  extern "C"
  void Copy_LatticeToDevice(double *d_f, double *f, const int size) {
    cudaErrorCheck(cudaMemcpy(d_f, f, size, cudaMemcpyHostToDevice));
  }

}
