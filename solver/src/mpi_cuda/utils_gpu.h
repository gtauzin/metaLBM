#ifndef UTILS_GPU_CUH
#define UTILS_GPU_CUH

#include <iostream>
#include <string>
#include <cuda_runtime_api.h>

#define BLOCK_SIZE 32
#define cudaErrorCheck(ans){ cudaAssert((ans), __FILE__, __LINE__); }

namespace lbm {

  inline void cudaAssert(cudaError_t code, char *file, int line) {
    if (code != cudaSuccess){
      std::cout << "CUDA Assert: " << cudaGetErrorString(code) << ", "
                << std::string(file) << ", " << line << std::endl;
      exit(code);
    }
  }


  //extern "C"
    void Init_DeviceLattice(double *f, double *d_f, const int size);

  //extern "C"
    void Free_DeviceLattice(double *d_f);

  //extern "C"
    void Copy_LatticeFromDevice(double *d_f, double *f, const int size);

  //extern "C"
    void Copy_LatticeToDevice(double *d_f, double *f, const int size);

}

#endif
