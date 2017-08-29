#ifndef COMMONS_H
#define COMMONS_H

#include <iostream>

#ifdef __CUDA_ARCH__
  #include <cuda_runtime_api.h>
  #include <cuda.h>

  #define CUDA_CALL(call) {
    cudaError_t error = call;
    if(error != cudaSuccess) {
      std::cout << "Cuda failure " << __FILE__ << ":" << __LINE__
                << ": '" << cudaGetErrorString(error) << "'" << std::endl;
      exit(0);
    }
  }

  #define HOST __host__
  #define DEVICE __device__
  #define GLOBAL(function) {
    __global__ function
  }
  #define RESTRICT

#else
  #define CUDA_CALL(call)

  #define HOST
  #define DEVICE
  #define GLOBAL(function)
  #define RESTRICT __restrict__

  #ifdef PROFILE_SCOREP
    #include <scorep/SCOREP_User.h>

    #define SCOREP_INSTRUMENT_ON(function) {
      SCOREP_USER_REGION(function,SCOREP_USER_REGION_TYPE_FUNCTION)
    }

  #else
    #define SCOREP_INSTRUMENT_ON(function)
  #endif

#endif

#define SCOREP_INSTRUMENT_OFF(function)


namespace lbm {

}

#endif // COMMONS_H
