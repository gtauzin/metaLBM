#ifndef COMMONS_H
#define COMMONS_H

#include <iostream>

#ifdef __NVCC__
#include <cuda_runtime_api.h>
#include <cuda.h>

  #define CUDA_CALL(call) {                                                     \
    cudaError_t error = call;                                                   \
    if(error != cudaSuccess) {                                                  \
      std::cerr << "Cuda failure " << __FILE__ << ":" << __LINE__               \
                << ": '" << cudaGetErrorString(error) << "'" << std::endl;      \
      exit(0);                                                                  \
    }                                                                           \
  }

  #define HOST __host__
  #define DEVICE __device__
  #define CONSTANT __constant__
  #define GLOBAL __global__
  #define RESTRICT

  #define SCOREP_INSTRUMENT_ON(function)

#else
  #define CUDA_CALL(call)

  #define HOST
  #define DEVICE
  #define CONSTANT
  #define GLOBAL
  #define RESTRICT __restrict__

  #ifdef PROFILE_SCOREP
    #include <scorep/SCOREP_User.h>

    #define SCOREP_INSTRUMENT_ON(function) {                                    \
      SCOREP_USER_REGION(function,SCOREP_USER_REGION_TYPE_FUNCTION)             \
    }

  #else
    #define SCOREP_INSTRUMENT_ON(function)
  #endif

#endif

#define SCOREP_INSTRUMENT_OFF(function)


namespace lbm {

}

#endif // COMMONS_H
