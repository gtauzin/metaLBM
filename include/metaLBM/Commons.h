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
  #define RESTRICT //__restrict__

  #ifdef USE_NVTX
    #include <nvToolsExt.h>

    class Tracer {
      static constexpr int numberColors = 7;
      const uint32_t colors[numberColors] = { 0x0000ff00, 0x000000ff, 0x00ffff00, 0x00ff00ff,
                                          0x0000ffff, 0x00ff0000, 0x00ffffff };

    public:
      //HOST DEVICE
      Tracer(const char* name, int colorID) {
        nvtxEventAttributes_t eventAttribute = {0};
        eventAttribute.version = NVTX_VERSION;
        eventAttribute.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
        eventAttribute.colorType = NVTX_COLOR_ARGB;
        eventAttribute.color = colors[colorID%numberColors];
        eventAttribute.messageType = NVTX_MESSAGE_TYPE_ASCII;
        eventAttribute.message.ascii = name;
        nvtxRangePushEx(&eventAttribute);
      }

      //HOST DEVICE
      ~Tracer() {
        nvtxRangePop();
      }
    };

    #define INSTRUMENT_ON(name,colorID) Tracer uniq_name_using_macros(name, colorID);

  #endif // USE_NVTX

#else // __NVCC__
  #define CUDA_CALL(call)

  #define HOST
  #define DEVICE //#pragma omp declare simd
  #define CONSTANT
  #define GLOBAL
  #define RESTRICT __restrict__

#endif // __NVXX__


#ifdef USE_SCOREP
  #include <scorep/SCOREP_User.h>

  #define INSTRUMENT_ON(name,colorID) {                                          \
      SCOREP_USER_REGION(name,SCOREP_USER_REGION_TYPE_FUNCTION)                  \
    }

#else // USE_SCOREP
  #ifndef USE_NVTX
    #define INSTRUMENT_ON(name,colorID)

  #endif // USE_NVTX

#endif  // USE_SCOREP

#define INSTRUMENT_OFF(name,colorID)

namespace lbm {

}

#endif // COMMONS_H
