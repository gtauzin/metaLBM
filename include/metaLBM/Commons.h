#ifndef COMMONS_H
#define COMMONS_H

#define RESTRICT __restrict__

#ifdef __CUDA_ARCH__

#define cudaCheckError() {                                          \
 cudaError_t e=cudaGetLastError();                                 \
 if(e!=cudaSuccess) {                                              \
   printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
   exit(0); \
 }                                                                 \
}

#define _HOST_ __host__
#define _DEVICE_ __device__

#else

#define _HOST_
#define _DEVICE_

#ifdef PROFILE_SCOREP
#include <scorep/SCOREP_User.h>
#define SCOREP_INSTRUMENT(function) SCOREP_USER_REGION(function,SCOREP_USER_REGION_TYPE_FUNCTION)

#else
#define SCOREP_INSTRUMENT(function)
#endif

#endif


namespace lbm {

}

#endif // COMMONS_H
