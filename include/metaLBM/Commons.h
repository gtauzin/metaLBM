#ifndef COMMONS_H
#define COMMONS_H

#define RESTRICT __restrict__

#ifdef __CUDACC__

#define _HOST_ __host__
#define _DEVICE_ __device__

#else

#define _HOST_
#define _DEVICE_

#endif

#ifdef PROFILE_SCOREP
#include <scorep/SCOREP_User.h>
#define SCOREP_INSTRUMENT(function) SCOREP_USER_REGION(function,SCOREP_USER_REGION_TYPE_FUNCTION)

#else
#define SCOREP_INSTRUMENT(function)
#endif

#define CACHE_LINE 64

namespace lbm {

}

#endif // COMMONS_H
