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

#define CACHE_LINE 64

namespace lbm {

}

#endif // COMMONS_H
