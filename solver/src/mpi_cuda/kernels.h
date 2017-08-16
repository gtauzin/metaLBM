#ifndef KERNELS_CUH
#define KERNELS_CUH

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace lbm {

  extern "C++"
  void dummyPropagateDevice(double * nxt,
                            double * prv,
                            const  int startX,
                            const  int   endX,
                            const  int startY,
                            const  int   endY);

}

#endif
