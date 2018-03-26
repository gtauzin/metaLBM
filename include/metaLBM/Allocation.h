#ifndef ALLOCATION_H
#define ALLOCATION_H

#include "Options.h"
#include "Helpers.h"
#include "Commons.h"
#include "Lattice.h"
#include "MathVector.h"

#ifdef USE_FFTW
  #include <fftw3-mpi.h>
#endif

namespace alloc {
  constexpr int globalLengthInt[3] = {::lbm::globalLengthX,
                                     ::lbm::L::dimD>1 ? ::lbm::globalLengthY: 1,
                                     ::lbm::L::dimD>2 ? ::lbm::globalLengthZ: 1};

  constexpr ptrdiff_t globalLengthPtrdiff_t[3] = {::lbm::globalLengthX,
                                                 ::lbm::L::dimD>1 ? ::lbm::globalLengthY: 1,
                                                 ::lbm::L::dimD>2 ? ::lbm::globalLengthZ: 1};
  HOST
  int getNumberElements() {
    #ifdef USE_FFTW
    ptrdiff_t lX_fftw;
    ptrdiff_t startX_fftw;

    return 2*fftw_mpi_local_size(::lbm::L::dimD,
                                 globalLengthPtrdiff_t,
                                 MPI_COMM_WORLD,
                                 &lX_fftw, &startX_fftw);
    #else
    return ::lbm::lSD::pVolume();
    #endif
  }

  const int numberElements = getNumberElements();

} // namespace alloc

#endif // ALLOCATION_H
