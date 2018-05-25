#pragma once

#include <fftw3-mpi.h>
#include <mpi.h>
#include "Lattice.h"
#include "MathVector.h"
#include "Domain.h"

namespace lbm {

  /// RAII container for launching FFTW
  template <int numThreadsAtCompileTime>
  struct FFTWInitializer {
    static unsigned int numberElements;

    FFTWInitializer() {
      #ifdef USE_FFTW
      fftw_init_threads();
      fftw_mpi_init();
      fftw_plan_with_nthreads(numThreadsAtCompileTime);

      ptrdiff_t lX_fftw;
      ptrdiff_t startX_fftw;
      numberElements = 2 * fftw_mpi_local_size(L::dimD,
                Cast<unsigned int, ptrdiff_t, 3>::Do(gSD::sLength()).data(),
                MPI_COMM_WORLD, &lX_fftw, &startX_fftw);
      #endif
    }

    /// Finalizes FFTW
    ~FFTWInitializer() {
      #ifdef USE_FFTW
      fftw_mpi_cleanup();
      #endif
    }

  };  // end class FFTWInitializer

  using FFTWInit = FFTWInitializer<numThreads>;

  template<> unsigned int FFTWInit::numberElements = lSD::pVolume();

}  // end namespace lbm
