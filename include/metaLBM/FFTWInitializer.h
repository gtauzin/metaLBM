#ifndef METALBM_FFTWINITIALIZER_H
#define METALBM_FFTWINITIALIZER_H

#include <fftw3-mpi.h>

namespace lbm {

/// RAII container for launching FFTW
template<int numThreadsAtCompileTime>
struct FFTWInitializer {

  FFTWInitializer() {
      fftw_init_threads();
      fftw_mpi_init();
      fftw_plan_with_nthreads(numThreadsAtCompileTime);
  }

  FFTWInitializer(FFTWInitializer const&) = delete;

  FFTWInitializer& operator=(FFTWInitializer const&) = delete;

  FFTWInitializer(FFTWInitializer&&) = delete;

  FFTWInitializer& operator=(FFTWInitializer&&) = delete;

  /// Finalizes FFTW
  ~FFTWInitializer() { fftw_mpi_cleanup(); }

};  // end class FFTWInitializer

}  // end namespace lbm

#endif  // METALBM_FFTWINITIALIZER_H
