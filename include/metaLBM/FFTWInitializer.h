#ifndef METALBM_FFTWINITIALIZER_H
#define METALBM_FFTWINITIALIZER_H

#include <fftw3-mpi.h>
#include <mpi.h>
#include "MathVector.h"
#include "Lattice.h"

namespace lbm {

/// RAII container for launching FFTW
template<int numThreadsAtCompileTime>
struct FFTWInitializer {

  FFTWInitializer() {
      fftw_init_threads();
      fftw_mpi_init();
      fftw_plan_with_nthreads(numThreadsAtCompileTime);
  }

  /// FFTW assumes control over the number of elements assigned to an MPI node
  unsigned int numElements() const noexcept {
      ptrdiff_t lX_fftw;
      ptrdiff_t startX_fftw;
      unsigned int numberElements =
          2 * fftw_mpi_local_size(L::dimD, Cast<unsigned int, ptrdiff_t, 3>::Do(gSD::sLength()).data(),
                                  MPI_COMM_WORLD, &lX_fftw, &startX_fftw);
      return numberElements;
  }

  /// Finalizes FFTW
  ~FFTWInitializer() { fftw_mpi_cleanup(); }

};  // end class FFTWInitializer

}  // end namespace lbm

#endif  // METALBM_FFTWINITIALIZER_H
