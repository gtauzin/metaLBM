#include <mpi.h>

#include "Input.in"
#include "metaLBM/Commons.h"
#include "metaLBM/Computation.h"
#include "metaLBM/Lattice.h"
#include "metaLBM/MPIInitializer.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Routine.h"

#ifdef USE_FFTW
#include <fftw3-mpi.h>
#include "metaLBM/FFTWInitializer.h"
#endif

#include "metaLBM/FourierDomain.h"

int main(int argc, char* argv[]) {
  using namespace lbm;
  INSTRUMENT_ON("main", 0)

  MPIInitializer<numProcs> mpiLauncher{argc, argv};
#ifndef USE_FFTW
  std::size_t numberElements = lSD::pVolume();
#else
  FFTWInitializer<numThreads> fftwLauncher{};
  ptrdiff_t lX_fftw;
  ptrdiff_t startX_fftw;
  std::size_t numberElements =
      2 * fftw_mpi_local_size(
              L::dimD,
              Cast<unsigned int, ptrdiff_t, 3>::Do(gSD::sLength()).data(),
              MPI_COMM_WORLD, &lX_fftw, &startX_fftw);
#endif
  auto sizeMPI = MathVector<int, 3>{mpiLauncher.numProcs(), 1, 1};
  auto rankMPI = MathVector<int, 3>{mpiLauncher.procRank(), 0, 0};
  auto routine = Routine<dataT, Architecture::CPU, implementationT>(
      rankMPI, sizeMPI, mpiLauncher.hostName(), numberElements);
  routine.compute();
}
