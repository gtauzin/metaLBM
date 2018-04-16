#include "Input.in"
#include "metaLBM/Commons.h"
#include "metaLBM/MPIInitializer.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Routine.h"

#ifdef USE_FFTW
#include "metaLBM/FFTWInitializer.h"
#endif

int main(int argc, char* argv[]) {
  using namespace lbm;
  LBM_SCOREP_INSTRUMENT_ON("main", 0)

  auto mpiLauncher = MPIInitializer<numProcs>{argc, argv};
#ifndef USE_FFTW
  unsigned int numberElements = lSD::pVolume();
#else
  auto fftwLauncher = FFTWInitializer<numThreads>{};
  auto numberElements = fftwLauncher.numElements();
#endif
  auto sizeMPI = MathVector<int, 3>{mpiLauncher.numProcs(), 1, 1};
  auto rankMPI = MathVector<int, 3>{mpiLauncher.procRank(), 0, 0};
  auto routine = Routine<dataT, Architecture::CPU, implementationT>(
      rankMPI, sizeMPI, mpiLauncher.hostName(), numberElements);

  routine.compute();
}
