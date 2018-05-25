#include "Input.in"
#include "metaLBM/Commons.h"
#include "metaLBM/MPIInitializer.h"
#include "metaLBM/FFTWInitializer.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Routine.h"

int main(int argc, char* argv[]) {
  using namespace lbm;
  LBM_INSTRUMENT_ON("main", 0)

  auto mpiLauncher = MPIInitializer<numProcs>{argc, argv};
  auto fftwLauncher = FFTWInitializer<numThreads>{};

  Routine<dataT, Architecture::CPU, implementationT> routine;

  routine.compute();
}
