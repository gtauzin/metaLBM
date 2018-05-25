#include "Input.in"
#include "metaLBM/Computation.cuh"
#include "metaLBM/Event.cuh"
#include "metaLBM/CUDAInitializer.h"
#include "metaLBM/MPIInitializer.h"
#include "metaLBM/FFTWInitializer.h"
#include "metaLBM/Commons.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Routine.h"


int main(int argc, char* argv[]) {
  using namespace lbm;
  LBM_INSTRUMENT_ON("main",0)

  auto mpiLauncher = MPIInitializer<numProcs>{argc, argv};
  auto cudaLauncher = CUDAInitializer{};
  auto fftwLauncher = FFTWInitializer<numThreads>{};

  Routine<dataT, Architecture::GPU, implementationT> routine;

  routine.compute();
}
