#include "Input.in"

constexpr ::lbm::Architecture architectureT = ::lbm::Architecture::CPU;

#include "metaLBM/Commons.h"
#include "metaLBM/MPIInitializer.h"
#include "metaLBM/FFTWInitializer.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Routine.h"

int main(int argc, char* argv[]) {
  using namespace lbm;
  LBM_INSTRUMENT_ON("main", 0)

  auto mpiLauncher = MPIInitializer<numProcs, architectureT>{argc, argv};
  auto fftwLauncher = FFTWInitializer<numThreads>{};

  Routine<dataT, algorithmT, architectureT, memoryL, partitionningT,
          communicationT, overlappingT> routine;

  routine.compute();
}
