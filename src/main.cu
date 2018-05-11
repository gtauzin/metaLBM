#include "Input.in"
#include "metaLBM/Computation.cuh"
#include "metaLBM/Event.cuh"
#include "metaLBM/CUDAInitializer.h"
#include "metaLBM/MPIInitializer.h"
#include "metaLBM/Commons.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Routine.h"

#ifdef USE_FFTW
  #include "metaLBM/FFTWInitializer.h"
#endif

#ifdef USE_NVSHMEM
  #include <shmem.h>
  #include <shmemx.h>
#endif

int main(int argc, char* argv[]) {
  using namespace lbm;
  { LBM_INSTRUMENT_ON("main",0) }

  auto mpiLauncher = MPIInitializer<numProcs>{argc, argv};
  auto cudaLauncher = CUDAInitializer{};

#ifndef USE_FFTW
  auto numberElements = lSD::pVolume();
#else
  auto fftwLauncher = FFTWInitializer<numThreads>{};
  auto numberElements = fftwLauncher.numElements();
#endif

#ifdef USE_NVSHMEM
  MPI_Comm commMPI;
  shmemx_init_attr_t attributeSHMEM;
  comm_MPI = MPI_COMM_WORLD;

  attributeSHMEM.mpi_comm = &commMPI;
  shmemx_init_attr (SHMEMX_INIT_WITH_MPI_COMM, &attributeSHMEM);
#endif

  auto sizeMPI = MathVector<int, 3>{mpiLauncher.numProcs(), 1, 1};
  auto rankMPI = MathVector<int, 3>{mpiLauncher.procRank(), 0, 0};
  Routine<dataT, Architecture::GPU, implementationT>
    routine(rankMPI, sizeMPI, mpiLauncher.hostName(), numberElements);

  routine.compute();
}
