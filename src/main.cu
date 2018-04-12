#include <mpi.h>
#include <iostream>
#include <string>

#ifdef USE_FFTW
  #include <fftw3-mpi.h>
  #include "metaLBM/FFTWInitializer.h"
#endif

#ifdef USE_NVSHMEM
  #include <shmem.h>
  #include <shmemx.h>
#endif

#include "Input.in"
#include "metaLBM/MPIInitializer.h"
#include "metaLBM/Lattice.h"
#include "metaLBM/Commons.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Stream.cuh"
#include "metaLBM/Computation.cuh"
#include "metaLBM/DynamicArray.cuh"
#include "metaLBM/Routine.h"

int main(int argc, char* argv[]) {
  using namespace lbm;
  { INSTRUMENT_ON("main",0) }
  
  MPIInitializer<numProcs> mpiLauncher {argc, argv};
  #ifndef USE_FFTW
      unsigned int numberElements = lSD::pVolume();

  #else
      FFTWInitializer<numThreads> fftwLauncher{};
      ptrdiff_t lX_fftw;
      ptrdiff_t startX_fftw;
      unsigned int numberElements
        = (unsigned int) 2*fftw_mpi_local_size(L::dimD,
                                               Cast<unsigned int,
                                               ptrdiff_t,
                                               3>::Do(gSD::sLength()).data(),
                                               MPI_COMM_WORLD,
                                               &lX_fftw, &startX_fftw);

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


  // MPI_Comm localComm;
  // MPI_Info info;

  // MPI_Info_create(&info);
  // MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rankMPI[d::X],
  //                     info, &localComm);

  // int localRank = 0;

  // MPI_Comm_rank(localComm, &localRank);
  int localRank = 0;
  int numberDevices = 0;
  CUDA_CALL( cudaGetDeviceCount(&numberDevices); )
  CUDA_CALL( cudaSetDevice(localRank % numberDevices); )

  std::cout << "Number of GPUs available: " << numberDevices << std::endl;

  // MPI_Comm_free(&localComm);
  // MPI_Info_free(&info);

  Routine<dataT, Architecture::GPU, implementationT> routine(rankMPI, sizeMPI,
                                                             mpiLauncher.hostName(),
                                                             numberElements);

  routine.compute();    
  cudaDeviceReset();
}
