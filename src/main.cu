#include <mpi.h>
#include <iostream>
#include <string>

#ifdef USE_FFTW
  #include <fftw3-mpi.h>
#endif

#ifdef USE_NVSHMEM
  #include <shmem.h>
  #include <shmemx.h>
#endif

#include "./Input.in"
#include "metaLBM/Lattice.h"
#include "metaLBM/Commons.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Computation.cuh"
#include "metaLBM/DynamicArray.cuh"
#include "metaLBM/Routine.h"

using namespace lbm;

int main(int argc, char* argv[]) {
  { INSTRUMENT_ON("main",0) }
  {
  #ifndef USE_FFTW
    MPI_Init(&argc, &argv);
    {
      unsigned int numberElements = lSD::pVolume();

  #else
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    {
      int useThreadsFFTW = (provided >= MPI_THREAD_FUNNELED);

      if(useThreadsFFTW) useThreadsFFTW = fftw_init_threads();
      fftw_mpi_init();

      if (useThreadsFFTW) fftw_plan_with_nthreads(NTHREADS);

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

  MathVector<int, 3> sizeMPI{1, 1, 1};
  MPI_Comm_size(MPI_COMM_WORLD, &sizeMPI[d::X]);
  if (sizeMPI[d::X] != numProcs) {
    std::cout << "Compile-time and runtime number of process don't match\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MathVector<int, 3> rankMPI{0, 0, 0};
  MPI_Comm_rank(MPI_COMM_WORLD, &rankMPI[d::X]);

  int  hostnameLength;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(hostname, &hostnameLength);

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

  // MPI_Comm_free(&localComm);
  // MPI_Info_free(&info);

  Routine<dataT, Architecture::GPU, implementationT> routine(rankMPI, sizeMPI,
                                                             std::string(hostname),
                                                             numberElements);

  routine.compute();

  #ifdef USE_FFTW
    fftw_mpi_cleanup();
  #endif

    }
    MPI_Finalize();
    }

    cudaDeviceReset();

    return EXIT_SUCCESS;
}
