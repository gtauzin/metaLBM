#define NTHREADS 1
#define NPROCS 1
#define _AOS
#define DATA_TYPE double

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

#include "Input.h"
#include "metaLBM/Lattice.h"
#include "metaLBM/Commons.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Computation.cuh"
#include "metaLBM/DynamicArray.cuh"

namespace lbm {
  constexpr Architecture arch = Architecture::GPU;
}

#include "metaLBM/Routine.h"
using namespace lbm;

int main(int argc, char* argv[]) {
  INSTRUMENT_ON("main",0)

  MPI_Init(&argc, &argv);
  #ifndef USE_FFTW
    MPI_Init(&argc, &argv);
  #else
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int useThreadsFFTW = (provided >= MPI_THREAD_FUNNELED);

    if(useThreadsFFTW) useThreadsFFTW = fftw_init_threads();
    fftw_mpi_init();

    if (useThreadsFFTW) fftw_plan_with_nthreads(NTHREADS);

    ptrdiff_t lengthX_l_fftw;
    ptrdiff_t startX_l_fftw;
    ptrdiff_t lengthX_g_fftw =
    fftw_mpi_local_size_many(L::dimD,
                             Cast<unsigned int, ptrdiff_t, 3>::Do(gSD::length()).data(),
                             1, FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD,
                             &lengthX_l_fftw, &startX_l_fftw);

    std::cout << "Initializing FFTW" << std::endl;
    std::cout << "-- global length x: " << lengthX_g_fftw << std::endl;
    std::cout << "-- local length x: " << lengthX_l_fftw << std::endl;
  #endif

  #ifdef USE_NVSHMEM
  MPI_Comm comm_MPI;
  shmemx_init_attr_t attribute_SHMEM;
  comm_MPI = MPI_COMM_WORLD;

  attribute_SHMEM.mpi_comm = &comm_MPI;
  shmemx_init_attr (SHMEMX_INIT_WITH_MPI_COMM, &attribute_SHMEM);
  #endif

  MathVector<int, 3> sizeMPI{1, 1, 1};
  MPI_Comm_size(MPI_COMM_WORLD, &sizeMPI[d::X]);

  MathVector<int, 3> rankMPI{0, 0, 0};
  MPI_Comm_rank(MPI_COMM_WORLD, &rankMPI[d::X]);

  int  hostnameLength;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(hostname, &hostnameLength);

  MPI_Comm localComm;
  MPI_Info info;

  MPI_Info_create(&info);
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rankMPI[d::X],
                      info, &localComm);

  int localRank = 0;

  MPI_Comm_rank(localComm, &localRank);

  std::cout << "Local rank: " << localRank << std::endl;
  int numberDevices = 0;
  CUDA_CALL( cudaGetDeviceCount(&numberDevices); )
  CUDA_CALL( cudaSetDevice(localRank % numberDevices); )

  MPI_Comm_free(&localComm);
  MPI_Info_free(&info);

  Routine<dataT, arch, implementationT,
          inputOutputType> routine(rankMPI, sizeMPI,
                                   std::string(hostname));

  routine.compute();

  MPI_Finalize();
  cudaDeviceReset();

  return EXIT_SUCCESS;
}
