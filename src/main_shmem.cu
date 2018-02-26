#include <mpi.h>
#include <shmem.h>
#include <shmemx.h>

#include <iostream>
#include <string>

#define NTHREADS 1
#define NPROCS 1
#define _AOS
#define DATA_TYPE double

#include "Input.h"
#include "metaLBM/Lattice.h"
#include "metaLBM/Commons.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Computation.sh"

#ifdef USE_FFTW
  #include <fftw3-mpi.h>
  #include "metaLBM/FourierDomain.h"
#else
  #include "metaLBM/Domain.h"
#endif


#include "metaLBM/Routine.h"
using namespace lbm;

int main(int argc, char* argv[]) {
  INSTRUMENT_ON("main",0)

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

  MPI_Comm comm_MPI;
  shmemx_init_attr_t attribute_SHMEM;
  comm_MPI = MPI_COMM_WORLD;

  attribute_SHMEM.mpi_comm = &comm_MPI;
  shmemx_init_attr (SHMEMX_INIT_WITH_MPI_COMM, &attribute_SHMEM);
  int sizePE_SHMEM = shmem_n_pes();
  int rankPE_SHMEM = shmem_my_pe();

  MathVector<int, 3> size_MPI{1, 1, 1};
  MPI_Comm_size(mpiComm, &size_MPI[d::X]);

  MathVector<int, 3> rank_MPI{0, 0, 0};
  MPI_Comm_rank(mpiComm, &rank_MPI[d::X]);

  int  processorNameLength_MPI;
  char processorName_MPI[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(processorName_MPI, &processorNameLength_MPI);

  MPI_Comm localComm_MPI;
  MPI_Info info_MPI;

  MPI_Info_create(&info_MPI);
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank_MPI[d::X],
  info, &localCommMPI);

  int localRank = 0;
  MPI_Comm_rank(localComm_MPI, &localRank_MPI);

  std::cout << "Local rank: " << localRank_MPI << std::endl;
  int numberDevices = 0;
  CUDA_CALL( cudaGetDeviceCount(&numberDevices); )
  CUDA_CALL( cudaSetDevice(localRank_MPI % numberDevices); )

  MPI_Comm_free(&localComm_MPI);
  MPI_Info_free(&info_MPI);


  Routine<dataT, Architecture::SHMEM, inputOutputType> routine(rank_MPI, size_MPI,
                                              std::string(processorName_MPI));

  routine.compute();

  MPI_Finalize();
  cudaDeviceReset();

  return EXIT_SUCCESS;
}
