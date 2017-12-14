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
#include "metaLBM/Computation.h"
#include "metaLBM/Computation.cuh"


#include "metaLBM/Routine.h"
using namespace lbm;

int main(int argc, char* argv[]) {
  INSTRUMENT_ON("main",0)

  MPI_Init(&argc, &argv);

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

  //MPI_Info_create(&info_MPI);
  //MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank_MPI[d::X],
  //info, &localCommMPI);

  int localRank = 0;

//MPI_Comm_rank(localComm_MPI, &localRank_MPI);

  std::cout << "Local rank: " << localRank_MPI << std::endl;
  int numberDevices = 0;
  CUDA_CALL( cudaGetDeviceCount(&numberDevices); )
  CUDA_CALL( cudaSetDevice(localRank_MPI % numberDevices); )

  MPI_Comm_free(&localComm_MPI);
  MPI_Info_free(&info_MPI);


  Routine<dataT, Architecture::SHMEM> routine(rank_MPI, size_MPI,
                                              std::string(processorName_MPI));

  routine.compute();

  MPI_Finalize();

  return EXIT_SUCCESS;
}
