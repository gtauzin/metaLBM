#include <mpi.h>
#include <iostream>
#include <string>

#define NTHREADS 1
#define NPROCS 1
#define _AOS
#define USE_PTX
#define DATA_TYPE double


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




  Routine<dataT, arch> routine(rankMPI, sizeMPI,
                               std::string(hostname));

  routine.compute();

  MPI_Finalize();

  return EXIT_SUCCESS;
}
