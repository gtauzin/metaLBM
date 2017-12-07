#define NTHREADS 1
#define NPROCS 1
#define _AOS
#define DATA_TYPE double

#include <mpi.h>
#include <iostream>
#include <string>

#include "Input.h"

namespace lbm {
  constexpr Architecture arch = Architecture::CPU;
}

#include "metaLBM/Commons.h"
#include "metaLBM/Routine.h"
#include "metaLBM/MathVector.h"

using namespace lbm;

int main(int argc, char* argv[]) {
  //  SCOREP_INSTRUMENT_ON("main")

  MPI_Init(&argc, &argv);

  MathVector<int, 3> sizeMPI = {1, 1, 1};
  MPI_Comm_size(MPI_COMM_WORLD, &sizeMPI[d::X]);

  MathVector<int, 3> rankMPI{0, 0, 0};
  MPI_Comm_rank(MPI_COMM_WORLD, &rankMPI[d::X]);

  int  hostnameLength;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(hostname, &hostnameLength);

  Routine<dataT, arch> routine(rankMPI, sizeMPI,
                                            std::string(hostname));

  routine.compute();

  MPI_Finalize();

  return EXIT_SUCCESS;
}