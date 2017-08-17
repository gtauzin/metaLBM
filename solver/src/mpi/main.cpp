#include <mpi.h>
#include <iostream>
#include <string>

#include "metaLBM/Routine.h"
#include "metaLBM/MathVector.h"

using namespace lbm;

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  MathVector<int, 3> sizeMPI = {1, 1, 1};
  MPI_Comm_size(MPI_COMM_WORLD, &sizeMPI[d::X]);

  MathVector<int, 3> rankMPI{0, 0, 0};
  MPI_Comm_rank(MPI_COMM_WORLD, &rankMPI[d::X]);

  int  hostnameLength;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(hostname, &hostnameLength);

  Routine_ routine(rankMPI, sizeMPI,
                   std::string(hostname));

  routine.compute();

  MPI_Finalize();

  return EXIT_SUCCESS;
}
