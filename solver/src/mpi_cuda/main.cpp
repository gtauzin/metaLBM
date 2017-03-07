#include <mpi.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "commons.h"
#include "input.h"
#include "init.h"
#include "compute.h"

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

using namespace lbm;

int main() {
  std::cout << lengthX_g << std::endl;
  //  static_assert(lengthX_g%NPROCS == 0, "lengthX_g must be a multiple of NPROCS");

  MPI_Init(NULL, NULL);

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  initLogging(mpi_rank);

  BOOST_LOG_TRIVIAL(debug) << "Logging for debug starts.";

  Init init = init_Simulation(mpi_rank);

  compute(init, mpi_rank);

  MPI_Finalize();

  return EXIT_SUCCESS;
}
