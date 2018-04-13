#ifndef METALBM_MPIINITIALIZER_H
#define METALBM_MPIINITIALIZER_H

#include <mpi.h>
#include <iostream>

namespace lbm {

/// RAII container for launching MPI

template <int numProcsAtCompileTile>
struct MPIInitializer {
  /// Launch MPI
  MPIInitializer(int argc, char** argv) {
    int providedThreadSupport;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &providedThreadSupport);

    if (numProcs() != numProcsAtCompileTile) {
      std::cout << "Compile-time and runtime number of process don't match\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  int numProcs() const noexcept {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
  }

  int procRank() const noexcept {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
  }

  std::string hostName() const noexcept {
    int hostnameLength;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(hostname, &hostnameLength);
    return std::string(hostname);
  }

  /// Finalizes MPI
  ~MPIInitializer() { MPI_Finalize(); }

};  // end class MPIInitializer

}  // end namespace lbm

#endif  // METALBM_MPIINITIALIZER_H
