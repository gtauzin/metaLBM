#pragma once

#include <mpi.h>
#include <iostream>

#ifdef USE_NVSHMEM
  #include <shmem.h>
  #include <shmemx.h>
#endif

namespace lbm {

  /// RAII container for launching MPI
  template <int numProcsAtCompileTile>
  struct MPIInitializer {
    static std::string hostName;
    static MathVector<int, 3> size;
    static MathVector<int, 3> rank;
    static int rankLeft;
    static int rankRight;
    static int rankTop;
    static int rankBottom;
    static int rankFront;
    static int rankBack;

    /// Launch MPI
    MPIInitializer(int argc, char** argv) {
      int providedThreadSupport;
      MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &providedThreadSupport);

      int hostNameLength;
      char hostNameChar[MPI_MAX_PROCESSOR_NAME];
      MPI_Get_processor_name(hostNameChar, &hostNameLength);
      hostName = std::string(hostNameChar);

      MPI_Comm_size(MPI_COMM_WORLD, &size[d::X]);

      if (size[d::X] != numProcsAtCompileTile) {
        std::cout << "Compile-time and runtime number of process don't match\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      #ifdef USE_NVSHMEM
        MPI_Comm comm;
        shmemx_init_attr_t attribute;
        comm = MPI_COMM_WORLD;

        attribute.mpi_comm = &comm;
        shmemx_init_attr(SHMEMX_INIT_WITH_MPI_COMM, &attribute);

        rank[d::X] = shmem_my_pe();
      #else
        MPI_Comm_rank(MPI_COMM_WORLD, &rank[d::X]);
      #endif

      rankLeft = (rank[d::X] + size[d::X] - 1) % size[d::X];
      rankRight = (rank[d::X] + 1) % size[d::X];
    }

    /// Finalizes MPI
    ~MPIInitializer() {
      MPI_Finalize();
    }

  };  // end class MPIInitializer

  using MPIInit = MPIInitializer<numProcs>;

  template<> std::string MPIInit::hostName = "";
  template<> MathVector<int, 3> MPIInit::size = MathVector<int, 3>{{0}};
  template<> MathVector<int, 3> MPIInit::rank = MathVector<int, 3>{{0}};
  template<> int MPIInit::rankLeft = 0;
  template<> int MPIInit::rankRight = 0;
  template<> int MPIInit::rankTop = 0;
  template<> int MPIInit::rankBottom = 0;
  template<> int MPIInit::rankFront = 0;
  template<> int MPIInit::rankBack = 0;

}  // end namespace lbm
