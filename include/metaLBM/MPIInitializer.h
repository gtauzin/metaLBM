#pragma once

#include <mpi.h>
#include <iostream>

#ifdef USE_NVSHMEM
  #include <shmem.h>
  #include <shmemx.h>
#endif

namespace lbm {

  /// RAII container for launching MPI
  template <int numProcsAtCompileTile, Architecture architecture>
  class MPIInitializer {};

  template <int numProcsAtCompileTile>
  class MPIInitializer<numProcsAtCompileTile, Architecture::Generic> {
  public:
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

      MPI_Comm_rank(MPI_COMM_WORLD, &rank[d::X]);

      rankLeft = (rank[d::X] + size[d::X] - 1) % size[d::X];
      rankRight = (rank[d::X] + 1) % size[d::X];
    }

    /// Finalizes MPI
    ~MPIInitializer() {
      MPI_Finalize();
    }

  };

  template <int numProcsAtCompileTile>
  class MPIInitializer<numProcsAtCompileTile, Architecture::CPU>
    : public MPIInitializer<numProcsAtCompileTile, Architecture::Generic> {
  private:
    using Base = MPIInitializer<numProcsAtCompileTile, Architecture::Generic>;

  public:
    using Base::MPIInitializer;
  };

  template <int numProcsAtCompileTile>
  class MPIInitializer<numProcsAtCompileTile, Architecture::GPU>
    : public MPIInitializer<numProcsAtCompileTile, Architecture::Generic> {
  private:
    using Base = MPIInitializer<numProcsAtCompileTile, Architecture::Generic>;

  public:
    using Base::MPIInitializer;

  };

#ifdef USE_NVSHMEM
  template <int numProcsAtCompileTile>
  class MPIInitializer<numProcsAtCompileTile, Architecture::GPU_SHMEM>
    : public MPIInitializer<numProcsAtCompileTile, Architecture::Generic> {
  private:
    using Base = MPIInitializer<numProcsAtCompileTile, Architecture::Generic>;

  public:
    /// Launch MPI
    MPIInitializer(int argc, char** argv)
      : Base(argc, argv)
    {
      MPI_Comm comm;
      shmemx_init_attr_t attribute;
      comm = MPI_COMM_WORLD;

      attribute.mpi_comm = &comm;
      shmemx_init_attr(SHMEMX_INIT_WITH_MPI_COMM, &attribute);

      Base::rank[d::X] = shmem_my_pe();
      Base::rankLeft = (Base::rank[d::X] + Base::size[d::X] - 1) % Base::size[d::X];
      Base::rankRight = (Base::rank[d::X] + 1) % Base::size[d::X];
    }

    /// Finalizes MPI
    ~MPIInitializer() {
      shmem_finalize();
    }

  };

#endif // USE_NVSHMEM

  using MPIInit = MPIInitializer<numProcs, architectureT>;

  using MPIInitGeneric = MPIInitializer<numProcs, Architecture::Generic>;
  template<> std::string MPIInitGeneric::hostName = "";
  template<> MathVector<int, 3> MPIInitGeneric::size = MathVector<int, 3>{{0}};
  template<> MathVector<int, 3> MPIInitGeneric::rank = MathVector<int, 3>{{0}};
  template<> int MPIInitGeneric::rankLeft = 0;
  template<> int MPIInitGeneric::rankRight = 0;
  template<> int MPIInitGeneric::rankTop = 0;
  template<> int MPIInitGeneric::rankBottom = 0;
  template<> int MPIInitGeneric::rankFront = 0;
  template<> int MPIInitGeneric::rankBack = 0;

}  // end namespace lbm
