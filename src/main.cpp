#define NTHREADS 2
#define NPROCS 1
#define _AOS
#define DATA_TYPE double

#include <mpi.h>
#include <iostream>
#include <string>

#include "Input.h"
#include "metaLBM/Lattice.h"
#include "metaLBM/Commons.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Computation.h"

#ifdef USE_FFTW
  #include <fftw3-mpi.h>
  #include "metaLBM/FourierDomain.h"
#else
  #include "metaLBM/Domain.h"
#endif


namespace lbm {
  constexpr Architecture arch = Architecture::CPU;

  typedef Computation<arch, L::dimD> Computation_;
}

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

  MathVector<int, 3> sizeMPI = {1, 1, 1};
  MPI_Comm_size(MPI_COMM_WORLD, &sizeMPI[d::X]);

  MathVector<int, 3> rankMPI{0, 0, 0};
  MPI_Comm_rank(MPI_COMM_WORLD, &rankMPI[d::X]);

  int  hostnameLength;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(hostname, &hostnameLength);

  Routine<dataT, arch, inputOutputType> routine(rankMPI, sizeMPI,
                                                std::string(hostname));

  routine.compute();

  #ifdef USE_FFTW
    fftw_mpi_cleanup();
  #endif
  MPI_Finalize();

  return EXIT_SUCCESS;
}
