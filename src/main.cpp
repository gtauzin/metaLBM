#include <mpi.h>
#include <iostream>
#include <string>

#include "./Input.in"
#include "metaLBM/Commons.h"
#include "metaLBM/Computation.h"
#include "metaLBM/Lattice.h"
#include "metaLBM/MathVector.h"
#include "metaLBM/Routine.h"

#ifdef USE_FFTW
#include <fftw3-mpi.h>
#endif

#include "metaLBM/FourierDomain.h"

using namespace lbm;

int main(int argc, char* argv[]) {
  INSTRUMENT_ON("main", 0)

#ifndef USE_FFTW
  MPI_Init(&argc, &argv);
  {
    unsigned int numberElements = lSD::pVolume();

#else
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  {
    int useThreadsFFTW = (provided >= MPI_THREAD_FUNNELED);

    if (useThreadsFFTW)
      useThreadsFFTW = fftw_init_threads();
    fftw_mpi_init();

    if (useThreadsFFTW)
      fftw_plan_with_nthreads(NTHREADS);

    ptrdiff_t lX_fftw;
    ptrdiff_t startX_fftw;

    unsigned int numberElements =
        (unsigned int)2 *
        fftw_mpi_local_size(
            L::dimD,
            Cast<unsigned int, ptrdiff_t, 3>::Do(gSD::sLength()).data(),
            MPI_COMM_WORLD, &lX_fftw, &startX_fftw);
#endif

    MathVector<int, 3> sizeMPI = {1, 1, 1};
    MPI_Comm_size(MPI_COMM_WORLD, &sizeMPI[d::X]);
    if (sizeMPI[d::X] != numProcs) {
      std::cout << "Compile-time and runtime number of process don't match\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MathVector<int, 3> rankMPI{0, 0, 0};
    MPI_Comm_rank(MPI_COMM_WORLD, &rankMPI[d::X]);

    int hostnameLength;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(hostname, &hostnameLength);

    Routine<dataT, Architecture::CPU, implementationT> routine(
        rankMPI, sizeMPI, std::string(hostname), numberElements);

    routine.compute();

#ifdef USE_FFTW
    fftw_mpi_cleanup();
#endif
  }
  MPI_Finalize();

  return EXIT_SUCCESS;
}
