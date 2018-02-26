#ifndef ANALYZER_H
#define ANALYZER_H

#include <cmath>

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "Lattice.h"
#include "Communication.h"
#include "Computation.h"

// TODO: Functors

namespace lbm {

  template<class T, AnalyzerType AnalyzerType,
           Architecture architecture, PartitionningType partitionningType>
  class Analyzer {};

  template<class T, Architecture architecture, PartitionningType partitionningType>
  class Analyzer<T, AnalyzerType::Generic, architecture, partitionningType>
  {};

#ifdef USE_FFTW
  template<>
  class Analyzer<double, AnalyzerType::FFT,
                 Architecture::CPU, PartitionningType::Serial>
    : public Analyzer<double, AnalyzerType::Generic,
                      Architecture::CPU, PartitionningType::Serial> {
  public:
    Analyzer()
      : Analyzer<double, AnalyzerType::Generic,
                 Architecture::CPU, PartitionningType::Serial>()
    {}

    HOST
    static inline void forwardFFT(double * RESTRICT spaceGlobalArray,
                                  fftw_complex * RESTRICT fourierGlobalArray) {

      fftw_plan planForward =
        fftw_plan_dft_r2c(L::dimD,
                          Cast<unsigned int, int, 3>::Do(gSD::length()).data(),
                          spaceGlobalArray, fourierGlobalArray,
                          FFTW_ESTIMATE);

      fftw_execute(planForward);
      fftw_destroy_plan(planForward);
    }

    HOST
    static inline void calculateBackwardFFT(fftw_complex * RESTRICT fourierGlobalArray,
                                            double * RESTRICT spaceGlobalArray) {

      fftw_plan planBackward =
        fftw_plan_dft_c2r(L::dimD,
                          Cast<unsigned int, int, 3>::Do(gSD::length()).data(),
                          fourierGlobalArray, spaceGlobalArray,
                          FFTW_ESTIMATE);

      fftw_execute(planBackward);
      fftw_destroy_plan(planBackward);
    }

  };

  template<>
  class Analyzer<double, AnalyzerType::FFT,
                 Architecture::CPU, PartitionningType::OneD>
    : public Analyzer<double, AnalyzerType::FFT,
                      Architecture::CPU, PartitionningType::Serial> {
  public:
    HOST
    static inline void calculateForwardFFT(double * RESTRICT spaceLocalArray,
                                           fftw_complex * RESTRICT fourierLocalArray) {

      fftw_plan planForward =
        fftw_mpi_plan_dft_r2c(L::dimD,
                              Cast<unsigned int, ptrdiff_t, 3>::Do(lSD::length()).data(),
                              spaceLocalArray, fourierLocalArray,
                              MPI_COMM_WORLD, FFTW_ESTIMATE);

      fftw_execute(planForward);
    }

    HOST
    static inline void calculateBackwardFFT(fftw_complex * RESTRICT fourierLocalArray,
                                            double * RESTRICT spaceLocalArray) {

      fftw_plan planBackward =
        fftw_mpi_plan_dft_c2r(L::dimD,
                              Cast<unsigned int, ptrdiff_t, 3>::Do(lSD::length()).data(),
                              fourierLocalArray, spaceLocalArray,
                              MPI_COMM_WORLD, FFTW_ESTIMATE);

      fftw_execute(planBackward);
    }

  };
#endif

}

#endif // ANALYZER_H
