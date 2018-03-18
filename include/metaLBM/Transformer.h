#ifndef TRANSFORMER_H
#define TRANSFORMER_H

#include <mpi.h>

#ifdef USE_FFTW
  #include <fftw3-mpi.h>
#endif

#include "Commons.h"
#include "MathVector.h"
#include "Computation.h"
#include "Domain.h"


namespace lbm {

  template<class T, Architecture architecture, PartitionningType partitionningType,
           unsigned int Dimension, unsigned int NumberComponents>
  class FourierTransformer {};

  template<unsigned int Dimension, unsigned int NumberComponents>
  class FourierTransformer<double, Architecture::CPU, PartitionningType::OneD,
                           Dimension, NumberComponents> {
  private:
    fftw_plan planForward[NumberComponents];
    fftw_plan planBackward[NumberComponents];
    Computation<Architecture::CPU, L::dimD> localComputation;
  public:
    FourierTransformer(double * localPtr_in[NumberComponents],
                       const ptrdiff_t globalLength_in[3])
      : localComputation(lSD::sStart(), lSD::sEnd())
    {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        planForward[iC] = fftw_mpi_plan_dft_r2c(Dimension,
                                                globalLength_in,
                                                (double *) localPtr_in[iC],
                                                (fftw_complex *) localPtr_in[iC],
                                                MPI_COMM_WORLD, FFTW_ESTIMATE);
        planBackward[iC] = fftw_mpi_plan_dft_c2r(Dimension,
                                                 globalLength_in,
                                                 (fftw_complex *) localPtr_in[iC],
                                                 (double *) localPtr_in[iC],
                                                 MPI_COMM_WORLD, FFTW_ESTIMATE);
      }
    }

    ~FourierTransformer()
    {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        fftw_destroy_plan(planForward[iC]);
        fftw_destroy_plan(planBackward[iC]);
      }
    }

    HOST
    inline void executeForwardFFT(double * localPtr[NumberComponents]) {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        fftw_execute(planForward[iC]);
        localComputation.Do([=] HOST (const Position& iP) {
            localPtr[iC][lSD::getIndex(iP)] /= gSD::sVolume();
          });
      }
    }

    HOST
    inline void executeBackwardFFT() {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        fftw_execute(planBackward[iC]);
      }
    }

  };

}

#endif // TRANSFORMER_H
