#pragma once

#include <mpi.h>

#include <fftw3-mpi.h>

#include "Commons.h"
#include "Computation.h"
#include "FourierDomain.h"
#include "MathVector.h"

namespace lbm {

  template <class T, Architecture architecture, PartitionningType partitionningType,
            unsigned int Dimension, unsigned int NumberComponents>
  class ForwardFFT {};

  template <class T, Architecture architecture, PartitionningType partitionningType,
            unsigned int Dimension, unsigned int NumberComponents>
  class BackwardFFT {};

  template <unsigned int Dimension, unsigned int NumberComponents>
  class ForwardFFT<double, Architecture::CPU, PartitionningType::OneD,
                   Dimension, NumberComponents> {
  public:
    double* localSpacePtr;
    double* localFourierPtr;

 private:
    fftw_plan planForward[NumberComponents];

  public:
    ForwardFFT(double* localSpacePtr_in, double* localFourierPtr_in,
               const ptrdiff_t globalLength_in[3])
      : localSpacePtr(localSpacePtr_in), localFourierPtr(localFourierPtr_in)
    {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        planForward[iC] = fftw_mpi_plan_dft_r2c(
          Dimension, globalLength_in, localSpacePtr + FFTWInit::numberElements * iC,
          (fftw_complex*)(localFourierPtr + FFTWInit::numberElements * iC),
          MPI_COMM_WORLD, FFTW_ESTIMATE);
      }
    }

    ForwardFFT(double* localSpacePtr_in, const ptrdiff_t globalLength_in[3])
      : ForwardFFT(localSpacePtr_in, localSpacePtr_in, globalLength_in)
    {}

    ~ForwardFFT() {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        fftw_destroy_plan(planForward[iC]);
      }
    }

    LBM_HOST
    inline void execute() {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        fftw_execute(planForward[iC]);
      }
    }
  };

  template <unsigned int Dimension, unsigned int NumberComponents>
  class BackwardFFT<double, Architecture::CPU, PartitionningType::OneD,
                    Dimension, NumberComponents> {
  public:
    double* localFourierPtr;
    double* localSpacePtr;

  private:
    fftw_plan planBackward[NumberComponents];
    Computation<Architecture::CPU, L::dimD> computationLocal;

  public:
    BackwardFFT(double* localFourierPtr_in, double* localSpacePtr_in,
                const ptrdiff_t globalLength_in[3])
      : localFourierPtr(localFourierPtr_in)
      , localSpacePtr(localSpacePtr_in)
      , computationLocal(lSD::sStart(), lSD::sEnd())
    {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        planBackward[iC] = fftw_mpi_plan_dft_c2r(
          Dimension, globalLength_in,
          (fftw_complex*)(localFourierPtr + FFTWInit::numberElements * iC),
          localSpacePtr + FFTWInit::numberElements * iC, MPI_COMM_WORLD,
          FFTW_ESTIMATE);
      }
    }

    BackwardFFT(double* localSpacePtr_in, const ptrdiff_t globalLength_in[3])
      : BackwardFFT(localSpacePtr_in, localSpacePtr_in, globalLength_in)
    {}

    ~BackwardFFT() {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        fftw_destroy_plan(planBackward[iC]);
      }
    }

    LBM_HOST
    inline void execute() {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        fftw_execute(planBackward[iC]);
        computationLocal.Do([=] LBM_HOST(const Position& iP) {
            (localSpacePtr + FFTWInit::numberElements * iC)[lSD::getIndex(iP)] /=
              gSD::sVolume();
        });
        computationLocal.synchronize();

      }
    }
  };

  template <class T, Architecture architecture, PartitionningType partitionningType,
            unsigned int Dimension, unsigned int NumberComponents>
  class Curl {};

  template <>
  class Curl<double, Architecture::CPU, PartitionningType::OneD, 2, 2> {
  private:
    ForwardFFT<double, Architecture::CPU, PartitionningType::OneD, 2, 2>
    forwardIn;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, 2, 2>
    backwardIn;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, 2, 1>
    backwardOut;

    const Position& offset;
    Computation<Architecture::CPU, L::dimD> computationFourier;

  public:
    Curl(double* localSpaceInPtr_in, double* localSpaceOutPtr_in,
         const ptrdiff_t globalLength_in[3], const Position& offset_in)
      : forwardIn(localSpaceInPtr_in, localSpaceInPtr_in, globalLength_in)
      , backwardIn(localSpaceInPtr_in, localSpaceInPtr_in, globalLength_in)
      , backwardOut(localSpaceOutPtr_in, localSpaceOutPtr_in, globalLength_in)
      , offset(offset_in)
      , computationFourier(lFD::start(), lFD::end())
    {}

    LBM_HOST
    inline void executeFourier() {
      computationFourier.Do([=] LBM_HOST(const Position& iFP) {
          auto index = lFD::getIndex(iFP);

          WaveNumber iK{{0}};
          iK[d::X] = iFP[d::X] + offset[d::X] <= gSD::sLength()[d::X] / 2
            ? iFP[d::X] + offset[d::X]
            : iFP[d::X] + offset[d::X] - gSD::sLength()[d::X];
          iK[d::Y] = iFP[d::Y] <= gSD::sLength()[d::Y] / 2
            ? iFP[d::Y] : iFP[d::Y] - gSD::sLength()[d::Y];

          ((fftw_complex*)(backwardOut.localFourierPtr))[index][p::Re] =
            -iK[d::X] * ((fftw_complex*)(forwardIn.localFourierPtr +
                            FFTWInit::numberElements * (d::Y)))[index][p::Im] +
            iK[d::Y] * ((fftw_complex*)(forwardIn.localFourierPtr +
                           FFTWInit::numberElements * (d::X)))[index][p::Im];

          ((fftw_complex*)(backwardOut.localFourierPtr))[index][p::Im] =
            iK[d::X] * ((fftw_complex*)(forwardIn.localFourierPtr +
                          FFTWInit::numberElements * (d::Y)))[index][p::Re] -
            iK[d::Y] * ((fftw_complex*)(forwardIn.localFourierPtr +
                          FFTWInit::numberElements * (d::X)))[index][p::Re];
      });
      computationFourier.synchronize();

      backwardOut.execute();
    }

    LBM_HOST
    inline void executeSpace() {
      forwardIn.execute();
      executeFourier();
      backwardIn.execute();
    }
  };

  template <>
  class Curl<double, Architecture::CPU, PartitionningType::OneD, 3, 3> {
  private:
    ForwardFFT<double, Architecture::CPU, PartitionningType::OneD, 3, 3>
    forwardIn;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, 3, 3>
    backwardIn;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, 3, 3>
    backwardOut;

    const Position& offset;
    Computation<Architecture::CPU, L::dimD> computationFourier;

  public:
    Curl(double* localSpaceInPtr_in, double* localSpaceOutPtr_in,
         const ptrdiff_t globalLength_in[3], const Position& offset_in)
      : forwardIn(localSpaceInPtr_in, globalLength_in)
      , backwardIn(localSpaceInPtr_in, globalLength_in)
      , backwardOut(localSpaceOutPtr_in, globalLength_in)
      , offset(offset_in)
      , computationFourier(lFD::start(), lFD::end())
    {}

    LBM_HOST
    inline void executeFourier() {
      computationFourier.Do([=] LBM_HOST(const Position& iFP) {
        auto index = lFD::getIndex(iFP);

        WaveNumber iK{{0}};
        iK[d::X] = iFP[d::X] + offset[d::X] <= gSD::sLength()[d::X] / 2
          ? iFP[d::X] + offset[d::X]
          : iFP[d::X] + offset[d::X] - gSD::sLength()[d::X];
        iK[d::Y] = iFP[d::Y] <= gSD::sLength()[d::Y] / 2
            ? iFP[d::Y]
          : iFP[d::Y] - gSD::sLength()[d::Y];
        iK[d::Z] = iFP[d::Z] <= gSD::sLength()[d::Z] / 2
          ? iFP[d::Z]
          : iFP[d::Z] - gSD::sLength()[d::Z];

        ((fftw_complex*)(backwardOut.localFourierPtr +
                         FFTWInit::numberElements * (d::X)))[index][p::Re] =
            -iK[d::Y] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                         FFTWInit::numberElements * (d::Z)))[index][p::Im] +
          iK[d::Z] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                      FFTWInit::numberElements * (d::Y)))[index][p::Im];

        ((fftw_complex*)(backwardOut.localFourierPtr +
                         FFTWInit::numberElements * (d::X)))[index][p::Im] =
          iK[d::Y] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                      FFTWInit::numberElements * (d::Z)))[index][p::Re] -
          iK[d::Z] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                      FFTWInit::numberElements * (d::Y)))[index][p::Re];

        ((fftw_complex*)(backwardOut.localFourierPtr +
                         FFTWInit::numberElements * (d::Y)))[index][p::Re] =
          -iK[d::Z] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                       FFTWInit::numberElements * (d::X)))[index][p::Im] +
          iK[d::X] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                      FFTWInit::numberElements * (d::Z)))[index][p::Im];

        ((fftw_complex*)(backwardOut.localFourierPtr +
                         FFTWInit::numberElements * (d::Y)))[index][p::Im] =
          iK[d::Z] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                      FFTWInit::numberElements * (d::X)))[index][p::Re] -
          iK[d::X] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                      FFTWInit::numberElements * (d::Z)))[index][p::Re];

        ((fftw_complex*)(backwardOut.localFourierPtr +
                         FFTWInit::numberElements * (d::Z)))[index][p::Re] =
          -iK[d::X] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                       FFTWInit::numberElements * (d::Y)))[index][p::Im] +
          iK[d::Y] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                      FFTWInit::numberElements * (d::X)))[index][p::Im];

        ((fftw_complex*)(backwardOut.localFourierPtr +
                         FFTWInit::numberElements * (d::Z)))[index][p::Im] =
          iK[d::X] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                      FFTWInit::numberElements * (d::Y)))[index][p::Re] -
          iK[d::Y] * ((fftw_complex*)(forwardIn.localFourierPtr +
                                      FFTWInit::numberElements * (d::X)))[index][p::Re];
      });
      computationFourier.synchronize();

      backwardOut.execute();
    }

    LBM_HOST
    inline void executeSpace() {
      forwardIn.execute();
      executeFourier();
      backwardIn.execute();
    }
  };

  template <class T, Architecture architecture, PartitionningType partitionningType,
            unsigned int Dimension>
  class MakeIncompressible {};

  template <>
  class MakeIncompressible<double, Architecture::CPU, PartitionningType::OneD, 2> {
  private:
    fftw_complex* localFourierInPtr;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, 2, 2>
    backwardOut;
    const Position& offset;
    Computation<Architecture::CPU, L::dimD> computationFourier;

  public:
    MakeIncompressible(double* localSpaceInPtr_in, double* localSpaceOutPtr_in,
                       const ptrdiff_t globalLength_in[3], const Position& offset_in)
      : localFourierInPtr((fftw_complex*)localSpaceInPtr_in)
      , backwardOut(localSpaceOutPtr_in, localSpaceOutPtr_in, globalLength_in)
      , offset(offset_in)
      , computationFourier(lFD::start(), lFD::end())
    {}

    LBM_HOST
    inline void executeFourier() {
      computationFourier.Do([=] LBM_HOST(const Position& iFP) {
          WaveNumber iK{{0}};
          Position iFP_symmetric{{0}};
          WaveNumber iK_symmetric{{0}};

          auto index = lFD::getIndex(iFP);

          iK[d::X] = iFP[d::X] + offset[d::X] <= gSD::sLength()[d::X] / 2
            ? iFP[d::X] + offset[d::X] : iFP[d::X] + offset[d::X] - gSD::sLength()[d::X];
          iK[d::Y] = iFP[d::Y] <= gSD::sLength()[d::Y] / 2
            ? iFP[d::Y] : iFP[d::Y] - gSD::sLength()[d::Y];

          ((fftw_complex*)(backwardOut.localFourierPtr +
                       FFTWInit::numberElements * (d::X)))[index][p::Re]
            = -iK[d::Y] * localFourierInPtr[index][p::Im];

          ((fftw_complex*)(backwardOut.localFourierPtr +
                           FFTWInit::numberElements * (d::X)))[index][p::Im]
            = iK[d::Y] * localFourierInPtr[index][p::Re];

          ((fftw_complex*)(backwardOut.localFourierPtr +
                       FFTWInit::numberElements * (d::Y)))[index][p::Re]
            = iK[d::X] * localFourierInPtr[index][p::Im];

          ((fftw_complex*)(backwardOut.localFourierPtr +
                           FFTWInit::numberElements * (d::Y)))[index][p::Im]
            = -iK[d::X] * localFourierInPtr[index][p::Re];

          if (iFP[d::Y] == 0) {
            iFP_symmetric[d::X] = gSD::sLength()[d::X] - iFP[d::X];
            if (iFP_symmetric[d::X] < lFD::length()[d::X]) {
              iFP_symmetric[d::Y] = 0;

              iK_symmetric[d::X] =
                iFP_symmetric[d::X] + offset[d::X] <= gSD::sLength()[d::X] / 2
                ? iFP_symmetric[d::X] + offset[d::X]
                : iFP_symmetric[d::X] + offset[d::X] - gSD::sLength()[d::X];
              iK_symmetric[d::Y] =
                iFP_symmetric[d::Y] <= gSD::sLength()[d::Y] / 2
                ? iFP_symmetric[d::Y]
                : iFP_symmetric[d::Y] - gSD::sLength()[d::Y];

              auto index_symmetric = lFD::getIndex(iFP_symmetric);
              ((fftw_complex*)(backwardOut.localFourierPtr +
                               FFTWInit::numberElements * (d::X)))[index_symmetric][p::Re] =
                -iK_symmetric[d::Y] * localFourierInPtr[index_symmetric][p::Im];

              ((fftw_complex*)(backwardOut.localFourierPtr +
                               FFTWInit::numberElements * (d::X)))[index_symmetric][p::Im] =
                iK_symmetric[d::Y] * localFourierInPtr[index_symmetric][p::Re];

              ((fftw_complex*)(backwardOut.localFourierPtr +
                               FFTWInit::numberElements * (d::Y)))[index_symmetric][p::Re] =
            iK_symmetric[d::X] * localFourierInPtr[index_symmetric][p::Im];

              ((fftw_complex*)(backwardOut.localFourierPtr +
                               FFTWInit::numberElements * (d::Y)))[index_symmetric][p::Im] =
                -iK_symmetric[d::X] * localFourierInPtr[index_symmetric][p::Re];
            }
          }
      });
      computationFourier.synchronize();

      backwardOut.execute();
    }
  };

  template <>
  class MakeIncompressible<double, Architecture::CPU, PartitionningType::OneD, 3>
    : public Curl<double, Architecture::CPU, PartitionningType::OneD, 3, 3> {
  private:
    using Base = Curl<double, Architecture::CPU, PartitionningType::OneD, 3, 3>;

  public:
    MakeIncompressible(double* localSpaceInPtr_in, double* localSpaceOutPtr_in,
                       const ptrdiff_t globalLength_in[3], const Position& offset_in)
      : Base(localSpaceInPtr_in, localSpaceOutPtr_in, globalLength_in, offset_in)
    {}

    using Base::executeFourier;
  };

}  // namespace lbm
