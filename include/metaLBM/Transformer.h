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
    double* spacePtr;
    double* fourierPtr;

 private:
    fftw_plan planForward[NumberComponents];

  public:
    ForwardFFT(double* spacePtr_in, double* fourierPtr_in,
               const ptrdiff_t globalLength_in[3])
      : spacePtr(spacePtr_in), fourierPtr(fourierPtr_in)
    {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        planForward[iC] = fftw_mpi_plan_dft_r2c(Dimension, globalLength_in,
                                                spacePtr + FFTWInit::numberElements * iC,
                                                (fftw_complex*)(fourierPtr + FFTWInit::numberElements * iC),
                                                MPI_COMM_WORLD, FFTW_ESTIMATE);
      }
    }

    ForwardFFT(double* spacePtr_in, const ptrdiff_t globalLength_in[3])
      : ForwardFFT(spacePtr_in, spacePtr_in, globalLength_in)
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
    double* fourierPtr;
    double* spacePtr;

  private:
    fftw_plan planBackward[NumberComponents];
    Computation<Architecture::CPU, L::dimD> computationLocal;

  public:
    BackwardFFT(double* fourierPtr_in, double* spacePtr_in,
                const ptrdiff_t globalLength_in[3])
      : fourierPtr(fourierPtr_in)
      , spacePtr(spacePtr_in)
      , computationLocal(lSD::sStart(), lSD::sEnd())
    {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        planBackward[iC] = fftw_mpi_plan_dft_c2r(Dimension, globalLength_in,
                                                 (fftw_complex*)(fourierPtr + FFTWInit::numberElements * iC),
                                                 spacePtr + FFTWInit::numberElements * iC, MPI_COMM_WORLD,
                                                 FFTW_ESTIMATE);
      }
    }

    BackwardFFT(double* spacePtr_in, const ptrdiff_t globalLength_in[3])
      : BackwardFFT(spacePtr_in, spacePtr_in, globalLength_in)
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
            (spacePtr + FFTWInit::numberElements * iC)[lSD::getIndex(iP)] /=
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

    const Position offset;
    Computation<Architecture::CPU, 2> computationFourier;
    Computation<Architecture::CPU, 2> computationLocal;

  public:
    Curl(double* spaceInPtr_in, double* spaceOutPtr_in,
         const ptrdiff_t globalLength_in[3], const Position& offset_in)
      : forwardIn(spaceInPtr_in, spaceInPtr_in, globalLength_in)
      , backwardIn(spaceInPtr_in, spaceInPtr_in, globalLength_in)
      , backwardOut(spaceOutPtr_in, spaceOutPtr_in, globalLength_in)
      , offset(offset_in)
      , computationFourier(lFD::start(), lFD::end())
      , computationLocal(lSD::sStart(), lSD::sEnd())
    {}

    LBM_HOST
    inline void executeFourier() {
      computationFourier.Do([=] LBM_HOST(const Position& iFP) {
          auto index = lFD::getIndex(iFP);

          WaveNumber iK{{0}};

          iK[d::X] = iFP[d::X] + offset[d::X] <= gSD::sLength()[d::X] / 2
            ? iFP[d::X] + offset[d::X] : iFP[d::X] + offset[d::X] - gSD::sLength()[d::X];
          iK[d::Y] = iFP[d::Y] <= gSD::sLength()[d::Y] / 2
            ? iFP[d::Y] : iFP[d::Y] - gSD::sLength()[d::Y];

          ((fftw_complex*)(backwardOut.fourierPtr))[index][p::Re] =
            - iK[d::X] * ((fftw_complex*)(forwardIn.fourierPtr
                           + FFTWInit::numberElements * (d::Y)))[index][p::Im]
            + iK[d::Y] * ((fftw_complex*)(forwardIn.fourierPtr
                           + FFTWInit::numberElements * (d::X)))[index][p::Im];

          ((fftw_complex*)(backwardOut.fourierPtr))[index][p::Im] =
            + 0*iK[d::X] * ((fftw_complex*)(forwardIn.fourierPtr
                           + FFTWInit::numberElements * (d::Y)))[index][p::Re]
            - 0*iK[d::Y] * ((fftw_complex*)(forwardIn.fourierPtr
                           + FFTWInit::numberElements * (d::X)))[index][p::Re];
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

    LBM_HOST
    inline void normalize() {
      computationLocal.Do([=] LBM_HOST(const Position& iP) {
          backwardOut.fourierPtr[lSD::getIndex(iP)] /= gSD::sVolume();
      });
      computationLocal.synchronize();
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
    Computation<Architecture::CPU, 3> computationFourier;
    Computation<Architecture::CPU, 3> computationLocal;

  public:
    Curl(double* spaceInPtr_in, double* spaceOutPtr_in,
         const ptrdiff_t globalLength_in[3], const Position& offset_in)
      : forwardIn(spaceInPtr_in, globalLength_in)
      , backwardIn(spaceInPtr_in, globalLength_in)
      , backwardOut(spaceOutPtr_in, globalLength_in)
      , offset(offset_in)
      , computationFourier(lFD::start(), lFD::end())
      , computationLocal(lSD::sStart(), lSD::sEnd())
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

        ((fftw_complex*)(backwardOut.fourierPtr +
                         FFTWInit::numberElements * (d::X)))[index][p::Re] =
            -iK[d::Y] * ((fftw_complex*)(forwardIn.fourierPtr +
                                         FFTWInit::numberElements * (d::Z)))[index][p::Im] +
          iK[d::Z] * ((fftw_complex*)(forwardIn.fourierPtr +
                                      FFTWInit::numberElements * (d::Y)))[index][p::Im];

        ((fftw_complex*)(backwardOut.fourierPtr +
                         FFTWInit::numberElements * (d::X)))[index][p::Im] =
          iK[d::Y] * ((fftw_complex*)(forwardIn.fourierPtr +
                                      FFTWInit::numberElements * (d::Z)))[index][p::Re] -
          iK[d::Z] * ((fftw_complex*)(forwardIn.fourierPtr +
                                      FFTWInit::numberElements * (d::Y)))[index][p::Re];

        ((fftw_complex*)(backwardOut.fourierPtr +
                         FFTWInit::numberElements * (d::Y)))[index][p::Re] =
          -iK[d::Z] * ((fftw_complex*)(forwardIn.fourierPtr +
                                       FFTWInit::numberElements * (d::X)))[index][p::Im] +
          iK[d::X] * ((fftw_complex*)(forwardIn.fourierPtr +
                                      FFTWInit::numberElements * (d::Z)))[index][p::Im];

        ((fftw_complex*)(backwardOut.fourierPtr +
                         FFTWInit::numberElements * (d::Y)))[index][p::Im] =
          iK[d::Z] * ((fftw_complex*)(forwardIn.fourierPtr +
                                      FFTWInit::numberElements * (d::X)))[index][p::Re] -
          iK[d::X] * ((fftw_complex*)(forwardIn.fourierPtr +
                                      FFTWInit::numberElements * (d::Z)))[index][p::Re];

        ((fftw_complex*)(backwardOut.fourierPtr +
                         FFTWInit::numberElements * (d::Z)))[index][p::Re] =
          -iK[d::X] * ((fftw_complex*)(forwardIn.fourierPtr +
                                       FFTWInit::numberElements * (d::Y)))[index][p::Im] +
          iK[d::Y] * ((fftw_complex*)(forwardIn.fourierPtr +
                                      FFTWInit::numberElements * (d::X)))[index][p::Im];

        ((fftw_complex*)(backwardOut.fourierPtr +
                         FFTWInit::numberElements * (d::Z)))[index][p::Im] =
          iK[d::X] * ((fftw_complex*)(forwardIn.fourierPtr +
                                      FFTWInit::numberElements * (d::Y)))[index][p::Re] -
          iK[d::Y] * ((fftw_complex*)(forwardIn.fourierPtr +
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

    LBM_HOST
    inline void normalize() {
      const unsigned int numberElements = FFTWInit::numberElements;
      computationLocal.Do([=] LBM_HOST(const Position& iP) {
          for (auto iC = 0; iC < 3; ++iC) {
            (backwardOut.fourierPtr + numberElements * iC)[lSD::getIndex(iP)] /=
              gSD::sVolume();
          }
        });
      computationLocal.synchronize();
    }
  };

  template <class T, Architecture architecture, PartitionningType partitionningType,
            unsigned int Dimension>
  class MakeIncompressible {};

  template <>
  class MakeIncompressible<double, Architecture::CPU, PartitionningType::OneD, 2> {
  private:
    fftw_complex* fourierInPtr;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, 2, 2>
    backwardOut;
    const Position& offset;
    Computation<Architecture::CPU, L::dimD> computationFourier;

  public:
    MakeIncompressible(double* spaceInPtr_in, double* spaceOutPtr_in,
                       const ptrdiff_t globalLength_in[3], const Position& offset_in)
      : fourierInPtr((fftw_complex*)spaceInPtr_in)
      , backwardOut(spaceOutPtr_in, spaceOutPtr_in, globalLength_in)
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

          ((fftw_complex*)(backwardOut.fourierPtr +
                       FFTWInit::numberElements * (d::X)))[index][p::Re]
            = -iK[d::Y] * fourierInPtr[index][p::Im];

          ((fftw_complex*)(backwardOut.fourierPtr +
                           FFTWInit::numberElements * (d::X)))[index][p::Im]
            = iK[d::Y] * fourierInPtr[index][p::Re];

          ((fftw_complex*)(backwardOut.fourierPtr +
                       FFTWInit::numberElements * (d::Y)))[index][p::Re]
            = iK[d::X] * fourierInPtr[index][p::Im];

          ((fftw_complex*)(backwardOut.fourierPtr +
                           FFTWInit::numberElements * (d::Y)))[index][p::Im]
            = -iK[d::X] * fourierInPtr[index][p::Re];

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
              ((fftw_complex*)(backwardOut.fourierPtr +
                               FFTWInit::numberElements * (d::X)))[index_symmetric][p::Re] =
                -iK_symmetric[d::Y] * fourierInPtr[index_symmetric][p::Im];

              ((fftw_complex*)(backwardOut.fourierPtr +
                               FFTWInit::numberElements * (d::X)))[index_symmetric][p::Im] =
                iK_symmetric[d::Y] * fourierInPtr[index_symmetric][p::Re];

              ((fftw_complex*)(backwardOut.fourierPtr +
                               FFTWInit::numberElements * (d::Y)))[index_symmetric][p::Re] =
            iK_symmetric[d::X] * fourierInPtr[index_symmetric][p::Im];

              ((fftw_complex*)(backwardOut.fourierPtr +
                               FFTWInit::numberElements * (d::Y)))[index_symmetric][p::Im] =
                -iK_symmetric[d::X] * fourierInPtr[index_symmetric][p::Re];
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
    MakeIncompressible(double* spaceInPtr_in, double* spaceOutPtr_in,
                       const ptrdiff_t globalLength_in[3], const Position& offset_in)
      : Base(spaceInPtr_in, spaceOutPtr_in, globalLength_in, offset_in)
    {}

    using Base::executeFourier;
  };


  template <class T, Architecture architecture, PartitionningType partitionningType,
            unsigned int Dimension>
  class StrainDiagonalCalculator {};


  template <unsigned int Dimension>
  class StrainDiagonalCalculator<double, Architecture::CPU, PartitionningType::OneD, Dimension> {
  private:
    ForwardFFT<double, Architecture::CPU, PartitionningType::OneD, Dimension, Dimension>
    forwardVelocity;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, Dimension, Dimension>
    backwardDiagonal;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, Dimension, Dimension>
    backwardVelocity;

    const Position& offset;
    Computation<Architecture::CPU, 3> computationFourier;
    Computation<Architecture::CPU, 3> computationLocal;

  public:
    StrainDiagonalCalculator(double* velocityPtr_in, double* strainDiagonalPtr_in,
                             const ptrdiff_t globalLength_in[3], const Position& offset_in)
      : forwardVelocity(velocityPtr_in, velocityPtr_in, globalLength_in)
      , backwardDiagonal(strainDiagonalPtr_in, strainDiagonalPtr_in, globalLength_in)
      , backwardVelocity(velocityPtr_in, velocityPtr_in, globalLength_in)
      , offset(offset_in)
      , computationFourier(lFD::start(), lFD::end())
      , computationLocal(lSD::sStart(), lSD::sEnd())
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

        for(unsigned int iD = 0; iD < Dimension; ++iD) {
          ((fftw_complex*)(backwardDiagonal.fourierPtr +
                           FFTWInit::numberElements * (iD)))[index][p::Re] =
            -iK[iD] * ((fftw_complex*)(forwardVelocity.fourierPtr +
                                       FFTWInit::numberElements * (iD)))[index][p::Im];

          ((fftw_complex*)(backwardDiagonal.fourierPtr +
                           FFTWInit::numberElements * (iD)))[index][p::Im] =
            iK[iD] * ((fftw_complex*)(forwardVelocity.fourierPtr +
                                      FFTWInit::numberElements * (iD)))[index][p::Re];
        }
      });
      computationFourier.synchronize();
      //normalize();
      backwardDiagonal.execute();
    }

    LBM_HOST
    inline void executeSpace() {
      forwardVelocity.execute();
      executeFourier();
      backwardVelocity.execute();
    }

    LBM_HOST
    inline void normalize() {
      const unsigned int numberElements = FFTWInit::numberElements;
      computationLocal.Do([=] LBM_HOST(const Position& iP) {
          for (auto iC = 0; iC < Dimension; ++iC) {
            (backwardDiagonal.fourierPtr + numberElements * iC)[lSD::getIndex(iP)] *=
              gSD::sVolume();
          }
        });
      computationLocal.synchronize();
    }
  };


  template <class T, Architecture architecture, PartitionningType partitionningType,
            unsigned int Dimension>
  class StrainSymmetricCalculator {};


  template <unsigned int Dimension>
  class StrainSymmetricCalculator<double, Architecture::CPU, PartitionningType::OneD, Dimension> {
  private:
    ForwardFFT<double, Architecture::CPU, PartitionningType::OneD, Dimension, Dimension>
      forwardVelocity;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, Dimension, 2*Dimension-3>
      backwardSymmetric;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, Dimension, Dimension>
      backwardVelocity;

    const Position& offset;
    Computation<Architecture::CPU, 3> computationFourier;
    Computation<Architecture::CPU, 3> computationLocal;

  public:
    StrainSymmetricCalculator(double* velocityPtr_in, double* strainSymmetricPtr_in,
                              const ptrdiff_t globalLength_in[3], const Position& offset_in)
      : forwardVelocity(velocityPtr_in, velocityPtr_in, globalLength_in)
      , backwardSymmetric(strainSymmetricPtr_in, strainSymmetricPtr_in, globalLength_in)
      , backwardVelocity(velocityPtr_in, velocityPtr_in, globalLength_in)
      , offset(offset_in)
      , computationFourier(lFD::start(), lFD::end())
      , computationLocal(lSD::sStart(), lSD::sEnd())
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

        ((fftw_complex*)(backwardSymmetric.fourierPtr +
                         FFTWInit::numberElements * (d::X)))[index][p::Re] =
          0.5 * (-iK[d::X] * ((fftw_complex*)(forwardVelocity.fourierPtr +
                                              FFTWInit::numberElements * (d::Y)))[index][p::Im]
                 -iK[d::Y] * ((fftw_complex*)(forwardVelocity.fourierPtr +
                                              FFTWInit::numberElements * (d::X)))[index][p::Im]);

        ((fftw_complex*)(backwardSymmetric.fourierPtr +
                         FFTWInit::numberElements * (d::X)))[index][p::Im] =
          0.5 * (iK[d::X] * ((fftw_complex*)(forwardVelocity.fourierPtr +
                                             FFTWInit::numberElements * (d::Y)))[index][p::Re] +
                 iK[d::Y] * ((fftw_complex*)(forwardVelocity.fourierPtr +
                                              FFTWInit::numberElements * (d::X)))[index][p::Re]);

        for(unsigned int iD = 1; iD < 2 * Dimension - 3; ++iD) {
          ((fftw_complex*)(backwardSymmetric.fourierPtr +
                           FFTWInit::numberElements * (d::X)))[index][p::Re] =
            0.5 * (-iK[iD-1] * ((fftw_complex*)(forwardVelocity.fourierPtr +
                                                FFTWInit::numberElements * (d::Z)))[index][p::Im]
                   -iK[d::Z] * ((fftw_complex*)(forwardVelocity.fourierPtr +
                                                FFTWInit::numberElements * (iD-1)))[index][p::Im]);

          ((fftw_complex*)(backwardSymmetric.fourierPtr +
                           FFTWInit::numberElements * (iD)))[index][p::Im] =
            0.5 * (iK[iD-1] * ((fftw_complex*)(forwardVelocity.fourierPtr +
                                                FFTWInit::numberElements * (d::Z)))[index][p::Re] +
                    iK[d::Z] * ((fftw_complex*)(forwardVelocity.fourierPtr +
                                                FFTWInit::numberElements * (iD-1)))[index][p::Re]);
        }

      });
      computationFourier.synchronize();
      //normalize();
     backwardSymmetric.execute();
    }

    LBM_HOST
    inline void executeSpace() {
      forwardVelocity.execute();
      executeFourier();
      backwardVelocity.execute();
    }

    LBM_HOST
    inline void normalize() {
      const unsigned int numberElements = FFTWInit::numberElements;
      computationLocal.Do([=] LBM_HOST(const Position& iP) {
          for (auto iC = 0; iC < 2*Dimension-3; ++iC) {
            (backwardSymmetric.fourierPtr + numberElements * iC)[lSD::getIndex(iP)] *=
              gSD::sVolume();
          }
        });
      computationLocal.synchronize();
    }
  };


}  // namespace lbm
