#ifndef TRANSFORMER_H
#define TRANSFORMER_H

#include <mpi.h>

#ifdef USE_FFTW
  #include <fftw3-mpi.h>
#endif

#include "Commons.h"
#include "MathVector.h"
#include "Computation.h"
#include "FourierDomain.h"


namespace lbm {

  template<class T, Architecture architecture, PartitionningType partitionningType,
           unsigned int Dimension, unsigned int NumberComponents>
  class ForwardFFT {};

  template<class T, Architecture architecture, PartitionningType partitionningType,
           unsigned int Dimension, unsigned int NumberComponents>
  class BackwardFFT {};


  template<unsigned int Dimension, unsigned int NumberComponents>
  class ForwardFFT<double, Architecture::CPU, PartitionningType::OneD,
                           Dimension, NumberComponents> {
  public:
    double * * localSpacePtr;
    double * * localFourierPtr;

  private:
    fftw_plan planForward[NumberComponents];
    Computation<Architecture::CPU, L::dimD> computationLocal;

  public:
    ForwardFFT(double * * localSpacePtr_in,
               double * * localFourierPtr_in,
               const ptrdiff_t globalLength_in[3])
      : localSpacePtr(localSpacePtr_in)
      , localFourierPtr(localFourierPtr_in)
      , computationLocal(lSD::sStart(), lSD::sEnd())
    {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        planForward[iC] = fftw_mpi_plan_dft_r2c(Dimension,
                                                globalLength_in,
                                                localSpacePtr[iC],
                                                (fftw_complex *) localFourierPtr[iC],
                                                MPI_COMM_WORLD, FFTW_ESTIMATE);
      }
    }

    ForwardFFT(double * * localSpacePtr_in,
               const ptrdiff_t globalLength_in[3])
      : ForwardFFT(localSpacePtr_in, localSpacePtr_in, globalLength_in)
    {}

    ~ForwardFFT()
    {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        fftw_destroy_plan(planForward[iC]);
      }
    }

    HOST
    inline void execute() {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        fftw_execute(planForward[iC]);
        computationLocal.Do([=] HOST (const Position& iP) {
            localSpacePtr[iC][lSD::getIndex(iP)] /= gSD::sVolume();
          });
      }
    }

  };


  template<unsigned int Dimension, unsigned int NumberComponents>
  class BackwardFFT<double, Architecture::CPU, PartitionningType::OneD,
                               Dimension, NumberComponents> {
  public:
    double * * localFourierPtr;
    double * * localSpacePtr;

  private:
    fftw_plan planBackward[NumberComponents];

  public:
    BackwardFFT(double * * localFourierPtr_in,
                double * * localSpacePtr_in,
                const ptrdiff_t globalLength_in[3])
      : localFourierPtr(localFourierPtr_in)
      , localSpacePtr(localSpacePtr_in)
    {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        planBackward[iC] = fftw_mpi_plan_dft_c2r(Dimension,
                                                 globalLength_in,
                                                 (fftw_complex *) localFourierPtr[iC],
                                                 localSpacePtr[iC],
                                                 MPI_COMM_WORLD, FFTW_ESTIMATE);
      }
    }

    BackwardFFT(double * * localSpacePtr_in,
                const ptrdiff_t globalLength_in[3])
      : BackwardFFT(localSpacePtr_in, localSpacePtr_in, globalLength_in)
    {}

    ~BackwardFFT()
    {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        fftw_destroy_plan(planBackward[iC]);
      }
    }

    HOST
    inline void execute() {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        fftw_execute(planBackward[iC]);
      }
    }

  };


  template<class T, Architecture architecture, PartitionningType partitionningType,
           unsigned int Dimension, unsigned int NumberComponents>
  class Curl {};

  template<>
  class Curl<double, Architecture::CPU, PartitionningType::OneD, 2, 2> {
  private:
    ForwardFFT<double, Architecture::CPU, PartitionningType::OneD, 2, 2> forwardIn;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, 2, 2> backwardIn;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, 2, 1> backwardOut;
    const Position& offset;
    Computation<Architecture::CPU, L::dimD> computationFourier;

  public:
    Curl(double * * localSpaceInPtr_in,
         double * * localSpaceOutPtr_in,
         const ptrdiff_t globalLength_in[3],
         const Position& offset_in)
      : forwardIn(localSpaceInPtr_in, localSpaceInPtr_in, globalLength_in)
      , backwardIn(localSpaceInPtr_in, localSpaceInPtr_in, globalLength_in)
      , backwardOut(localSpaceOutPtr_in, localSpaceOutPtr_in, globalLength_in)
      , offset(offset_in)
      , computationFourier(lFD::start(), lFD::end())
    {}

    HOST
    inline void executeFourier() {
      computationFourier.Do
        ([=] HOST (const Position& iFP) {
          auto index = lFD::getIndex(iFP);

          WaveNumber iK{{0}};
          iK[d::X] = iFP[d::X]+offset[d::X] <= gSD::sLength()[d::X]/2 ?
            iFP[d::X]+offset[d::X] : iFP[d::X]+offset[d::X]-gSD::sLength()[d::X];
          iK[d::Y] = iFP[d::Y] <= gSD::sLength()[d::Y]/2 ?
            iFP[d::Y] : iFP[d::Y]-gSD::sLength()[d::Y];

          ((fftw_complex *) (backwardOut.localFourierPtr[0]))[index][p::Re]
            = - iK[d::X] * ((fftw_complex *) (forwardIn.localFourierPtr[d::Y]))[index][p::Im]
            + iK[d::Y] * ((fftw_complex *) (forwardIn.localFourierPtr[d::X]))[index][p::Im];

          ((fftw_complex *) (backwardOut.localFourierPtr[0]))[index][p::Im]
            = iK[d::X] * ((fftw_complex *) (forwardIn.localFourierPtr[d::Y]))[index][p::Re]
            - iK[d::Y] * ((fftw_complex *) (forwardIn.localFourierPtr[d::X]))[index][p::Re];
        });

      backwardOut.execute();
    }

    HOST
    inline void executeSpace() {
      forwardIn.execute();
      executeFourier();
      backwardIn.execute();
    }


  };

  template<>
  class Curl<double, Architecture::CPU, PartitionningType::OneD, 3, 3> {
  private:
    ForwardFFT<double, Architecture::CPU, PartitionningType::OneD, 3, 3> forwardIn;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, 3, 3> backwardIn;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD, 3, 3> backwardOut;

    const Position& offset;
    Computation<Architecture::CPU, L::dimD> computationFourier;

  public:
    Curl(double * * localSpaceInPtr_in,
         double * * localSpaceOutPtr_in,
         const ptrdiff_t globalLength_in[3],
         const Position& offset_in)
      : forwardIn(localSpaceInPtr_in, globalLength_in)
      , backwardIn(localSpaceInPtr_in, globalLength_in)
      , backwardOut(localSpaceOutPtr_in, globalLength_in)
      , offset(offset_in)
      , computationFourier(lFD::start(), lFD::end())
    {}

    Curl(double * * localSpacePtr_in,
         const ptrdiff_t globalLength_in[3],
         const Position& offset_in)
      : Curl(localSpacePtr_in, localSpacePtr_in, globalLength_in, offset_in)
    {}

    HOST
    inline void executeFourier() {
      computationFourier.Do
        ([=] HOST (const Position& iFP) {
          auto index = lFD::getIndex(iFP);

          WaveNumber iK{{0}};
          iK[d::X] = iFP[d::X]+offset[d::X] <= gSD::sLength()[d::X]/2 ?
            iFP[d::X]+offset[d::X] : iFP[d::X]+offset[d::X]-gSD::sLength()[d::X];
          iK[d::Y] = iFP[d::Y] <= gSD::sLength()[d::Y]/2 ?
            iFP[d::Y] : iFP[d::Y]-gSD::sLength()[d::Y];
          iK[d::Z] = iFP[d::Z] <= gSD::sLength()[d::Z]/2 ?
            iFP[d::Z] : iFP[d::Z]-gSD::sLength()[d::Z];

          ((fftw_complex *) (backwardOut.localFourierPtr[d::X]))[index][p::Re]
            = - iK[d::Y] * ((fftw_complex *) (forwardIn.localFourierPtr[d::Z]))[index][p::Im]
            + iK[d::Z] * ((fftw_complex *) (forwardIn.localFourierPtr[d::Y]))[index][p::Im];

          ((fftw_complex *) (backwardOut.localFourierPtr[d::X]))[index][p::Im]
            = iK[d::Y] * ((fftw_complex *) (forwardIn.localFourierPtr[d::Z]))[index][p::Re]
            - iK[d::Z] * ((fftw_complex *) (forwardIn.localFourierPtr[d::Y]))[index][p::Re];


          ((fftw_complex *) (backwardOut.localFourierPtr[d::Y]))[index][p::Re]
            = - iK[d::Z] * ((fftw_complex *) (forwardIn.localFourierPtr[d::X]))[index][p::Im]
            + iK[d::X] * ((fftw_complex *) (forwardIn.localFourierPtr[d::Z]))[index][p::Im];

          ((fftw_complex *) (backwardOut.localFourierPtr[d::Y]))[index][p::Im]
            = iK[d::Z] * ((fftw_complex *) (forwardIn.localFourierPtr[d::X]))[index][p::Re]
            - iK[d::X] * ((fftw_complex *) (forwardIn.localFourierPtr[d::Z]))[index][p::Re];


            ((fftw_complex *) (backwardOut.localFourierPtr[d::Z]))[index][p::Re]
            = - iK[d::X] * ((fftw_complex *) (forwardIn.localFourierPtr[d::Y]))[index][p::Im]
            + iK[d::Y] * ((fftw_complex *) (forwardIn.localFourierPtr[d::X]))[index][p::Im];

          ((fftw_complex *) (backwardOut.localFourierPtr[d::Z]))[index][p::Im]
            = iK[d::X] * ((fftw_complex *) (forwardIn.localFourierPtr[d::Y]))[index][p::Re]
            - iK[d::Y] * ((fftw_complex *) (forwardIn.localFourierPtr[d::X]))[index][p::Re];
        });

      backwardOut.execute();
    }

    HOST
    inline void executeSpace() {
      forwardIn.execute();
      executeFourier();
      backwardIn.execute();
    }

  };



}

#endif // TRANSFORMER_H
