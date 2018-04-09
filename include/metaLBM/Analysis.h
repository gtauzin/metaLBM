#pragma once

#include <cmath>

#ifdef USE_FFTW
  #include <fftw3-mpi.h>
#endif

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "Lattice.h"
#include "Communication.h"
#include "Computation.h"

// TODO: Functors

namespace lbm {

  template<class T>
  class AnalysisScalar {
  public:
    T scalar;
    const unsigned int numberElements;

  protected:
    AnalysisScalar(const unsigned int numberElements_in)
      : scalar((T) 0)
      , numberElements(numberElements_in)
    {}

    inline void reset() {
      scalar = (T) 0;
    }

    inline void normalize() {
      scalar /= gSD::sVolume();
    }
  };

  template<class T>
  class TotalEnergy
    : public AnalysisScalar<T> {
  private:
    using Base = AnalysisScalar<T>;

    T& scalarRef;
    T * localDensityPtr;
    T * localVelocityPtr;

  public:
    static constexpr auto analysisName = "total_energy";

    TotalEnergy(T * localDensityPtr_in, T * localVelocityPtr_in,
                const unsigned int numberElements_in)
      : Base(numberElements_in)
      , scalarRef(Base::scalar)
      , localDensityPtr(localDensityPtr_in)
      , localVelocityPtr(localVelocityPtr_in)
    {}

    HOST
    void operator()(const Position& iP) {
      auto index = lSD::getIndex(iP);
      for(auto iD = 0; iD < L::dimD; ++iD) {
        scalarRef
          += 0.5*localDensityPtr[index]
          *(localVelocityPtr+iD*Base::numberElements)[index]
          *(localVelocityPtr+iD*Base::numberElements)[index];
      }
    }

    using Base::reset;
    using Base::normalize;
    using Base::scalar;
  };


  template<class T>
  class TotalEnstrophy
    : public AnalysisScalar<T> {
  private:
    using Base = AnalysisScalar<T>;

    T * localVorticityPtr;
    T& scalarRef;

  public:
    static constexpr auto analysisName = "total_enstrophy";

    TotalEnstrophy(T * localVorticityPtr_in,
                   const unsigned int numberElements_in)
      : Base(numberElements_in)
      , scalarRef(Base::scalar)
      , localVorticityPtr(localVorticityPtr_in)
    {}

    HOST
    void operator()(const Position& iP) {
      auto index = lSD::getIndex(iP);
      for(auto iD = 0; iD < 2*L::dimD-3; ++iD) {
        scalarRef
          += 0.5*(localVorticityPtr+iD*Base::numberElements)[index]
          *(localVorticityPtr+iD*Base::numberElements)[index];
      }
    }

    using Base::reset;
    using Base::normalize;
    using Base::scalar;
  };


  template<class T, unsigned int maxWaveNumber>
  class AnalysisSpectral {
  public:
    T spectra[maxWaveNumber];
    const unsigned int numberElements;

  protected:
    AnalysisSpectral(const unsigned int numberElements_in)
      : numberElements(numberElements_in)
    {
      reset();
    }

    inline void reset() {
      for(auto waveNumber = 0; waveNumber < maxWaveNumber; ++waveNumber){
        spectra[waveNumber] = (T) 0;
      }
    }

    inline void normalize() {
      for(auto waveNumber = 0; waveNumber < maxWaveNumber; ++waveNumber){
        spectra[waveNumber] /= gSD::sVolume();
      }
    }
  };

  template<class T, unsigned int MaxWaveNumber>
  class EnergySpectra
    : public AnalysisSpectral<T, MaxWaveNumber> {
  private:
    using Base = AnalysisSpectral<T, MaxWaveNumber>;

    T * localVelocityPtr;
    T (&spectraRef)[MaxWaveNumber];
    ForwardFFT<double, Architecture::CPU, PartitionningType::OneD,
               L::dimD, L::dimD> forward;
    BackwardFFT<double, Architecture::CPU, PartitionningType::OneD,
                L::dimD, L::dimD> backward;

  public:
    static constexpr auto analysisName = "energy_spectra";

    EnergySpectra(T * localVelocityPtr_in,
                  const unsigned int numberElements_in,
                  const ptrdiff_t globalLength_in[3])
      : Base(numberElements_in)
      , spectraRef(Base::spectra)
      , localVelocityPtr(localVelocityPtr_in)
      , forward(localVelocityPtr_in, numberElements_in, globalLength_in)
      , backward(localVelocityPtr_in, numberElements_in, globalLength_in)
    {}

    HOST
    void operator()(const Position& iFP, const unsigned int index,
                    const WaveNumber& iK, const unsigned int kNorm) {
      T coefficient = (T) 1;
      if(iK[L::dimD-1] == 0) {
        coefficient = (T) 0.5;
      }

      T energy = (T) 0;
      for(auto iD = 0; iD < L::dimD; ++iD) {
        fftw_complex * fourierVelocityPtr = ((fftw_complex *) (localVelocityPtr+iD*Base::numberElements));
        energy += fourierVelocityPtr[index][p::Re]*fourierVelocityPtr[index][p::Re];
        energy += fourierVelocityPtr[index][p::Im]*fourierVelocityPtr[index][p::Im];
      }

      if(kNorm < MaxWaveNumber) spectraRef[kNorm] += coefficient * energy;
    }

    void executeForward() {
      forward.execute();
    }

    void executeBackward() {
      backward.execute();
    }

    using Base::reset;
    using Base::normalize;
    using Base::spectra;
  };

}
