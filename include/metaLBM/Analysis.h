#ifndef ANALYSIS_H
#define ANALYSIS_H

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

  protected:
    AnalysisScalar()
      : scalar((T) 0)
    {}

    inline void reset() {
      scalar = (T) 0;
    }

    inline void normalize() {
      scalar /= lSD::sVolume();
    }
  };

  template<class T>
  class TotalEnergy
    : public AnalysisScalar<T> {
  private:
    using Base = AnalysisScalar<T>;

    T& scalarRef;
    T * * localDensityPtr;
    T * * localVelocityPtr;

  public:
    static constexpr auto analysisName = "total_energy";

    TotalEnergy(T * * localDensityPtr_in, T * * localVelocityPtr_in)
      : Base()
      , scalarRef(Base::scalar)
      , localDensityPtr(localDensityPtr_in)
      , localVelocityPtr(localVelocityPtr_in)
    {}

    HOST
    void operator()(const Position& iP) {
      auto index = lSD::getIndex(iP);
      for(auto iD = 0; iD < L::dimD; ++iD) {
        scalarRef
          += 0.5*localDensityPtr[0][index]
          *localVelocityPtr[iD][index]*localVelocityPtr[iD][index];
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

    T * * localVorticityPtr;
    T& scalarRef;

  public:
    static constexpr auto analysisName = "total_enstrophy";

    TotalEnstrophy(T * localVorticityPtr_in[L::dimD])
      : Base()
      , scalarRef(Base::scalar)
      , localVorticityPtr(localVorticityPtr_in)
    {}

    HOST
    void operator()(const Position& iP) {
      auto index = lSD::getIndex(iP);
      for(auto iD = 0; iD < 2*L::dimD-3; ++iD) {
        scalarRef
          += 0.5*localVorticityPtr[iD][index]*localVorticityPtr[iD][index];
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

  protected:
    AnalysisSpectral()
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

    T * * localVelocityPtr;
    T (&spectraRef)[MaxWaveNumber];

  public:
    static constexpr auto analysisName = "energy_spectra";

    EnergySpectra(T * * localVelocityPtr_in)
      : Base()
      , spectraRef(Base::spectra)
      , localVelocityPtr(localVelocityPtr_in)
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
        fftw_complex * fourierForcePtr = ((fftw_complex *) localVelocityPtr[iD]);

        energy += fourierForcePtr[index][p::Re]*fourierForcePtr[index][p::Re];
        energy += fourierForcePtr[index][p::Im]*fourierForcePtr[index][p::Im];
      }

      spectraRef[kNorm] += coefficient * energy;
    }

    using Base::reset;
    using Base::normalize;
    using Base::spectra;
  };



}

#endif // ANALYSIS_H
