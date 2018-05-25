#pragma once

#include <cmath>

#ifdef USE_FFTW
#include <fftw3-mpi.h>
#endif

#include "Commons.h"
#include "Communication.h"
#include "Computation.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Options.h"

namespace lbm {

template <class T>
class AnalysisScalar {
 public:
  T scalar;

 protected:
  AnalysisScalar()
    : scalar((T)0)
  {}

  inline void reset() { scalar = (T)0; }

  inline void normalize() { scalar /= gSD::sVolume(); }
};

template <class T>
class TotalEnergy : public AnalysisScalar<T> {
 private:
  using Base = AnalysisScalar<T>;

  T& scalarRef;
  T* localDensityPtr;
  T* localVelocityPtr;

 public:
  static constexpr auto analysisName = "total_energy";

  TotalEnergy(T* localDensityPtr_in,
              T* localVelocityPtr_in)
    : Base()
    , scalarRef(Base::scalar)
    , localDensityPtr(localDensityPtr_in)
    , localVelocityPtr(localVelocityPtr_in)
  {}

  LBM_HOST
  void operator()(const Position& iP) {
    auto index = lSD::getIndex(iP);
    for (auto iD = 0; iD < L::dimD; ++iD) {
      scalarRef += 0.5 * localDensityPtr[index] *
                   (localVelocityPtr + iD * FFTWInit::numberElements)[index] *
                   (localVelocityPtr + iD * FFTWInit::numberElements)[index];
    }
  }

  using Base::normalize;
  using Base::reset;
  using Base::scalar;
};

template <class T>
class TotalEnstrophy : public AnalysisScalar<T> {
 private:
  using Base = AnalysisScalar<T>;

  T* localVorticityPtr;
  T& scalarRef;

 public:
  static constexpr auto analysisName = "total_enstrophy";

  TotalEnstrophy(T* localVorticityPtr_in)
    : Base()
    , scalarRef(Base::scalar)
    , localVorticityPtr(localVorticityPtr_in)
  {}

  LBM_HOST
  void operator()(const Position& iP) {
    auto index = lSD::getIndex(iP);
    for (auto iD = 0; iD < 2 * L::dimD - 3; ++iD) {
      scalarRef += 0.5 *
                   (localVorticityPtr + iD * FFTWInit::numberElements)[index] *
                   (localVorticityPtr + iD * FFTWInit::numberElements)[index];
    }
  }

  using Base::normalize;
  using Base::reset;
  using Base::scalar;
};

template <class T, unsigned int maxWaveNumber>
class AnalysisSpectral {
 public:
  T spectra[maxWaveNumber];

 protected:
  AnalysisSpectral() { reset(); }

  inline void reset() {
    for (auto waveNumber = 0; waveNumber < maxWaveNumber; ++waveNumber) {
      spectra[waveNumber] = (T)0;
    }
  }

  inline void normalize() {
    for (auto waveNumber = 0; waveNumber < maxWaveNumber; ++waveNumber) {
      spectra[waveNumber] /= gSD::sVolume();
    }
  }
};

template <class T, unsigned int MaxWaveNumber>
class PowerSpectra : public AnalysisSpectral<T, MaxWaveNumber> {
 private:
  using Base = AnalysisSpectral<T, MaxWaveNumber>;

  T* localArrayPtr;
  T (&spectraRef)[MaxWaveNumber];
  ForwardFFT<double, Architecture::CPU, PartitionningType::OneD,
             L::dimD, L::dimD> forward;
  BackwardFFT<double, Architecture::CPU, PartitionningType::OneD,
              L::dimD, L::dimD> backward;

 public:
  const std::string analysisName;

  PowerSpectra(T* localArrayPtr_in,
               const ptrdiff_t globalLength_in[3],
               const std::string& analysisName_in)
    : Base()
    , spectraRef(Base::spectra)
    , localArrayPtr(localArrayPtr_in)
    , forward(localArrayPtr_in, globalLength_in)
    , backward(localArrayPtr_in, globalLength_in)
    , analysisName(analysisName_in)
  {}

  LBM_HOST
  void operator()(const Position& iFP, const unsigned int index,
                  const WaveNumber& iK, const unsigned int kNorm) {
    T coefficient = (T)1;
    if (iK[L::dimD - 1] == 0) {
      coefficient = (T)0.5;
    }

    T energy = (T)0;
    for (auto iD = 0; iD < L::dimD; ++iD) {
      fftw_complex* fourierArrayPtr =
          ((fftw_complex*)(localArrayPtr + iD * FFTWInit::numberElements));
      energy +=
          fourierArrayPtr[index][p::Re] * fourierArrayPtr[index][p::Re];
      energy +=
          fourierArrayPtr[index][p::Im] * fourierArrayPtr[index][p::Im];
    }

    if (kNorm < MaxWaveNumber)
      spectraRef[kNorm] += coefficient * energy;
  }

  void executeForward() { forward.execute(); }

  void executeBackward() { backward.execute(); }

  using Base::normalize;
  using Base::reset;
  using Base::spectra;
};

}  // namespace lbm
