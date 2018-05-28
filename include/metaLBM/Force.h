#pragma once

#include <cmath>

#ifdef USE_FFTW
#include "Transformer.h"
#endif

#include "Commons.h"
#include "DynamicArray.h"
#include "FourierDomain.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Options.h"

namespace lbm {

template <class T, ForceType forceType, unsigned int Dimension = 0>
class Force {};

template <class T>
class Force<T, ForceType::Generic> {
 protected:
  MathVector<T, L::dimD> amplitude;

  Force(const MathVector<T, 3>& amplitude_in) : amplitude{(T)0} {
    for (auto iD = 0; iD < L::dimD; ++iD) {
      amplitude[iD] = amplitude_in[iD];
    }
  }

  LBM_DEVICE LBM_HOST LBM_INLINE void setForce(T* localForceArray,
                                               const Position& iP,
                                               MathVector<T, L::dimD>& force,
                                               const unsigned int numberElements) {
    auto index = lSD::getIndex(iP);
#pragma unroll
    for (auto iD = 0; iD < L::dimD; ++iD) {
      force[iD] = (localForceArray + iD * numberElements)[index];
    }
  }

  LBM_DEVICE LBM_HOST inline void update(const unsigned int iteration) {}
};

template <class T>
class Force<T, ForceType::GenericTimeDependent>
    : public Force<T, ForceType::Generic> {
 private:
  using Base = Force<T, ForceType::Generic>;

 protected:
  using Base::amplitude;

  using Base::Force;

  using Base::setForce;
  using Base::update;
};

template <class T>
class Force<T, ForceType::GenericTimeIndependent>
    : public Force<T, ForceType::Generic> {
 private:
  using Base = Force<T, ForceType::Generic>;

 protected:
  using Base::amplitude;

  using Base::Force;

  using Base::setForce;
  using Base::update;
};

template <class T>
class Force<T, ForceType::None>
    : public Force<T, ForceType::GenericTimeIndependent> {
 private:
  using Base = Force<T, ForceType::GenericTimeIndependent>;

  using Base::amplitude;

 public:
  Force(const MathVector<T, 3>& amplitude_in,
        const MathVector<T, 3>& waveLength_in,
        const unsigned int kMin_in,
        const unsigned int kMax_in)
    : Base(MathVector<T, 3>{{0}}) {}

  LBM_DEVICE LBM_HOST inline void setForce(const Position& iP,
                                           MathVector<T, L::dimD>& force) {
    force = amplitude;
  }

  LBM_HOST
  inline void setLocalForceArray(double* localForcePtr,
                                 const Position& offset) {
  }

  using Base::setForce;
  using Base::update;
};


 template <class T>
class Force<T, ForceType::Constant>
    : public Force<T, ForceType::GenericTimeIndependent> {
 private:
  using Base = Force<T, ForceType::GenericTimeIndependent>;

  using Base::amplitude;

 public:
  Force(const MathVector<T, 3>& amplitude_in,
        const MathVector<T, 3>& waveLength_in,
        const unsigned int kMin_in,
        const unsigned int kMax_in)
      : Base(amplitude_in) {}

  LBM_DEVICE LBM_HOST inline void setForce(const Position& iP,
                                           MathVector<T, L::dimD>& force) {
    force = amplitude;
  }

  LBM_HOST
  inline void setLocalForceArray(double* localForcePtr,
                                 const Position& offset) {
    Computation<Architecture::CPU, L::dimD> computationLocal(lSD::sStart(),
                                                             lSD::sEnd());

    computationLocal.Do([=] LBM_HOST(const Position& iP) {
      auto index = lSD::getIndex(iP);

      MathVector<T, L::dimD> force;
      setForce(iP + offset, force);
      for (auto iD = 0; iD < L::dimD; ++iD) {
        (localForcePtr + iD * FFTWInit::numberElements)[index] = force[iD];
      }
    });
    computationLocal.synchronize();

  }

  using Base::setForce;
  using Base::update;
};

template <class T>
class Force<T, ForceType::Sinusoidal>
    : public Force<T, ForceType::GenericTimeIndependent> {
 private:
  using Base = Force<T, ForceType::GenericTimeIndependent>;

 protected:
  using Base::amplitude;
  MathVector<T, L::dimD> waveLength;

 public:
  Force(const MathVector<T, 3>& amplitude_in,
        const MathVector<T, 3>& waveLength_in,
        const unsigned int kMin_in,
        const unsigned int kMax_in)
      : Base(amplitude_in), waveLength{(T)0} {
    for (auto iD = 0; iD < L::dimD; ++iD) {
      waveLength[iD] = waveLength_in[iD];
    }
  }

  LBM_DEVICE LBM_HOST inline void setForce(const Position& iP,
                                           MathVector<T, L::dimD>& force) {
    for (auto iD = 0; iD < L::dimD; ++iD) {
      force[iD] = amplitude[iD] * sin(iP[iD] * 2 * M_PI / waveLength[iD]);
    }
  }

  LBM_HOST
  inline void setLocalForceArray(double* localForcePtr,
                                 const Position& offset) {
    Computation<Architecture::CPU, L::dimD> computationLocal(lSD::sStart(),
                                                             lSD::sEnd());
    computationLocal.Do([=] LBM_HOST(const Position& iP) {
      MathVector<T, L::dimD> force;

      auto index = lSD::getIndex(iP);
      setForce(iP + offset, force);
      for (auto iD = 0; iD < L::dimD; ++iD) {
        (localForcePtr + iD * FFTWInit::numberElements)[index] = force[iD];
      }
    });
    computationLocal.synchronize();

  }

  using Base::setForce;
  using Base::update;
};

template <class T>
class Force<T, ForceType::Kolmogorov>
    : public Force<T, ForceType::GenericTimeIndependent> {
 private:
  using Base = Force<T, ForceType::GenericTimeIndependent>;

  using Base::amplitude;
  MathVector<T, L::dimD> waveLength;

 public:
  Force(const MathVector<T, 3>& amplitude_in,
        const MathVector<T, 3>& waveLength_in,
        const unsigned int kMin_in,
        const unsigned int kMax_in)
      : Base(amplitude_in), waveLength{(T)0} {
    for (auto iD = 0; iD < L::dimD; ++iD) {
      waveLength[iD] = waveLength_in[iD];
    }
  }

  LBM_DEVICE LBM_HOST void setForce(const Position& iP,
                                    MathVector<T, L::dimD>& force) {
    force[d::X] = amplitude[d::X] * sin(iP[d::Y] * 2 * M_PI / waveLength[d::X]);
  }

  LBM_HOST
  inline void setLocalForceArray(double* localForcePtr,
                                 const Position& offset) {
    Computation<Architecture::CPU, L::dimD> computationLocal(lSD::sStart(),
                                                             lSD::sEnd());
    computationLocal.Do([=] LBM_HOST(const Position& iP) {
      MathVector<T, L::dimD> force;
      auto index = lSD::getIndex(iP);
      setForce(iP + offset, force);
      for (auto iD = 0; iD < L::dimD; ++iD) {
        (localForcePtr + iD * FFTWInit::numberElements)[index] = force[iD];
      }
    });
    computationLocal.synchronize();

  }

  using Base::setForce;
  using Base::update;
};

#ifdef USE_FFTW

template <>
class Force<double, ForceType::ConstantShell>
    : public Force<double, ForceType::GenericTimeIndependent> {
 private:
  using Base = Force<double, ForceType::GenericTimeIndependent>;

  const unsigned int kMin, kMax;

 protected:
  using Base::amplitude;

 public:
  Force(const MathVector<double, 3>& amplitude_in,
        const MathVector<double, 3>& waveLength_in,
        const unsigned int kMin_in,
        const unsigned int kMax_in)
      : Base(amplitude_in), kMin(kMin_in), kMax(kMax_in) {}

  inline void setLocalForceArray(double* localForcePtr,
                                 const Position& offset) {
    DynamicArray<double, Architecture::CPU> tempArray((2 * L::dimD - 3) *
                                                      FFTWInit::numberElements);
    double* spaceTempPtr = tempArray.data();

    Computation<Architecture::CPU, L::dimD> computationFourier(lFD::start(),
                                                               lFD::end());

    computationFourier.Do([=] LBM_HOST(const Position& iFP) {
      int kNormSquared;
      WaveNumber iK{{0}};
      Position iFP_symmetric{{0}};
      WaveNumber iK_symmetric{{0}};
      auto index = lFD::getIndex(iFP);

      iK[d::X] = iFP[d::X] + offset[d::X] <= gSD::sLength()[d::X] / 2
                     ? iFP[d::X] + offset[d::X]
                     : iFP[d::X] + offset[d::X] - gSD::sLength()[d::X];
      iK[d::Y] = iFP[d::Y] <= gSD::sLength()[d::Y] / 2
                     ? iFP[d::Y]
                     : iFP[d::Y] - gSD::sLength()[d::Y];
      iK[d::Z] = iFP[d::Z] <= gSD::sLength()[d::Z] / 2
                     ? iFP[d::Z]
                     : iFP[d::Z] - gSD::sLength()[d::Z];

      kNormSquared = iK.norm2();
      if (kNormSquared >= kMin * kMin && kNormSquared <= kMax * kMax) {
        for (auto iD = 0; iD < 2 * L::dimD - 3; ++iD) {
          fftw_complex* fourierTempPtr =
              (fftw_complex*)(spaceTempPtr + iD * FFTWInit::numberElements);

          fourierTempPtr[index][p::Re] = amplitude[iD];
          fourierTempPtr[index][p::Im] = 0;

          if (iFP[L::dimD - 1] == 0) {
            iFP_symmetric[d::X] = gSD::sLength()[d::X] - iFP[d::X];
            if (iFP_symmetric[d::X] < lFD::length()[d::X]) {
              iFP_symmetric[d::Y] = gSD::sLength()[d::Y] - iFP[d::Y];

              iK_symmetric[d::X] =
                  iFP_symmetric[d::X] + offset[d::X] <= gSD::sLength()[d::X] / 2
                      ? iFP_symmetric[d::X] + offset[d::X]
                      : iFP_symmetric[d::X] + offset[d::X] -
                            gSD::sLength()[d::X];
              iK_symmetric[d::Y] =
                  iFP_symmetric[d::Y] <= gSD::sLength()[d::Y] / 2
                      ? iFP_symmetric[d::Y]
                      : iFP_symmetric[d::Y] - gSD::sLength()[d::Y];

              auto index_symmetric = lFD::getIndex(iFP_symmetric);
              fourierTempPtr[index_symmetric][p::Re] = amplitude[iD];
              fourierTempPtr[index_symmetric][p::Im] = 0;
            }
          }
        }
      }
    });
    computationFourier.synchronize();


    MakeIncompressible<double, Architecture::CPU, PartitionningType::OneD,
                       L::dimD>
        makeIncompressible(tempArray.data(), localForcePtr, globalLengthPtrdiff_t, offset);
    makeIncompressible.executeFourier();
  }

  using Base::setForce;
  using Base::update;
};

#endif

typedef Force<dataT, forceT> Force_;

}  // namespace lbm
