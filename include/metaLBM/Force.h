#pragma once

#include <cmath>

#ifdef USE_FFTW
#include "Transformer.h"
#endif

#include "Commons.h"
#include "DynamicArray.h"
#include "FourierDomain.h"
#include "FieldList.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Options.h"

namespace lbm {

  template <class T, ForceType forceType, Architecture architecture,
            unsigned int Dimension = 0>
  class Force {};

  template <class T, Architecture architecture>
  class Force<T, ForceType::Generic, architecture> {
  protected:
    Position offset;
    MathVector<T, L::dimD> amplitude;

    Force(const Position& offset_in,
          const MathVector<T, 3>& amplitude_in)
      : offset(offset_in)
      , amplitude{(T)0}
    {
      for (auto iD = 0; iD < L::dimD; ++iD) {
        amplitude[iD] = amplitude_in[iD];
      }
    }

    LBM_DEVICE LBM_HOST LBM_INLINE
    void setForce(T* forcePtr, const Position& iP,
                  MathVector<T, L::dimD>& force,
                  const unsigned int numberElements) {
      auto index = lSD::getIndex(iP);
      #pragma unroll
      for (auto iD = 0; iD < L::dimD; ++iD) {
        force[iD] = (forcePtr + iD * numberElements)[index];
      }
    }


    LBM_HOST inline
    void update(double * forcePtr,
                FieldList<double, architecture>& fieldList,
                const unsigned int iteration) {}
  };

  template <class T, Architecture architecture>
  class Force<T, ForceType::GenericTimeDependent, architecture>
    : public Force<T, ForceType::Generic, architecture> {
  private:
    using Base = Force<T, ForceType::Generic, architecture>;

  protected:
    using Base::offset;
    using Base::amplitude;

    using Base::Force;

    LBM_HOST
    inline void reset(double * forcePtr) {
      Computation<Architecture::CPU, L::dimD> computationLocal(lSD::sStart(),
                                                               lSD::sEnd());
      computationLocal.Do([=] LBM_HOST(const Position& iP) {
          auto index = lSD::getIndex(iP);

          for (auto iD = 0; iD < L::dimD; ++iD) {
            (forcePtr + iD * FFTWInit::numberElements)[index] = (T) 0;
          }
        });
      computationLocal.synchronize();

    }

    using Base::setForce;
    using Base::update;
  };

  template <class T, Architecture architecture>
  class Force<T, ForceType::GenericTimeIndependent, architecture>
    : public Force<T, ForceType::Generic, architecture> {
  private:
    using Base = Force<T, ForceType::Generic, architecture>;

  protected:
    using Base::offset;
    using Base::amplitude;

    using Base::Force;

    using Base::setForce;
    using Base::update;
  };

  template <class T, Architecture architecture>
  class Force<T, ForceType::None, architecture>
    : public Force<T, ForceType::GenericTimeIndependent, architecture> {
  private:
    using Base = Force<T, ForceType::GenericTimeIndependent, architecture>;

    using Base::amplitude;

  public:
    Force(const Position& offset_in,
          const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in,
          const unsigned int kMin_in,
          const unsigned int kMax_in)
      : Base(offset_in, MathVector<T, 3>{{0}}) {}

    LBM_DEVICE LBM_HOST inline void setForce(T* forceArray,
                                             const Position& iP,
                                             MathVector<T, L::dimD>& force,
                                             const unsigned int numberElements) {
      force = amplitude;
    }

    LBM_HOST
    inline void setForceArray(double * forcePtr,
                              FieldList<double, architecture>& fieldList) {
    }

    using Base::setForce;
    using Base::update;
  };


  template <class T, Architecture architecture>
  class Force<T, ForceType::Constant, architecture>
    : public Force<T, ForceType::GenericTimeIndependent, architecture> {
  private:
    using Base = Force<T, ForceType::GenericTimeIndependent, architecture>;

    using Base::offset;
    using Base::amplitude;

  public:
    Force(const Position& offset_in,
          const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in,
          const unsigned int kMin_in,
          const unsigned int kMax_in)
      : Base(offset_in, amplitude_in) {}

    LBM_DEVICE LBM_HOST inline void setForce(T* forceArray,
                                             const Position& iP,
                                             MathVector<T, L::dimD>& force,
                                             const unsigned int numberElements) {
      force = amplitude;
    }

    LBM_HOST
    inline void setForceArray(double * forcePtr,
                              FieldList<double, architecture>& fieldList) {
      Computation<Architecture::CPU, L::dimD> computationLocal(lSD::sStart(),
                                                               lSD::sEnd());

      computationLocal.Do([=] LBM_HOST(const Position& iP) {
          auto index = lSD::getIndex(iP);

          MathVector<T, L::dimD> force;
          setForce(iP + offset, force);
          for (auto iD = 0; iD < L::dimD; ++iD) {
            (forcePtr + iD * FFTWInit::numberElements)[index] = force[iD];
          }
        });
      computationLocal.synchronize();

    }

    using Base::setForce;
    using Base::update;
  };

  template <class T, Architecture architecture>
  class Force<T, ForceType::Sinusoidal, architecture>
    : public Force<T, ForceType::GenericTimeIndependent, architecture> {
  private:
    using Base = Force<T, ForceType::GenericTimeIndependent, architecture>;

  protected:
    using Base::offset;
    using Base::amplitude;
    MathVector<T, L::dimD> waveLength;

  public:
    Force(const Position& offset_in,
          const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in,
          const unsigned int kMin_in,
          const unsigned int kMax_in)
      : Base(offset_in, amplitude_in), waveLength{(T)0}
    {
      for (auto iD = 0; iD < L::dimD; ++iD) {
        waveLength[iD] = waveLength_in[iD];
      }
    }

    LBM_DEVICE LBM_HOST inline void setForce(T* forceArray,
                                             const Position& iP,
                                             MathVector<T, L::dimD>& force,
                                             const unsigned int numberElements) {
      for (auto iD = 0; iD < L::dimD; ++iD) {
        force[iD] = amplitude[iD] * sin(iP[iD] * 2 * M_PI / waveLength[iD]);
      }
    }

    LBM_HOST
    inline void setForceArray(T * forcePtr,
                              FieldList<T, architecture>& fieldList) {
      Computation<Architecture::CPU, L::dimD> computationLocal(lSD::sStart(),
                                                               lSD::sEnd());
      computationLocal.Do([=] LBM_HOST(const Position& iP) {
          MathVector<T, L::dimD> force;

          auto index = lSD::getIndex(iP);
          setForce(iP + offset, force);
          for (auto iD = 0; iD < L::dimD; ++iD) {
            (forcePtr + iD * FFTWInit::numberElements)[index] = force[iD];
          }
        });
      computationLocal.synchronize();

    }

    using Base::setForce;
    using Base::update;
  };

  template <class T, Architecture architecture>
  class Force<T, ForceType::Kolmogorov, architecture>
    : public Force<T, ForceType::GenericTimeIndependent, architecture> {
  private:
    using Base = Force<T, ForceType::GenericTimeIndependent, architecture>;

    using Base::offset;
    using Base::amplitude;
    MathVector<T, L::dimD> waveLength;

  public:
    Force(const Position& offset_in,
          const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in,
          const unsigned int kMin_in,
          const unsigned int kMax_in)
      : Base(offset_in, amplitude_in), waveLength{(T)0}
    {
      for (auto iD = 0; iD < L::dimD; ++iD) {
        waveLength[iD] = waveLength_in[iD];
      }
    }

    LBM_DEVICE LBM_HOST void setForce(T* forcePtr,
                                      const Position& iP,
                                      MathVector<T, L::dimD>& force,
                                      const unsigned int numberElements) {
      force[d::X] = amplitude[d::X] * sin(iP[d::Y] * 2 * M_PI / waveLength[d::X]);
    }

    LBM_HOST
    inline void setForceArray(double * forcePtr,
                              FieldList<double, architecture>& fieldList) {
      Computation<Architecture::CPU, L::dimD> computationLocal(lSD::sStart(),
                                                               lSD::sEnd());

      unsigned int numberElements = FFTWInit::numberElements;

      computationLocal.Do([=] LBM_HOST(const Position& iP) {
          MathVector<T, L::dimD> force;
          auto index = lSD::getIndex(iP);
          setForce(forcePtr, iP + offset, force, numberElements);

          for (auto iD = 0; iD < L::dimD; ++iD) {
            (forcePtr + iD * numberElements)[index] = force[iD];
          }
        });
      computationLocal.synchronize();

    }

    using Base::setForce;
    using Base::update;
  };

#ifdef USE_FFTW

  template <Architecture architecture>
  class Force<double, ForceType::ConstantShell, architecture>
    : public Force<double, ForceType::GenericTimeIndependent, architecture> {
  private:
    using Base = Force<double, ForceType::GenericTimeIndependent, architecture>;

    const unsigned int kMin, kMax;
    DynamicArray<double, Architecture::CPU> tempArray;


  protected:
    using Base::offset;
    using Base::amplitude;

  public:
    Force(const Position& offset_in,
          const MathVector<double, 3>& amplitude_in,
          const MathVector<double, 3>& waveLength_in,
          const unsigned int kMin_in,
          const unsigned int kMax_in)
      : Base(offset_in, amplitude_in), kMin(kMin_in), kMax(kMax_in)
      , tempArray((2 * L::dimD - 3) * FFTWInit::numberElements)
    {
      initTempArray();
    }

    inline void setForceArray(double * forcePtr,
                              FieldList<double, architecture>& fieldList) {
      double* spaceTempPtr = tempArray.data();

      MakeIncompressible<double, Architecture::CPU, PartitionningType::OneD, L::dimD>
        makeIncompressible(spaceTempPtr, forcePtr, globalLengthPtrdiff_t, offset);

      makeIncompressible.executeFourier();
    }

    inline void initTempArray() {
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
            ? iFP[d::Y] : iFP[d::Y] - gSD::sLength()[d::Y];
          iK[d::Z] = iFP[d::Z] <= gSD::sLength()[d::Z] / 2
            ? iFP[d::Z] : iFP[d::Z] - gSD::sLength()[d::Z];

          for (auto iD = 0; iD < 2 * L::dimD - 3; ++iD) {
            fftw_complex* fourierTempPtr =
              (fftw_complex*)(spaceTempPtr + iD * FFTWInit::numberElements);

            fourierTempPtr[index][p::Re] = 0.0;
            fourierTempPtr[index][p::Im] = 0.0;

            if (iFP[L::dimD - 1] == 0) {
              iFP_symmetric[d::X] = gSD::sLength()[d::X] - iFP[d::X];
              if (iFP_symmetric[d::X] < lFD::length()[d::X]) {
                iFP_symmetric[d::Y] = gSD::sLength()[d::Y] - iFP[d::Y];
                iFP_symmetric[L::dimD - 1] = 0;

                iK_symmetric[d::X] =
                  iFP_symmetric[d::X] + offset[d::X] <= gSD::sLength()[d::X] / 2
                  ? iFP_symmetric[d::X] + offset[d::X]
                  : iFP_symmetric[d::X] + offset[d::X] - gSD::sLength()[d::X];
                iK_symmetric[d::Y] = iFP_symmetric[d::Y] <= gSD::sLength()[d::Y] / 2
                  ? iFP_symmetric[d::Y] : iFP_symmetric[d::Y] - gSD::sLength()[d::Y];

                auto index_symmetric = lFD::getIndex(iFP_symmetric);
                fourierTempPtr[index_symmetric][p::Re] = 0.0;
                fourierTempPtr[index_symmetric][p::Im] = 0.0;
              }
            }
          }

          kNormSquared = iK.norm2();
          if (kNormSquared >= kMin * kMin && kNormSquared <= kMax * kMax) {
            for (auto iD = 0; iD < 2 * L::dimD - 3; ++iD) {
              fftw_complex* fourierTempPtr =
                (fftw_complex*)(spaceTempPtr + iD * FFTWInit::numberElements);

              fourierTempPtr[index][p::Re] = amplitude[iD];
              fourierTempPtr[index][p::Im] = 0.0;

              if (iFP[L::dimD - 1] == 0) {
                iFP_symmetric[d::X] = gSD::sLength()[d::X] - iFP[d::X];
                if (iFP_symmetric[d::X] < lFD::length()[d::X]) {
                  iFP_symmetric[d::Y] = gSD::sLength()[d::Y] - iFP[d::Y];
                  iFP_symmetric[L::dimD - 1] = 0;

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
                  fourierTempPtr[index_symmetric][p::Im] = 0.0;
                }
              }
            }
          }
      });

      computationFourier.synchronize();
    }

    using Base::setForce;
    using Base::update;
  };


  template <Architecture architecture>
  class Force<double, ForceType::EnergyRemoval, architecture>
    : public Force<double, ForceType::GenericTimeDependent, architecture> {
  private:
    using Base = Force<double, ForceType::GenericTimeDependent, architecture>;

    const unsigned int kMin, kMax;

  protected:
    using Base::offset;
    using Base::amplitude;

  public:
    Force(const Position& offset_in,
          const MathVector<double, 3>& amplitude_in,
          const MathVector<double, 3>& waveLength_in,
          const unsigned int kMin_in,
          const unsigned int kMax_in)
      : Base(offset_in, amplitude_in), kMin(kMin_in), kMax(kMax_in) {}

    inline void setForceArray(double * forcePtr,
                              FieldList<double, architecture>& fieldList) {
      unsigned int numberElements = FFTWInit::numberElements;

      DynamicArray<double, Architecture::CPU> momentumArray(L::dimD * numberElements);
      double* momentumPtr = momentumArray.data();

      const double * velocityPtr = fieldList.velocity.getData(numberElements);
      const double * densityPtr = fieldList.density.getData(numberElements);

      Computation<Architecture::CPU, L::dimD> computationLocal(lSD::sStart(),
                                                               lSD::sEnd());

      computationLocal.Do([=] LBM_HOST(const Position& iP) {
          const auto index = lSD::getIndex(iP);

          momentumPtr[index] = 0;
          for (auto iD = 0; iD < L::dimD; ++iD) {
            momentumPtr[index] +=
              (velocityPtr + iD * numberElements)[index];
          }
          momentumPtr[index] *= densityPtr[index];
        });
      computationLocal.synchronize();

      ForwardFFT<double, Architecture::CPU, PartitionningType::OneD,
                 L::dimD, L::dimD> forward(momentumPtr, globalLengthPtrdiff_t);
      forward.execute();

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
            ? iFP[d::Y] : iFP[d::Y] - gSD::sLength()[d::Y];
          iK[d::Z] = iFP[d::Z] <= gSD::sLength()[d::Z] / 2
            ? iFP[d::Z] : iFP[d::Z] - gSD::sLength()[d::Z];

          for (auto iD = 0; iD < L::dimD; ++iD) {
            fftw_complex* fourierForcePtr_iD =
              (fftw_complex*)(forcePtr + iD * numberElements);

            fourierForcePtr_iD[index][p::Re] = 0.0;
            fourierForcePtr_iD[index][p::Im] = 0.0;

            if (iFP[L::dimD - 1] == 0) {
              iFP_symmetric[d::X] = gSD::sLength()[d::X] - iFP[d::X];
              if (iFP_symmetric[d::X] < lFD::length()[d::X]) {
                iFP_symmetric[d::Y] = gSD::sLength()[d::Y] - iFP[d::Y];
                iFP_symmetric[L::dimD - 1] = 0;

                iK_symmetric[d::X] =
                  iFP_symmetric[d::X] + offset[d::X] <= gSD::sLength()[d::X] / 2
                  ? iFP_symmetric[d::X] + offset[d::X]
                  : iFP_symmetric[d::X] + offset[d::X] - gSD::sLength()[d::X];
                iK_symmetric[d::Y] = iFP_symmetric[d::Y] <= gSD::sLength()[d::Y] / 2
                  ? iFP_symmetric[d::Y] : iFP_symmetric[d::Y] - gSD::sLength()[d::Y];

                auto index_symmetric = lFD::getIndex(iFP_symmetric);
                fourierForcePtr_iD[index_symmetric][p::Re] = 0.0;
                fourierForcePtr_iD[index_symmetric][p::Im] = 0.0;
              }
            }
          }

          kNormSquared = iK.norm2();
          if (kNormSquared >= kMin * kMin && kNormSquared <= kMax * kMax) {
            for (auto iD = 0; iD < L::dimD; ++iD) {
              fftw_complex* fourierForcePtr_iD =
                (fftw_complex*)(forcePtr + iD * numberElements);
              fftw_complex* fourierMomentumPtr_iD =
                (fftw_complex*)(momentumPtr + iD * numberElements);

              fourierForcePtr_iD[index][p::Re] =
                - amplitude[iD] * fourierMomentumPtr_iD[index][p::Re];
              fourierForcePtr_iD[index][p::Im] =
                - amplitude[iD] * fourierMomentumPtr_iD[index][p::Im];

              if (iFP[L::dimD - 1] == 0) {
                iFP_symmetric[d::X] = gSD::sLength()[d::X] - iFP[d::X];
                if (iFP_symmetric[d::X] < lFD::length()[d::X]) {
                  iFP_symmetric[d::Y] = gSD::sLength()[d::Y] - iFP[d::Y];
                  iFP_symmetric[L::dimD - 1] = 0;

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
                  fourierForcePtr_iD[index_symmetric][p::Re] =
                    - amplitude[iD] * fourierMomentumPtr_iD[index_symmetric][p::Re];
                  fourierForcePtr_iD[index_symmetric][p::Im] =
                    - amplitude[iD] * fourierMomentumPtr_iD[index_symmetric][p::Im];
                }
              }
            }
          }
        });
      computationFourier.synchronize();

      BackwardFFT<double, Architecture::CPU, PartitionningType::OneD,
                  L::dimD, L::dimD> backward(forcePtr, globalLengthPtrdiff_t);
      backward.execute();
    }

    LBM_HOST inline
    void update(double * forcePtr,
                FieldList<double, architecture>& fieldList,
                const unsigned int iteration) {
      setForceArray(forcePtr, fieldList);
    }

    using Base::setForce;
  };


  template <Architecture architecture>
  class Force<double, ForceType::Turbulent2D, architecture>
    : public Force<double, ForceType::GenericTimeDependent, architecture> {
  private:
    using Base = Force<double, ForceType::GenericTimeDependent, architecture>;

    using Base::amplitude;
    MathVector<double, L::dimD> waveLength;
    Force<double, ForceType::ConstantShell, architecture> injection;
    Force<double, ForceType::EnergyRemoval, architecture> removal;

  public:
    Force(const Position& offset_in, const MathVector<double, 3>& amplitude_in,
          const MathVector<double, 3>& waveLength_in,
          const unsigned int kMin_in, const unsigned int kMax_in)
      : Base(offset_in, MathVector<double, 3>{0})
      , injection(offset_in, forceAmplitude, forceWaveLength, forcekMin, forcekMax)
      , removal(offset_in, removalForceAmplitude, removalForceWaveLength,
                removalForcekMin, removalForcekMax)
    {}

    LBM_HOST inline
    void setForceArray(double * forcePtr,
                       FieldList<double, architecture>& fieldList) {
      Computation<Architecture::CPU, L::dimD> computationLocal(lSD::sStart(),
                                                               lSD::sEnd());
      unsigned int numberElements = FFTWInit::numberElements;

      injection.setForceArray(forcePtr, fieldList);

      DynamicArray<double, Architecture::CPU> removalForceArray(L::dimD * numberElements);
      double * removalForcePtr = removalForceArray.data();

      removal.setForceArray(removalForcePtr, fieldList);

      computationLocal.Do([=] LBM_HOST(const Position& iP) {
        auto index = lSD::getIndex(iP);

        for(auto iD = 0; iD < L::dimD; ++iD) {
          (forcePtr + iD * numberElements)[index]
            += (removalForcePtr + iD * numberElements)[index];
        }
      });
      computationLocal.synchronize();
    }

    LBM_HOST inline
    void update(double * forcePtr,
                FieldList<double, architecture>& fieldList,
                const unsigned int iteration) {
      setForceArray(forcePtr, fieldList);
    }

    using Base::setForce;
  };


#endif

  template <Architecture architecture>
  using Force_ = Force<dataT, forceT, architecture>;

}  // namespace lbm
