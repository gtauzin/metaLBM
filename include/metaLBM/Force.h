#ifndef FORCE_H
#define FORCE_H

#include <cmath>

#ifdef USE_FFTW
  #include "FourierDomain.h"
#else
  #include "Domain.h"
#endif


#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "StaticArray.h"
#include "Lattice.h"
#include "Transformer.h"

// TODO: Functors

namespace lbm {

  template<class T, ForceType forceType, unsigned int Dimension = 0>
  class Force {};

  template<class T>
  class Force<T, ForceType::Generic> {
  protected:
    MathVector<T, L::dimD> force;

    MathVector<T, L::dimD> amplitude;

    Force(const MathVector<T, 3>& amplitude_in)
      : force{(T) 0}
      , amplitude{(T) 0}
    {
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        amplitude[iD] = amplitude_in[iD];
      }
    }

  public:
    #pragma omp declare simd
    DEVICE HOST
    inline const MathVector<T, L::dimD>& getForce(){
      return force;
    }


  };


    template<class T>
  class Force<T, ForceType::GenericTimeDependent>
    : public Force<T, ForceType::Generic> {
  protected:
    using Force<T, ForceType::Generic>::force;
    using Force<T, ForceType::Generic>::amplitude;

  public:
    using Force<T, ForceType::Generic>::Force;

    #pragma omp declare simd
    DEVICE HOST
      inline void setForce(T * localForceArray[L::dimD],
                           const MathVector<unsigned int, 3>& iP,
                           const MathVector<unsigned int, 3>& offset) {
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        force[iD] = localForceArray[iD][lSD::getIndex(iP)];
      }
    }

    #pragma omp declare simd
    DEVICE HOST
    inline void update(const unsigned int iteration) {
    }

  };


  template<class T>
  class Force<T, ForceType::GenericTimeIndependent>
    : public Force<T, ForceType::Generic> {
  protected:
    using Force<T, ForceType::Generic>::force;
    using Force<T, ForceType::Generic>::amplitude;

  public:
    using Force<T, ForceType::Generic>::Force;

    #pragma omp declare simd
    DEVICE HOST
      inline void setForce(T * localForceArray[L::dimD],
                           const MathVector<unsigned int, 3>& iP,
                           const MathVector<unsigned int, 3>& offset) {
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        force[iD] = localForceArray[iD][hSD::getIndexLocal(iP)];
      }
    }

    #pragma omp declare simd
    DEVICE HOST
    inline void update(const unsigned int iteration) {
    }

  };

  template<class T>
  class Force<T, ForceType::Constant>
    : public Force<T, ForceType::GenericTimeIndependent> {
  private:
    using Force<T, ForceType::GenericTimeIndependent>::force;
    using Force<T, ForceType::GenericTimeIndependent>::amplitude;

  public:
    Force(const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in,
          const unsigned int kMin_in, const unsigned int kMax_in)
      : Force<T, ForceType::GenericTimeIndependent>(amplitude_in)
    {}

    #pragma omp declare simd
    DEVICE HOST
    inline void setForce(MathVector<unsigned int, 3> iP,
                         const MathVector<unsigned int, 3>& offset) {
      iP += offset;
      force = amplitude;
    }

    using Force<T, ForceType::GenericTimeIndependent>::setForce;
    using Force<T, ForceType::GenericTimeIndependent>::getForce;
    using Force<T, ForceType::GenericTimeIndependent>::update;

  };


  template<class T>
  class Force<T, ForceType::Sinusoidal>
    : public Force<T, ForceType::GenericTimeIndependent> {
  protected:
    using Force<T, ForceType::GenericTimeIndependent>::force;
    using Force<T, ForceType::GenericTimeIndependent>::amplitude;
    MathVector<T, L::dimD> waveLength;

  public:
    Force(const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in,
          const unsigned int kMin_in, const unsigned int kMax_in)
      : Force<T, ForceType::GenericTimeIndependent>(amplitude_in)
      , waveLength{(T) 0}
    {
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        waveLength[iD] =  waveLength_in[iD];
      }
    }

    #pragma omp declare simd
    DEVICE HOST
    inline void setForce(MathVector<unsigned int, 3> iP,
                         const MathVector<unsigned int, 3>& offset) {
      iP += offset;
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        force[iD] = amplitude[iD] * sin(iP[iD]*2*M_PI/waveLength[iD]);
      }
    }

    using Force<T, ForceType::GenericTimeIndependent>::setForce;
    using Force<T, ForceType::GenericTimeIndependent>::getForce;
    using Force<T, ForceType::GenericTimeIndependent>::update;
  };

  template<class T>
  class Force<T, ForceType::Kolmogorov>
    : public Force<T, ForceType::Sinusoidal> {
  private:
    using Force<T, ForceType::Sinusoidal>::force;

    using Force<T, ForceType::Sinusoidal>::amplitude;
    using Force<T, ForceType::Sinusoidal>::waveLength;

  public:
    Force(const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in,
          const unsigned int kMin_in, const unsigned int kMax_in)
      : Force<T, ForceType::Sinusoidal>(amplitude_in, waveLength_in,
                                        kMin_in, kMax_in)
    {}

    #pragma omp declare simd
    DEVICE HOST
    inline void setForce(MathVector<unsigned int, 3> iP,
                         const MathVector<unsigned int, 3>& offset) {
      iP += offset;
      force[d::X] = amplitude[d::X] * sin(iP[d::Y]*2*M_PI/waveLength[d::X]);
    }

    #pragma omp declare simd
    DEVICE HOST
    inline void setForce(const MathVector<unsigned int, 3>& iP) {
      force[d::X] = amplitude[d::X] * sin(iP[d::Y]*2*M_PI/waveLength[d::X]);
    }

    using Force<T, ForceType::Sinusoidal>::setForce;
    using Force<T, ForceType::Sinusoidal>::getForce;
    using Force<T, ForceType::Sinusoidal>::update;
  };


#ifdef USE_FFTW

  template<>
  class Force<double, ForceType::ConstantShell>
    : public Force<double, ForceType::GenericTimeIndependent> {
  private:
    using Base = Force<double, ForceType::GenericTimeIndependent>;

    const unsigned int kMin, kMax;

  protected:
    using Base::force;
    using Base::amplitude;

  public:
    Force(const MathVector<double, 3>& amplitude_in,
          const MathVector<double, 3>& waveLength_in,
          const unsigned int kMin_in, const unsigned int kMax_in)
      : Base(amplitude_in)
      , kMin(kMin_in)
      , kMax(kMax_in)
    {}

    #pragma omp declare simd
    DEVICE HOST
    inline void setForce(MathVector<unsigned int, 3> iP,
                         const MathVector<unsigned int, 3>& offset){
      iP += offset;
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        //force[iD] = amplitude[iD] * sin(iP[iD]*2*M_PI/waveLength[iD]);
      }
    }

    inline void setLocalForceArray(double * localForcePtr[L::dimD],
                                   const unsigned int offsetX) {


      Computation<Architecture::CPU, L::dimD> computation(lFD::start(), lFD::end());

      computation.Do([=] HOST (const Position& iFP) {
          int kNormSquared;
          WaveNumber iK{{0}};
          Position iFP_symmetric{{0}};
          WaveNumber iK_symmetric{{0}};

          iK[d::X] = iFP[d::X]+offsetX <= gSD::sLength()[d::X]/2 ?
            iFP[d::X]+offsetX : iFP[d::X]+offsetX-gSD::sLength()[d::X];
          iK[d::Y] = iFP[d::Y] <= gSD::sLength()[d::Y]/2 ?
            iFP[d::Y] : iFP[d::Y]-gSD::sLength()[d::Y];
          iK[d::Z] = iFP[d::Z] <= gSD::sLength()[d::Z]/2 ?
            iFP[d::Z] : iFP[d::Z]-gSD::sLength()[d::Z];

          kNormSquared = iK.norm2();
          if (kNormSquared >= kMin*kMin && kNormSquared <= kMax*kMax) {
            auto index = lFD::getIndex(iFP);

            for(auto iD = 0; iD < L::dimD; ++iD) {
              fftw_complex * fourierForcePtr = ((fftw_complex *) localForcePtr[iD]);

              fourierForcePtr[index][0] = 1;
              fourierForcePtr[index][1] = 0;

              if(iFP[L::dimD-1] == 0) {
                iFP_symmetric[d::X] = gSD::sLength()[d::X] - iFP[d::X];
                if(iFP_symmetric[d::X] < lFD::length()[d::X]) {
                  iFP_symmetric[d::Y] = gSD::sLength()[d::Y] - iFP[d::Y];

                  iK_symmetric[d::X]
                    = iFP_symmetric[d::X]+offsetX <= gSD::sLength()[d::X]/2 ?
                    iFP_symmetric[d::X]+offsetX :
                    iFP_symmetric[d::X]+offsetX - gSD::sLength()[d::X];
                  iK_symmetric[d::Y]
                    = iFP_symmetric[d::Y] <= gSD::sLength()[d::Y]/2 ?
                    iFP_symmetric[d::Y] : iFP_symmetric[d::Y]-gSD::sLength()[d::Y];

                  auto index_symmetric = lFD::getIndex(iFP_symmetric);
                  fourierForcePtr[index_symmetric][0] = 1;
                  fourierForcePtr[index_symmetric][1] = 0;
                }
              }
            }
          }
        });

      FourierTransformer<double, Architecture::CPU, PartitionningType::OneD,
                         L::dimD, L::dimD> fftTransformer(localForcePtr,
                                                          Cast<unsigned int,
                                                          ptrdiff_t, 3>::Do(gSD::sLength()).data());
    fftTransformer.executeBackwardFFT();

    }

    using Base::setForce;
    using Base::getForce;
    using Base::update;
  };



#endif

  typedef Force<dataT, forceT> Force_;

}

#endif // FORCE_H
