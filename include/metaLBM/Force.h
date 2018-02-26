#ifndef FORCE_H
#define FORCE_H

#include <cmath>

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "StaticArray.h"
#include "Lattice.h"
#include "Analyzer.h"

// TODO: Functors

namespace lbm {

  template<class T, ForceType forceType>
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
    inline void setForce(const T * RESTRICT localForceArray,
                         const MathVector<unsigned int, 3>& iP,
                         const MathVector<unsigned int, 3>& offset) {
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        force[iD] = localForceArray[lSD::getIndex(iP, iD)];
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
    inline void setForce(const T * RESTRICT localForceArray,
                         const MathVector<unsigned int, 3>& iP,
                         const MathVector<unsigned int, 3>& offset) {
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        force[iD] = localForceArray[lSD::getIndex(iP, iD)];
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
    inline void setForce(const MathVector<unsigned int, 3>& iP,
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
    inline void setForce(const MathVector<unsigned int, 3>& iP,
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
    inline void setForce(const MathVector<unsigned int, 3>& iP,
                         const MathVector<unsigned int, 3>& offset) {
      //iP += offset;
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
    const unsigned int kMin, kMax;

  protected:
    using Force<double, ForceType::GenericTimeIndependent>::force;
    using Force<double, ForceType::GenericTimeIndependent>::amplitude;

  public:
    Force(const MathVector<double, 3>& amplitude_in,
          const MathVector<double, 3>& waveLength_in,
          const unsigned int kMin_in, const unsigned int kMax_in)
      : Force<double, ForceType::GenericTimeIndependent>(amplitude_in)
      , kMin(kMin_in)
      , kMax(kMax_in)
    {}

    #pragma omp declare simd
    DEVICE HOST
    inline void setForce(const MathVector<unsigned int, 3>& iP,
                         const MathVector<unsigned int, 3>& offset){
      //iP += offset;
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        //force[iD] = amplitude[iD] * sin(iP[iD]*2*M_PI/waveLength[iD]);
      }
    }

    inline void setGlobalForceArray(double * RESTRICT globalForceArray) {
      DynamicArray<fftw_complex, Architecture::CPU>
        globalForceFourierArray(L::dimD * gFD::volume());

      int kNormSquared;
      MathVector<int, 3> iK{{0}};
      MathVector<unsigned int, 3> iFP{{0}};
      for(unsigned int iFX = gFD::start()[d::X]; iFX < gFD::end()[d::X]; ++iFX) {
        iFP[d::X] = iFX;
        iK[d::X] = gFD::start()[d::X]/2 ? iFX : iFX - gFD::start()[d::X];
        for(unsigned int iFY = gFD::start()[d::Y]; iFY < gFD::end()[d::Y]; ++iFY) {
          iFP[d::Y] = iFY;
          iK[d::Y] = gFD::start()[d::Y]/2 ? iFY : iFY - gFD::start()[d::Y];
          for(unsigned int iFZ = gFD::start()[d::Z]; iFZ < gFD::end()[d::Z]; ++iFZ) {
            iFP[d::Z] = iFZ;
            iK[d::Z] = gFD::start()[d::Z]/2 ? iFZ : iFZ - gFD::start()[d::Z];

            kNormSquared = iK.norm2();
            if (kNormSquared >= kMin*kMin && kNormSquared <= kMax*kMax) {
              for(unsigned int iD = 0; iD < L::dimD; ++iD) {
                globalForceFourierArray[gFD::getIndex(iFP, iD)] = {amplitude[iD], amplitude[iD]};
              }
            }

          }
        }
      }

      for(unsigned int iD = 0; iD < L::dimD; ++iD) {

      }

    }

    using Force<double, ForceType::GenericTimeIndependent>::setForce;
    using Force<double, ForceType::GenericTimeIndependent>::getForce;
    using Force<double, ForceType::GenericTimeIndependent>::update;
  };
#endif

  /* template<class T, unsigned int NumberForces> */
  /* class Forces { */
  /* private: */
  /*   MathVector<T, L::dimD> force; */
  /*   const StaticArray<Force<T, ForceType::Generic>, NumberForces> forceArray; */

  /* public: */
  /* Forces(const StaticArray<Force<T, ForceType::Generic>, NumberForces>& forceArray_in) */
  /*   : force{(T) 0} */
  /*   , forceArray(forceArray_in) */
  /*   {} */

  /*   #pragma omp declare simd */
  /*   inline void update() { */

  /*   } */

  /*   #pragma omp declare simd */
  /*   inline MathVector<T, L::dimD> getForce(const MathVector<int, 3>& iP) { */
  /*      force = MathVector<T, L::dimD>{{(T) 0}}; */

  /*     for(Force<T, ForceType::Generic> force : forceArray) { */
  /*       force += force.force(iP); */
  /*     } */

  /*     return force; */
  /*   } */

  /*   #pragma omp declare simd */
  /*   inline MathVector<T, L::dimD>& getForce() { */
  /*     return force; */
  /*   } */

  /* }; */

  typedef Force<dataT, forceT> Force_;

}

#endif // FORCE_H
