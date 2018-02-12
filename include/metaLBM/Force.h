#ifndef FORCE_H
#define FORCE_H

#include <cmath>

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "StaticArray.h"
#include "Lattice.h"

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
  class Force<T, ForceType::Constant> : public Force<T, ForceType::Generic> {
  private:
    using Force<T, ForceType::Generic>::force;
    using Force<T, ForceType::Generic>::amplitude;

  public:
  Force(const MathVector<T, 3>& amplitude_in,
        const MathVector<T, 3>& waveLength_in)
      : Force<T, ForceType::Generic>(amplitude_in)
    {}

    #pragma omp declare simd
    DEVICE HOST
    inline void setForce(const MathVector<unsigned int, 3>& iP) {
      force = amplitude;
    }

    using Force<T, ForceType::Generic>::getForce;

  };

  template<class T>
  class Force<T, ForceType::Sinusoidal> : public Force<T, ForceType::Generic> {
  protected:
    using Force<T, ForceType::Generic>::force;
    using Force<T, ForceType::Generic>::amplitude;
    MathVector<T, L::dimD> waveLength;

  public:
    Force(const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in)
      : Force<T, ForceType::Generic>(amplitude_in)
      , waveLength{(T) 0}
    {
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        waveLength[iD] =  waveLength_in[iD];
      }
    }

    #pragma omp declare simd
    DEVICE HOST
    inline void setForce(const MathVector<unsigned int, 3>& iP){
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        force[iD] = amplitude[iD] * sin(iP[iD]*2*M_PI/waveLength[iD]);
      }
    }

    using Force<T, ForceType::Generic>::getForce;

  };

  template<class T>
  class Force<T, ForceType::Kolmogorov> : public Force<T, ForceType::Sinusoidal> {
  private:
    using Force<T, ForceType::Sinusoidal>::force;

    using Force<T, ForceType::Sinusoidal>::amplitude;
    using Force<T, ForceType::Sinusoidal>::waveLength;

  public:
    Force(const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in)
      : Force<T, ForceType::Sinusoidal>(amplitude_in, waveLength_in)
    {}

    #pragma omp declare simd
    DEVICE HOST
    inline void setForce(const MathVector<unsigned int, 3>& iP){
      force[d::X] = amplitude[d::X] * sin(iP[d::Y]*2*M_PI/waveLength[d::X]);
    }

    using Force<T, ForceType::Generic>::getForce;


  };


  template<class T>
  class Force<T, ForceType::Shell> : public Force<T, ForceType::Generic> {
  protected:
    using Force<T, ForceType::Generic>::force;
    using Force<T, ForceType::Generic>::amplitude;

  public:
    Force(const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in)
      : Force<T, ForceType::Generic>(amplitude_in)
      , kMin(minWavenumber)
      , kMax(maxWavenumber)
    {
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        waveLength[iD] =  waveLength_in[iD];
      }
    }

    #pragma omp declare simd
    DEVICE HOST
    inline void setForce(const MathVector<unsigned int, 3>& iP){
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
        force[iD] = amplitude[iD] * sin(iP[iD]*2*M_PI/waveLength[iD]);
      }
    }

    using Force<T, ForceType::Generic>::getForce;

  };


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
