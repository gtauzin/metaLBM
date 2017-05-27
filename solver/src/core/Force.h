#ifndef FORCE_H
#define FORCE_H

#include <cmath>
#include <omp.h>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include "Options.h"
#include "MathVector.h"
#include "StaticArray.h"

// TODO: Functors

namespace lbm {

  template<class T, ForceType forceType>
  class Force {};

  template<class T>
  class Force<T, ForceType::Generic> {
  protected:
    MathVector<unsigned int, 3> offset;
    MathVector<T, L::dimD> amplitude;

    Force(const MathVector<T, 3>& amplitude_in)
    {
      UnrolledFor<0, L::dimD>::Do([&] (int iD) {
          amplitude[iD] = amplitude_in[iD];
      });
    }

  public:
    #pragma omp declare simd
    inline MathVector<T, L::dimD> force(const MathVector<int, 3>& iP) {
      // ABORT EXCEPTIONS
      return MathVector<T, L::dimD>{{(T)(0)}};
    }
  };


  template<class T>
  class Force<T, ForceType::Constant> : public Force<T, ForceType::Generic> {
  private:
    using Force<T, ForceType::Generic>::amplitude;

  public:
  Force(const MathVector<unsigned int, 3> offset_in = MathVector<unsigned int, 3>{{0}},
        const MathVector<T, 3>& amplitude_in)
      : Force<T, ForceType::Generic>(amplitude_in)
    {}

    #pragma omp declare simd
    inline MathVector<T, L::dimD> force(const MathVector<int, 3>& iP) {
      return amplitude;
    }
  };

  template<class T>
  class Force<T, ForceType::Sinusoidal> : public Force<T, ForceType::Generic> {
  private:
    using Force<T, ForceType::Generic>::amplitude;
    MathVector<T, L::dimD> waveLength;

  public:
    Force(const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in)
      : Force<T, ForceType::Generic>(amplitude_in)
      , waveLength{(T) 0}
    {
      UnrolledFor<0, L::dimD>::Do([&] (int iD) {
          waveLength[iD] =  waveLength_in[iD];
      });
    }

    #pragma omp declare simd
    inline MathVector<T, L::dimD> force(const MathVector<int, 3>& iP){
      MathVector<T, L::dimD> forceR{{ (T) 0}};

      UnrolledFor<0, L::dimD>::Do([&] (int iD) {
          forceR[iD] = amplitude[iD] * sin(iP[iD]*2*M_PI/waveLength[iD]);
      });

      return forceR;
    }
  };

  template<class T>
  class Force<T, ForceType::Kolmogorov> : public Force<T, ForceType::Sinusoidal> {
  private:
    using Force<T, ForceType::Sinusoidal>::amplitude;
    using Force<T, ForceType::Sinusoidal>::waveLength;

  public:

    Force(const MathVector<T, 3>& amplitude_in,
          const MathVector<T, 3>& waveLength_in)
      : Force<T, ForceType::Sinusoidal>(amplitude_in, waveLength_in)
    {}

    #pragma omp declare simd
    inline MathVector<T, L::dimD> force(const MathVector<int, 3>& iP){
      MathVector<T, L::dimD> forceR{{ (T) 0}};

      forceR[d::X] = amplitude[d::X] * sin(iP[d::Y]*2*M_PI/waveLength[d::X]);

      return forceR;
    }

  };

  template<class T>
  class Forces {
  private:
    MathVector<T, L::dimD> force = MathVector<T, L::dimD>{{(T) 0}};
    const StaticArray<Force<T, ForceType::Generic>, numberForces> forcesArray;
    const MathVector<unsigned int, 3> offset;

  public:
  Forces(const StaticArray<Force<T, ForceType::Generic>, numberForces>& forcesArray_in,
         const MathVector<unsigned int, 3>& offset_in = MathVector<unsigned int, 3>{{0}})
      : forcesArray(forcesArray_in)
      , offset(offset_in)
    {}


    #pragma omp declare simd
    inline MathVector<T, L::dimD> force(const MathVector<int, 3>& iP) {
      BOOST_LOG_TRIVIAL(debug) << " - Computing force.";
       force = MathVector<T, L::dimD>{{(T) 0}};

      for(Force<T, ForceType::Generic> force : forcesArray) {
        force += force.force(offset + iP);
      }

      return force;
    }

    inline MathVector<T, L::dimD> getForce() {
      return force;
    }

  };

}

#endif // FORCE_H
