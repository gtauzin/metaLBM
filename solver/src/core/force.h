#ifndef FORCE_H
#define FORCE_H

#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include <math.h>
#include <omp.h>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include "input.h"
#include "structure.h"
#include "commons.h"
#include "calculate.h"

namespace lbm {

  template<class T>
    class Force {

  public:
    Force()
      {}

#pragma omp declare simd
    virtual MathVector<T, L::dimD> force(const MathVector<int, 3>& iP) = 0;

  };


  template<class T>
    class ConstantForce : public Force<T> {
  private:
    MathVector<T, L::dimD> amplitude;

  public:

  ConstantForce(const MathVector<double, 3>& amplitude_in)
    : Force<T>()
    , amplitude{(T) 0}
    {
      UnrolledFor<0, L::dimD>::Do([&] (int d) {
          amplitude[d] = (T) amplitude_in[d];
      });
    }

#pragma omp declare simd
    inline MathVector<T, L::dimD> force(const MathVector<int, 3>& iP) {
      return amplitude;
    }
  };


  template<class T>
    class SinusoidalForce : public Force<T> {
  private:
    MathVector<T, L::dimD> amplitude;
    MathVector<T, L::dimD> waveLength;

  public:

  SinusoidalForce(const MathVector<double, 3>& amplitude_in,
                  const MathVector<double, 3>& waveLength_in)
    : Force<T>()
    , amplitude{(T) 0}
    , waveLength{(T) 0}
    {
      UnrolledFor<0, L::dimD>::Do([&] (int d) {
          amplitude[d] = (T) amplitude_in[d];
          waveLength[d] = (T) waveLength_in[d];
      });
    }

#pragma omp declare simd
    inline MathVector<T, L::dimD> force(const MathVector<int, 3>& iP){
      MathVector<T, L::dimD> forceR{{(T)(0)}};

      UnrolledFor<0, L::dimD>::Do([&] (int d) {
          forceR[d] = amplitude[d] * sin(iP[d]*2*M_PI/waveLength[d]);
      });

      return forceR;
    }

  };

  template<class T>
  class KolmogorovForce : public Force<T> {
  private:
    MathVector<T, L::dimD> amplitude;
    MathVector<T, L::dimD> waveLength;

  public:

  KolmogorovForce(const MathVector<double, 3>& amplitude_in,
                  const MathVector<double, 3>& waveLength_in)
    : Force<T>()
    , amplitude{(T) 0}
    , waveLength{(T) 0}
    {
      UnrolledFor<0, L::dimD>::Do([&] (int d) {
          amplitude[d] = (T) amplitude_in[d];
          waveLength[d] = (T) waveLength_in[d];
      });
    }

#pragma omp declare simd
    inline MathVector<T, L::dimD> force(const MathVector<int, 3>& iP){
      MathVector<T, L::dimD> forceR{{(T)(0)}};

      forceR[d::X] = amplitude[d::X] * sin(iP[d::Y]*2*M_PI/waveLength[d::X]);

      return forceR;
    }

  };


  template<class T>
    std::shared_ptr<Force<T>> Create(const ForceType& forceType,
                                        const MathVector<double, 3>& amplitude,
                                        const MathVector<double, 3>& waveLength) {

    switch(forceType){
    case ForceType::constant : {
      return std::shared_ptr<Force<T>>(new ConstantForce<T>(amplitude));
    }
    case ForceType::sinusoidal : {
      return std::shared_ptr<Force<T>>(new SinusoidalForce<T>(amplitude,
                                                                    waveLength));
    }
    case ForceType::kolmogorov : {
      return std::shared_ptr<Force<T>>(new KolmogorovForce<T>(amplitude,
                                                                    waveLength));
    }

    default:{
      BOOST_LOG_TRIVIAL(error) << "Unknown type of force.";
      return nullptr;
    }
    }
  }


  template<class T>
    class Forces {
  private:
    std::array<std::shared_ptr<Force<T>>, numberForces> forcesArray;

  public:
  Forces(std::array<std::shared_ptr<Force<T>>, numberForces> forcesArray_in)
    : forcesArray(forcesArray_in)
    {}


#pragma omp declare simd
    inline MathVector<T, L::dimD> force(const MathVector<int, 3>& iP) {
      MathVector<T, L::dimD> forceR = MathVector<T, L::dimD>{{(T) 0}};

      for(std::shared_ptr<Force<T>> bodyForce : forcesArray) {
        forceR += bodyForce->force(iP);
      }
      return forceR;
    }

  };

}

#endif // FORCE_H
