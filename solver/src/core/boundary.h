#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <iostream>
#include <array>
#include <vector>
#include <memory>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include "input.h"
#include "structure.h"
#include "commons.h"
#include "distribution.h"

namespace lbm {

  template<class T>
    class BoundaryCondition{
  public:
    const BoundaryPosition boundaryPosition;
    int start;
    int end;
    std::array<int, 6> present;
    std::array<int, 3> missing;

    inline int opposite(const int i) {
      if(i == 0) return 0;
      else return (i <= 4)?(i+4):((i+4)%8);
    }

  BoundaryCondition(const BoundaryPosition& boundaryPosition_in,
                    const int start_in, const int end_in)
    : boundaryPosition(boundaryPosition_in)
      , start(start_in)
      , end(end_in)
      , present()
      , missing()
      {}

    virtual void apply(Distribution<T>& lattice) = 0;

  private:
    void initBoundary(const BoundaryPosition& boundaryPosition);
  };

  template<class T>
    class BounceBack_HalfWay : public BoundaryCondition<T> {
  private:

  public:
  BounceBack_HalfWay(const BoundaryPosition& boundaryPosition_in,
                     const int start_in, const int end_in)
    : BoundaryCondition<T>(boundaryPosition_in, start_in, end_in)
      {}

    void apply(Distribution<T>& lattice)
    {}
  };

  template<class T>
    class Pressure_ZouHe : public BoundaryCondition<T> {
  private:
    T pressure;
    T density;

  public:
  Pressure_ZouHe(const BoundaryPosition& boundaryPosition_in,
                 const int start_in, const int end_in,
                 const T pressure_in)
    : BoundaryCondition<T>(boundaryPosition_in, start_in, end_in)
      , pressure(pressure_in)
      , density(pressure_in*L::inv_cs2)
      {}

    void apply(Distribution<T>& lattice)
    {}
  };

  template<class T>
    class Velocity_ZouHe : public BoundaryCondition<T> {
  private:
    const MathVector<T, L::dimD> velocity;

  public:
  Velocity_ZouHe(const BoundaryPosition& boundaryPosition_in,
                 const int start_in, const int end_in,
                 const MathVector<T, L::dimD>& velocity)
    : BoundaryCondition<T>(boundaryPosition_in, start_in, end_in)
      , velocity(velocity)
    {}

    void apply(Distribution<T>& lattice)
    {}
  };

  template<class T>
    class Entropic : public BoundaryCondition<T> {
  private:
    const T pressure;
    const T density;
    const MathVector<T, L::dimD> velocity;
    const T velocity2;

  public:
  Entropic(const BoundaryPosition& boundaryPosition_in,
           const int start_in, const int end_in,
           const T pressure_in,
           const MathVector<T, L::dimD>& velocity_in)
    : BoundaryCondition<T>(boundaryPosition_in, start_in, end_in)
      , pressure(pressure_in)
      , density(pressure_in*L::inv_cs2)
      , velocity(velocity_in)
      , velocity2(velocity_in.norm2())
      {}

    void apply(Distribution<T>& lattice)
    {}
  };

  template<class T>
    class Corner : public BoundaryCondition<T> {
  private:
    const T pressure;
    const T density;
    const MathVector<T, L::dimD> velocity;

  public:
  Corner(const BoundaryPosition& boundaryPosition_in,
         const T pressure_in,
         const MathVector<T, L::dimD>& velocity_in)
    : BoundaryCondition<T>(boundaryPosition_in, -1, -1)
      , pressure(pressure_in)
      , density(pressure_in*L::inv_cs2)
      , velocity(velocity_in)
    {}

    void apply(Distribution<T>& lattice)
    {}
  };

  template<class T>
    std::shared_ptr<BoundaryCondition<T>> Create(const BoundaryType& boundaryType,
                                                    const BoundaryPosition& boundaryPosition,
                                                    const int start, const int end,
                                                    const T pressure,
                                                    const MathVector<T, L::dimD>& velocity) {
    switch(boundaryType){
    case BoundaryType::bounceBack_halfWay:{
      return std::shared_ptr<BoundaryCondition<T>>(new BounceBack_HalfWay<T>(boundaryPosition,
                                                                                   start, end));
    }
    case BoundaryType::pressure_ZouHe:{
      return std::shared_ptr<BoundaryCondition<T>>(new Pressure_ZouHe<T>(boundaryPosition,
                                                                               start, end, pressure));
    }
    case BoundaryType::velocity_ZouHe:{
      return std::shared_ptr<BoundaryCondition<T>>(new Velocity_ZouHe<T>(boundaryPosition,
                                                                               start, end, velocity));
    }
    case BoundaryType::entropic:{
      return std::shared_ptr<BoundaryCondition<T>>(new Entropic<T>(boundaryPosition,
                                                                         start, end,
                                                                         pressure, velocity));
    }
    case BoundaryType::corner:{
      return std::shared_ptr<BoundaryCondition<T>>(new Corner<T>(boundaryPosition,
                                                                       pressure, velocity));
    }
    default:
      BOOST_LOG_TRIVIAL(error) << "Error: unknown type of boundary condition.";
      return nullptr;
    }
  }

  template<class T>
    class BoundaryConditions {
  private:
    std::vector<std::shared_ptr<BoundaryCondition<T>> > boundaryConditionsVector;

  public:
    BoundaryConditions(const std::vector<std::shared_ptr<BoundaryCondition<T>>>&
                       boundaryConditionsVector_in)
      : boundaryConditionsVector(boundaryConditionsVector_in)
    {}

    inline void apply(Distribution<T>& l){
      for(auto boundaryCondition : boundaryConditionsVector){
        boundaryCondition->apply(l);
      }
    }
  };

}

#endif // BOUNDARY_H
