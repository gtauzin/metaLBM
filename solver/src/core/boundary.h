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
#include "lattice.h"

namespace lbm {

  template<class T, LatticeType L>
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

    virtual void apply(Lattice<T, L>& lattice) = 0;

  private:
    void initBoundary(const BoundaryPosition& boundaryPosition);
  };

  template<class T, LatticeType L>
    class BounceBack_HalfWay : public BoundaryCondition<T, L> {
  private:

  public:
  BounceBack_HalfWay(const BoundaryPosition& boundaryPosition_in,
                     const int start_in, const int end_in)
    : BoundaryCondition<T, L>(boundaryPosition_in, start_in, end_in)
      {}

    void apply(Lattice<T, L>& lattice)
    {}
  };

  template<class T, LatticeType L>
    class Pressure_ZouHe : public BoundaryCondition<T, L> {
  private:
    T pressure;
    T density;

  public:
  Pressure_ZouHe(const BoundaryPosition& boundaryPosition_in,
                 const int start_in, const int end_in,
                 const T pressure_in)
    : BoundaryCondition<T, L>(boundaryPosition_in, start_in, end_in)
      , pressure(pressure_in)
      , density(pressure_in*inv_cs2<T, L>())
      {}

    void apply(Lattice<T, L>& lattice)
    {}
  };

  template<class T, LatticeType L>
    class Velocity_ZouHe : public BoundaryCondition<T, L> {
  private:
    const MathVector<T, dimD<T, L>()> velocity;

  public:
  Velocity_ZouHe(const BoundaryPosition& boundaryPosition_in,
                 const int start_in, const int end_in,
                 const MathVector<T, dimD<T, L>()>& velocity)
    : BoundaryCondition<T, L>(boundaryPosition_in, start_in, end_in)
      , velocity(velocity)
    {}

    void apply(Lattice<T, L>& lattice)
    {}
  };

  template<class T, LatticeType L>
  class Entropic : public BoundaryCondition<T, L> {
  private:
    const T pressure;
    const T density;
    const MathVector<T, dimD<T, L>()> velocity;
    const T velocity2;

  public:
  Entropic(const BoundaryPosition& boundaryPosition_in,
           const int start_in, const int end_in,
           const T pressure_in,
           const MathVector<T, dimD<T, L>()>& velocity_in)
    : BoundaryCondition<T, L>(boundaryPosition_in, start_in, end_in)
      , pressure(pressure_in)
      , density(pressure_in*inv_cs2<T, L>())
      , velocity(velocity_in)
      , velocity2(velocity_in.norm2())
      {}

    void apply(Lattice<T, L>& lattice)
    {}
  };

    template<class T, LatticeType L>

  class Corner : public BoundaryCondition<T, L> {
  private:
    const T pressure;
    const T density;
    const MathVector<T, dimD<T, L>()> velocity;

  public:
  Corner(const BoundaryPosition& boundaryPosition_in,
         const T pressure_in,
         const MathVector<T, dimD<T, L>()>& velocity_in)
    : BoundaryCondition<T, L>(boundaryPosition_in, -1, -1)
      , pressure(pressure_in)
      , density(pressure_in*inv_cs2<T, L>())
      , velocity(velocity_in)
    {}

    void apply(Lattice<T, L>& lattice)
    {}
  };

  template<class T, LatticeType L>
    std::shared_ptr<BoundaryCondition<T, L>> Create(const BoundaryType& boundaryType,
                                              const BoundaryPosition& boundaryPosition,
                                              const int start, const int end,
                                              const T pressure,
                                              const MathVector<T, dimD<T, L>()>& velocity) {
    switch(boundaryType){
    case BoundaryType::bounceBack_halfWay:{
      return std::shared_ptr<BoundaryCondition<T, L>>(new BounceBack_HalfWay<T, L>(boundaryPosition,
                                                                                   start, end));
    }
    case BoundaryType::pressure_ZouHe:{
      return std::shared_ptr<BoundaryCondition<T, L>>(new Pressure_ZouHe<T, L>(boundaryPosition,
                                                                               start, end, pressure));
    }
    case BoundaryType::velocity_ZouHe:{
      return std::shared_ptr<BoundaryCondition<T, L>>(new Velocity_ZouHe<T, L>(boundaryPosition,
                                                                               start, end, velocity));
    }
    case BoundaryType::entropic:{
      return std::shared_ptr<BoundaryCondition<T, L>>(new Entropic<T, L>(boundaryPosition,
                                                                         start, end,
                                                                         pressure, velocity));
    }
    case BoundaryType::corner:{
      return std::shared_ptr<BoundaryCondition<T, L>>(new Corner<T, L>(boundaryPosition,
                                                                       pressure, velocity));
    }
    default:
      BOOST_LOG_TRIVIAL(error) << "Error: unknown type of boundary condition.";
      return nullptr;
    }
  }

  template<class T, LatticeType L>
    class BoundaryConditions{
  private:
    std::vector<std::shared_ptr<BoundaryCondition<T, L>> > boundaryConditionsVector;

  public:
    BoundaryConditions(const std::vector<std::shared_ptr<BoundaryCondition<T, L>>>&
                       boundaryConditionsVector_in);

    inline void apply(Lattice<T, L>& l){
      for(auto boundaryCondition : boundaryConditionsVector){
        boundaryCondition->apply(l);
      }
    }
  };

}

#endif // BOUNDARY_H
