#pragma once

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "Stream.h"

namespace lbm {

  template<Architecture architecture, unsigned int Dimension>
  class Computation {
  protected:
    const Position start;
    const Position end;
    const Position length;

  public:
    Computation(const Position& start_in,
                const Position& end_in)
      : start(start_in)
      , end(end_in)
      , length(end_in-start_in)
    {}

  };

  template<>
  class Computation<Architecture::CPU, 1>
    : public Computation<Architecture::CPU, 0> {
  private:
   using Base = Computation<Architecture::CPU, 0>;

  public:
   using Base::Computation;
    //using Base::Do;

    template<typename Callback, typename... Arguments>
    void Do(Stream<Architecture::CPU> stream,
            Callback function, const Arguments... arguments) {
      { INSTRUMENT_OFF("Computation<Architecture::CPU, 2>::Do<Callback>",1) }

      Position iP{{0}};
      for(auto iX = Base::start[d::X]; iX < Base::end[d::X]; ++iX) {
        iP[d::X] = iX;
        function(iP, arguments...);
      }
    }

    template<typename Callback, typename... Arguments>
    void Do(Callback function, const Arguments... arguments) {
      Do<Callback>(DefaultStream<Architecture::CPU>(), function, arguments...);
    }

  };

  template<>
  class Computation<Architecture::CPU, 2>
    : public Computation<Architecture::CPU, 0> {
  private:
   using Base = Computation<Architecture::CPU, 0>;

  public:
   using Base::Computation;
    //using Base::Do;

    template<typename Callback, typename... Arguments>
    void Do(Stream<Architecture::CPU> stream,
            Callback function, const Arguments... arguments) {
      { INSTRUMENT_OFF("Computation<Architecture::CPU, 2>::Do<Callback>",2) }

      Position iP{{0}};
      for(auto iX = Base::start[d::X]; iX < Base::end[d::X]; ++iX) {
        iP[d::X] = iX;
        for(auto iY = Base::start[d::Y]; iY < Base::end[d::Y]; ++iY) {
          iP[d::Y] = iY;
          function(iP, arguments...);
        }
      }
    }

    template<typename Callback, typename... Arguments>
    void Do(Callback function, const Arguments... arguments) {
      Do<Callback>(DefaultStream<Architecture::CPU>(), function, arguments...);
    }

  };

  template<>
  class Computation<Architecture::CPU, 3>
    : public Computation<Architecture::CPU, 0> {
  private:
   using Base = Computation<Architecture::CPU, 0>;

  public:
   using Base::Computation;
    //using Base::Do;

    template<typename Callback, typename... Arguments>
    void Do(Stream<Architecture::CPU> stream,
            Callback function, const Arguments... arguments) {
      { INSTRUMENT_OFF("Computation<Architecture::CPU, 2>::Do<Callback>",3) }

      Position iP{{0}};
      for(auto iX = Base::start[d::X]; iX < Base::end[d::X]; ++iX) {
        iP[d::X] = iX;
        for(auto iY = Base::start[d::Y]; iY < Base::end[d::Y]; ++iY) {
          iP[d::Y] = iY;
          for(auto iZ = Base::start[d::Z]; iZ < Base::end[d::Z]; ++iZ) {
            iP[d::Z] = iZ;
            function(iP, arguments...);
          }
        }
      }
    }

    template<typename Callback, typename... Arguments>
    void Do(Callback function, const Arguments... arguments) {
      Do<Callback>(DefaultStream<Architecture::CPU>(), function, arguments...);
    }

  };

}
