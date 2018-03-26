#ifndef COMPUTATION_H
#define COMPUTATION_H

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"

namespace lbm {

  template<Architecture architecture, unsigned int Dimension>
  class Computation {
  public:
    Computation(const Position& start_in,
                const Position& end_in) {}

    template<typename Callback, typename... Arguments>
    void Do(Callback function, const Arguments&... arguments) {}
  };

  template<unsigned int Dimension>
  class Computation<Architecture::Generic, Dimension> {
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
    : public Computation<Architecture::Generic, 1> {
  private:
    using Computation<Architecture::Generic, 1>::start;
    using Computation<Architecture::Generic, 1>::end;
    using Computation<Architecture::Generic, 1>::length;

  public:
    using Computation<Architecture::Generic, 1>::Computation;

    template<typename Callback, typename... Arguments>
    void Do(Callback function, const Arguments&... arguments) {
      INSTRUMENT_OFF("Computation<Architecture::CPU, 2>::Do<Callback>",1)

      Position iP{{0}};
      for(auto iX = start[d::X]; iX < end[d::X]; ++iX) {
        iP[d::X] = iX;
        function(iP, arguments...);
      }
    }
  };

  template<>
  class Computation<Architecture::CPU, 2>
    : public Computation<Architecture::Generic, 2> {
  private:
    using Computation<Architecture::Generic, 2>::start;
    using Computation<Architecture::Generic, 2>::end;
    using Computation<Architecture::Generic, 2>::length;

  public:
    using Computation<Architecture::Generic, 2>::Computation;

    template<typename Callback, typename... Arguments>
    void Do(Callback function, const Arguments&... arguments) {
    INSTRUMENT_OFF("Computation<Architecture::CPU, 2>::Do<Callback>",2)

    Position iP{{0}};
    for(auto iX = start[d::X]; iX < end[d::X]; ++iX) {
      iP[d::X] = iX;
      for(auto iY = start[d::Y]; iY < end[d::Y]; ++iY) {
        iP[d::Y] = iY;
        function(iP, arguments...);
      }
    }
  }
  };

  template<>
  class Computation<Architecture::CPU, 3>
    : public Computation<Architecture::Generic, 3> {
  private:
    using Computation<Architecture::Generic, 3>::start;
    using Computation<Architecture::Generic, 3>::end;
    using Computation<Architecture::Generic, 3>::length;

  public:
    using Computation<Architecture::Generic, 3>::Computation;

    template<typename Callback, typename... Arguments>
    void Do(Callback function, const Arguments&... arguments) {
      INSTRUMENT_OFF("Computation<Architecture::CPU, 2>::Do<Callback>",3)

      Position iP{{0}};
      for(auto iX = start[d::X]; iX < end[d::X]; ++iX) {
        iP[d::X] = iX;
        for(auto iY = start[d::Y]; iY < end[d::Y]; ++iY) {
          iP[d::Y] = iY;
          for(auto iZ = start[d::Z]; iZ < end[d::Z]; ++iZ) {
            iP[d::Z] = iZ;
            function(iP, arguments...);
          }
        }
      }
    }

  };

}

#endif // COMPUTATION_H
