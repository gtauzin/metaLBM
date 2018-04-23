#pragma once

#include "Commons.h"
#include "MathVector.h"
#include "Options.h"
#include "Stream.h"

namespace lbm {

template <Architecture architecture, unsigned int Dimension>
class Computation {
 protected:
  const Position start;
  const Position end;
  const Position length;
  const Position dir;

 public:
  Computation(const Position& start_in,
              const Position& end_in,
              const Position& dir_in = {{d::X, d::Y, d::Z}})
      : start(start_in), end(end_in), length(end_in - start_in), dir(dir_in) {}
};

template <>
class Computation<Architecture::CPU, 1>
    : public Computation<Architecture::CPU, 0> {
 private:
  using Base = Computation<Architecture::CPU, 0>;

 public:
  using Base::Computation;

  template <typename Callback, typename... Arguments>
  void Do(Callback function, const Arguments... arguments) {
    {LBM_INSTRUMENT_OFF(
        "Computation<Architecture::CPU, 2>::Do<Callback>", 1)}

    Position iP = start;
    for (auto i0 = Base::start[Base::dir[0]]; i0 < Base::end[Base::dir[0]];
         ++i0) {
      iP[Base::dir[0]] = i0;
      function(iP, arguments...);
    }
  }

  template <typename Callback, typename... Arguments>
  void Do(const Stream<Architecture::CPU>& stream,
          Callback function,
          const Arguments... arguments) {
    Do(function, arguments...);
  }
};

template <>
class Computation<Architecture::CPU, 2>
    : public Computation<Architecture::CPU, 0> {
 private:
  using Base = Computation<Architecture::CPU, 0>;

 public:
  using Base::Computation;

  template <typename Callback, typename... Arguments>
  void Do(Callback function, const Arguments... arguments) {
    {LBM_INSTRUMENT_OFF(
        "Computation<Architecture::CPU, 2>::Do<Callback>", 2)}

    Position iP = start;
    for (auto i0 = Base::start[Base::dir[0]]; i0 < Base::end[Base::dir[0]];
         ++i0) {
      iP[Base::dir[0]] = i0;
      for (auto i1 = Base::start[Base::dir[1]]; i1 < Base::end[Base::dir[1]];
           ++i1) {
        iP[Base::dir[1]] = i1;
        function(iP, arguments...);
      }
    }
  }

  template <typename Callback, typename... Arguments>
  void Do(const Stream<Architecture::CPU>& stream,
          Callback function,
          const Arguments... arguments) {
    Do(function, arguments...);
  }
};

template <>
class Computation<Architecture::CPU, 3>
    : public Computation<Architecture::CPU, 0> {
 private:
  using Base = Computation<Architecture::CPU, 0>;

 public:
  using Base::Computation;

  template <typename Callback, typename... Arguments>
  void Do(Callback function, const Arguments... arguments) {
    {LBM_INSTRUMENT_OFF(
        "Computation<Architecture::CPU, 2>::Do<Callback>", 3)}

    Position iP = start;
    for (auto i0 = Base::start[Base::dir[0]]; i0 < Base::end[Base::dir[0]];
         ++i0) {
      iP[Base::dir[0]] = i0;
      for (auto i1 = Base::start[Base::dir[1]]; i1 < Base::end[Base::dir[1]];
           ++i1) {
        iP[Base::dir[1]] = i1;
        for (auto i2 = Base::start[Base::dir[2]]; i2 < Base::end[Base::dir[2]];
             ++i2) {
          iP[Base::dir[2]] = i2;
          function(iP, arguments...);
        }
      }
    }
  }

  template <typename Callback, typename... Arguments>
  void Do(const Stream<Architecture::CPU>& stream,
          Callback function,
          const Arguments... arguments) {
    Do(function, arguments...);
  }
};

}  // namespace lbm
