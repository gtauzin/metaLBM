#ifndef COMPUTATION_H
#define COMPUTATION_H

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"

namespace lbm {

  template<Architecture architecture, unsigned int Dimension>
  class Computation {
  public:
    Computation(const MathVector<unsigned int, 3>& start_in,
                const MathVector<unsigned int, 3>& end_in) {}

    template<typename Callback, typename... Arguments>
    void Do(Callback function, const Arguments&... arguments) {}
  };

  template<unsigned int Dimension>
  class Computation<Architecture::Generic, Dimension> {
  protected:
    const MathVector<unsigned int, 3> start;
    const MathVector<unsigned int, 3> end;
    const MathVector<unsigned int, 3> length;

  public:
      Computation(const MathVector<unsigned int, 3>& start_in,
                  const MathVector<unsigned int, 3>& end_in)
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

      MathVector<unsigned int, 3> iP{{0}};
      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
      #pragma omp simd
      for(iP[d::X] = start[d::X]; iP[d::X] < end[d::X]; ++iP[d::X]) {
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

    MathVector<unsigned int, 3> iP{{0}};
    #pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(iP[d::X] = start[d::X]; iP[d::X] < end[d::X]; ++iP[d::X]) {
      #pragma omp simd
      for(iP[d::Y] = start[d::Y]; iP[d::Y] < end[d::Y]; ++iP[d::Y]) {
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

      MathVector<unsigned int, 3> iP{{0}};
    #pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(iP[d::X] = start[d::X]; iP[d::X] < end[d::X]; ++iP[d::X]) {
      for(iP[d::Y] = start[d::Y]; iP[d::Y] < end[d::Y]; ++iP[d::Y]) {
        #pragma omp simd
        for(iP[d::Z] = start[d::Z]; iP[d::Z] < end[d::Z]; ++iP[d::Z]) {
          function(iP, arguments...);
          }
        }
      }
    }

  };

}

#endif // COMPUTATION_H
