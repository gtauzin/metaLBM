#ifndef COMPUTATION_H
#define COMPUTATION_H

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"

namespace lbm {

  template<Architecture architecture, unsigned int Dimension>
  struct Computation {
    template<typename Callback>
    void Do(const MathVector<unsigned int, 3> start,
            const MathVector<unsigned int, 3> end,
            const MathVector<unsigned int, 3> length,
            Callback function);
  };


  template<>
  struct Computation<Architecture::CPU, 1> {
    template<typename Callback>
    void Do(const MathVector<unsigned int, 3> start,
            const MathVector<unsigned int, 3> end,
            const MathVector<unsigned int, 3> length,
            Callback function) {
      SCOREP_INSTRUMENT_ON("Computation<Callback>::Do")

      MathVector<unsigned int, 3> iP;
      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
      #pragma omp simd
      for(iP[d::X] = start[d::X]; iP[d::X] < end[d::X]; ++iP[d::X]) {
        function(iP);
      }
    }
  };

  template<>
  struct Computation<Architecture::CPU, 2> {
  template<typename Callback>
    void Do(const MathVector<unsigned int, 3> start,
            const MathVector<unsigned int, 3> end,
            const MathVector<unsigned int, 3> length,
            Callback function) {
    SCOREP_INSTRUMENT_ON("Computation<Callback>::Do")

    MathVector<unsigned int, 3> iP;
    #pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(iP[d::X] = start[d::X]; iP[d::X] < end[d::X]; ++iP[d::X]) {
      #pragma omp simd
      for(iP[d::Y] = start[d::Y]; iP[d::Y] < end[d::Y]; ++iP[d::Y]) {
        function(iP);
      }
    }
  }
  };

  template<>
  struct Computation<Architecture::CPU, 3> {
  template<typename Callback>
    void Do(const MathVector<unsigned int, 3> start,
            const MathVector<unsigned int, 3> end,
            const MathVector<unsigned int, 3> length,
            Callback function) {
    SCOREP_INSTRUMENT_ON("Computation<Callback>::Do")

    MathVector<unsigned int, 3> iP;
    #pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(iP[d::X] = start[d::X]; iP[d::X] < end[d::X]; ++iP[d::X]) {
      for(iP[d::Y] = start[d::Y]; iP[d::Y] < end[d::Y]; ++iP[d::Y]) {
        #pragma omp simd
        for(iP[d::Z] = start[d::Z]; iP[d::Z] < end[d::Z]; ++iP[d::Z]) {
          function(iP);
          }
        }
      }
    }

  };

}


#endif // COMPUTATION_H
