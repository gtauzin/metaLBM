#ifndef COMPUTATION_H
#define COMPUTATION_H

#include <omp.h>

#include "Options.h"
#include "Domain.h"

namespace lbm {

  template<Architecture architecture, unsigned int Dimension>
  struct Computation {};


  template<>
  struct Computation<Architecture::CPU, 1> {
  template<typename Callback>
  static void Do(const MathVector<unsigned int, 3>& start,
                 const MathVector<unsigned int, 3>& end,
                 Callback function) {
    SCOREP_INSTRUMENT_ON("Computation<Callback>::Do")

    MathVector<unsigned int, 3> iP;
    #pragma omp parallel for schedule(static) num_threads(NTHREADS)
    #pragma omp simd
    for(unsigned int iX = start[d::X]; iX < end[d::X]; ++iX) {
      iP = {iX, start[d::Y], start[d::Z]};
      function(iP);
      }
    }
  };


  template<>
  struct Computation<Architecture::CPU, 2> {
  template<typename Callback>
  static void Do(const MathVector<unsigned int, 3>& start,
                 const MathVector<unsigned int, 3>& end,
                 Callback function) {
    SCOREP_INSTRUMENT_ON("Computation<Callback>::Do")

      std::cout << "Do called!" << std::endl;
      std::cout << "start   : " << start << std::endl;
      std::cout << "end     : " << end << std::endl;


    MathVector<unsigned int, 3> iP;
    #pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(unsigned int iX = start[d::X]; iX < end[d::X]; ++iX) {
      #pragma omp simd
      for(unsigned int iY = start[d::Y]; iY < end[d::Y]; ++iY) {
        iP = {iX, iY, start[d::Z]};
        function(iP);
      }
    }
  }

  };


  template<>
  struct Computation<Architecture::CPU, 3> {
  template<typename Callback>
  static void Do(const MathVector<unsigned int, 3>& start,
                 const MathVector<unsigned int, 3>& end,
                 Callback function) {
    SCOREP_INSTRUMENT_ON("Computation<Callback>::Do")

    MathVector<unsigned int, 3> iP;
      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(unsigned int iX = start[d::X]; iX < end[d::X]; ++iX) {
      for(unsigned int iY = start[d::Y]; iY < end[d::Y]; ++iY) {
          #pragma omp simd
        for(unsigned int iZ = start[d::Z]; iZ < end[d::Z]; ++iZ) {
            iP = {iX, iY, iZ};
            function(iP);
          }
        }
      }
    }

  };

  typedef Computation<architecture, L::dimD> Computation_;


}


#endif // COMPUTATION_H
