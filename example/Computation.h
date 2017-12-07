#ifndef COMPUTATION_H
#define COMPUTATION_H

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
      MathVector<unsigned int, 3> iP;
      for(iP[0] = start[0]; iP[0] < end[0]; ++iP[0]) {
        function(iP);
      }
    }
  };


}


#endif // COMPUTATION_H
