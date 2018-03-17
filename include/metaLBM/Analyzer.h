#ifndef ANALYZER_H
#define ANALYZER_H

#include <cmath>

#ifdef USE_FFTW
  #include <fftw3-mpi.h>
#endif

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "Lattice.h"
#include "Communication.h"
#include "Computation.h"

// TODO: Functors

namespace lbm {

  template<class T, AnalyzerType AnalyzerType, Architecture architecture,
    Implementation implementation, PartitionningType partitionningType>
  class Analyzer {};

  template<class T, Architecture architecture, Implementation implementation,
           PartitionningType partitionningType>
  class Analyzer<T, AnalyzerType::Generic, architecture, implementation,
                 partitionningType>
  {};

}

#endif // ANALYZER_H
