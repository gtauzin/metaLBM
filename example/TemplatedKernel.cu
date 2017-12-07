#include <iostream>
#include <string>

namespace lbm {
  enum class Architecture {Generic, CPU, GPU};
}

#include "Computation.cuh"
#include "Functor.h"

using namespace lbm;


int main (int argc, char** argv) {
  cudaSetDevice ( 0 );

  Functor<Architecture::GPU> testFunctor;

  testFunctor.compute();

  cudaThreadExit();

  return 0;
}
