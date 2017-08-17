#include <iostream>
#include <string>
#include <unistd.h>

#define NTHREADS 1
#define NPROCS 1
#define _SOA
#define DATA_TYPE double

#include "metaLBM/Routine.h"
#include "metaLBM/MathVector.h"

using namespace lbm;

int main(int argc, char* argv[]) {

  char hostname[1024];
  gethostname(hostname, sizeof(hostname)-1);

  MathVector<int, 3> sizeMPI = {1, 1, 1};

  MathVector<int, 3> rankMPI{0, 0, 0};

  Routine_ routine(rankMPI, sizeMPI,
                   std::string(hostname));

  routine.compute();

  return 0;
}
