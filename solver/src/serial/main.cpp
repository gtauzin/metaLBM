#include <iostream>
#include <string>
#include <unistd.h>

#include "Routine.h"
#include "MathVector.h"

using namespace lbm;

int main(int argc, char* argv[]) {

  char hostname[1024];
  gethostname(hostname, sizeof(hostname)-1);

  MathVector<int, 3> sizeMPI = {1, 1, 1};

  MathVector<int, 3> rankMPI{0, 0, 0};

  Routine_ routine(rankMPI, sizeMPI,
                   std::string(hostname));

  routine.computeLBM();

  return 0;
}
