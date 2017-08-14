#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unistd.h>

#include "Algorithm.h"


using namespace lbm;

int main(int argc, char* argv[]) {

  char hostname[1024];
  gethostname(hostname, sizeof(hostname)-1);

  Algorithm_ algorithm(MathVector<unsigned int, 3>{{ 0 }},
                       MathVector<unsigned int, 3>{{ 0 }},
                       std::string(hostname));

  algorithm.computeLBM();

  return 0;
}
