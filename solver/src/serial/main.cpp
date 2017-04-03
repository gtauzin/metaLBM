#define NPROCS 1
#define NTHREADS 1

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include "input.h"
#include "lattice.h"
typedef lbm::Lattice<lbm::valueType, lbm::latticeType> L;

#include "init.h"
#include "compute.h"
#include "commons.h"

using namespace lbm;

int main() {
  initLogging(0);

  Init<valueType> init = init_Simulation<double>(0);

  compute<valueType>(init);

  return 0;

}
