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
#include "init.h"
#include "compute.h"
#include "commons.h"

using namespace lbm;

int main() {
  initLogging(0);

  BOOST_LOG_TRIVIAL(debug) << "Logging for debug starts.";

  Init<double, latticeType> init = init_Simulation<double, latticeType>(0);

  compute<double, latticeType>(init);

  return 0;

}
