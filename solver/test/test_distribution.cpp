#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
namespace tt = boost::test_tools;

#define NPROCS 1
#include "commons.h"
#include "distribution.h"

#include <math.h>
#include <iostream>

using namespace lbm;

BOOST_AUTO_TEST_SUITE(Distribution)


BOOST_AUTO_TEST_SUITE_END()
