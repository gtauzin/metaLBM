#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
namespace tt = boost::test_tools;

#include <iostream>

#define NPROCS 1
#include "Lattice.h"
typedef double valueType;
typedef lbm::Lattice<valueType, lbm::LatticeType::D1Q3> L;
#include "Force.h"
#include "MathVector.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestForce)

BOOST_AUTO_TEST_CASE(TestConstant) {
  constexpr ForceType forceType = ForceType::Constant;
  constexpr MathVector<valueType, 3> amplitude{{(valueType) 1}};
  constexpr MathVector<unsigned int, 3> iP{{1}};

  Force<valueType, forceType> force(amplitude, amplitude);
  force.setForce(iP);

  BOOST_TEST(force.getForce()[d::X] == (valueType) 1,
             tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_SUITE_END()
