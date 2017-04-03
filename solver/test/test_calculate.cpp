#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include <math.h>

#define NPROCS 1
#include "input.h"
#include "commons.h"
#include "calculate.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(Calculate)

BOOST_AUTO_TEST_CASE(doubleD1Q3_Density) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D1Q3;
  typedef double valueType;
  typedef Parameters<valueType, latticeType> Param;

  std::vector<valueType> fEq(Param::dimQ * Param::lX_g * Param::lY_g * Param::lZ_g, (valueType)0);

  const int idx = Param::lX_g * Param::lY_g * Param::lZ_g - 1;

  const auto density_ = computeDensity<valueType>(fEq.data(), idx);
  BOOST_CHECK_EQUAL(density_, (valueType)0);

}

BOOST_AUTO_TEST_SUITE_END()
