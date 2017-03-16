#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include <math.h>
#include <array>

#include "input.h"
#include "calculate.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(EquilibriumDistribution)

BOOST_AUTO_TEST_CASE(doubleD1Q3_Density) {
  constexpr ::lbm::LatticeType L = ::lbm::LatticeType::D1Q3;
  typedef double T;

  typedef Parameters<T, L> Param;

  std::vector<T> fEq(Param::dimQ * s_g<T, L>(), (T)(0.0));

  const int idx = s_g<T, L>() - 1;

  const auto density_ = computeDensity<T, L>(fEq.data(), idx);
  BOOST_CHECK_EQUAL(density_, (T)0.0);

}

BOOST_AUTO_TEST_SUITE_END()
