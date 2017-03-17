#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include <math.h>

#include "output.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(Output)

BOOST_AUTO_TEST_CASE(OutputVTR) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D1Q3;

  typedef Parameters<T, L> Param;

  std::vector<T> fEq(Param::dimQ * Param::lX_g * Param::lY_g * Param::lZ_g, (T)0);

  const int idx = Param::lX_g * Param::lY_g * Param::lZ_g - 1;

  const auto density_ = computeDensity<T, L>(fEq.data(), idx);
  BOOST_CHECK_EQUAL(density_, (T)0);

}

BOOST_AUTO_TEST_SUITE_END()
