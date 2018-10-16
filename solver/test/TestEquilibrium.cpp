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
#include "Equilibrium.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestEquilibrium)

BOOST_AUTO_TEST_CASE(TestGuo) {
  MathVector<valueType, L::dimD> velocity;
  valueType density = 0;

  Equilibrium<valueType, lbm::LatticeType::D1Q3, EquilibriumType::Basic> equi(density, velocity);

  for(int iQ = 0; iQ < L::dimQ; ++iQ) {
    BOOST_TEST(0 == (valueType) 0,
               tt::tolerance(1e-15));
  }
}

BOOST_AUTO_TEST_SUITE_END()
