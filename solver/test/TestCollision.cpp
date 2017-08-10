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
#include "Collision.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestCollision)

BOOST_AUTO_TEST_CASE(TestBGK) {
  constexpr CollisionType collisionType = CollisionType::BGK;
  valueType f[100]
    Collision<valueType, collisionType> collision(0.6, f);


  BOOST_TEST(collision.getAlpha() == (valueType) 2,
             tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_SUITE_END()
