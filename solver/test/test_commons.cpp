#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "input.h"
#include "commons.h"

#include <math.h>


using namespace lbm;

BOOST_AUTO_TEST_SUITE(LatticeType)

BOOST_AUTO_TEST_CASE(doubleD1Q3) {
  constexpr ::lbm::LatticeType L = ::lbm::LatticeType::D1Q3;
  typedef double T;

  const auto dimD_ = dimD<T, L>();
  BOOST_CHECK_EQUAL(dimD_, 1);
  const auto dimQ_ = dimQ<T, L>();
  BOOST_CHECK_EQUAL(dimQ_, 3);

  const auto lX_g_ = lX_g<T, L>();
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = lY_g<T, L>();
  BOOST_CHECK_EQUAL(lY_g_, 1);
  const auto lZ_g_ = lZ_g<T, L>();
  BOOST_CHECK_EQUAL(lZ_g_, 1);

  const auto hX_ = hX<T, L>();
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = hY<T, L>();
  BOOST_CHECK_EQUAL(hY_, 0);
  const auto hZ_ = hZ<T, L>();
  BOOST_CHECK_EQUAL(hZ_, 0);

  const auto inv_cs2_ = inv_cs2<T, L>();
  BOOST_CHECK_EQUAL(inv_cs2_, (T)(3.0));
  const auto cs2_ = cs2<T, L>();
  BOOST_CHECK_EQUAL(cs2_, (T)(1.0/3.0));

  const auto w0_ = weight<T, L>(0);
  BOOST_CHECK_EQUAL(w0_, 2.0/3.0);
  const auto w1_ = weight<T, L>(1);
  BOOST_CHECK_EQUAL(w1_, 1.0/3.0);
  const auto w2_ = weight<T, L>(2);
  BOOST_CHECK_EQUAL(w2_, 1.0/3.0);

  const auto c0_0_ = celerity<T, L>(0, 0);
  BOOST_CHECK_EQUAL(c0_0_, (T)0.0);
  const auto c1_0_ = celerity<T, L>(1, 0);
  BOOST_CHECK_EQUAL(c1_0_, (T)-1.0);
  const auto c2_0_ = celerity<T, L>(2, 0);
  BOOST_CHECK_EQUAL(c2_0_, (T)1.0);

}

BOOST_AUTO_TEST_CASE(doubleD2Q9) {
  constexpr ::lbm::LatticeType L = ::lbm::LatticeType::D2Q9;
  typedef double T;

  const auto dimD_ = dimD<T, L>();
  BOOST_CHECK_EQUAL(dimD_, 2);
  const auto dimQ_ = dimQ<T, L>();
  BOOST_CHECK_EQUAL(dimQ_, 9);

  const auto lX_g_ = lX_g<T, L>();
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = lY_g<T, L>();
  BOOST_CHECK_EQUAL(lY_g_, lengthY_g);
  const auto lZ_g_ = lZ_g<T, L>();
  BOOST_CHECK_EQUAL(lZ_g_, 1);

  const auto hX_ = hX<T, L>();
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = hY<T, L>();
  BOOST_CHECK_EQUAL(hY_, 1);
  const auto hZ_ = hZ<T, L>();
  BOOST_CHECK_EQUAL(hZ_, 0);

  const auto inv_cs2_ = inv_cs2<T, L>();
  BOOST_CHECK_EQUAL(inv_cs2_, (T)(3.0));
  const auto cs2_ = cs2<T, L>();
  BOOST_CHECK_EQUAL(cs2_, (T)(1.0/3.0));

  const auto w0_ = weight<T, L>(0);
  BOOST_CHECK_EQUAL(w0_, (T)4.0/9.0);
  const auto w1_ = weight<T, L>(1);
  BOOST_CHECK_EQUAL(w1_, (T)1.0/36.0);
  const auto w2_ = weight<T, L>(2);
  BOOST_CHECK_EQUAL(w2_, (T)1.0/9.0);
  const auto w3_ = weight<T, L>(3);
  BOOST_CHECK_EQUAL(w3_, (T)1.0/36.0);
  const auto w4_ = weight<T, L>(4);
  BOOST_CHECK_EQUAL(w4_, (T)1.0/9.0);
  const auto w5_ = weight<T, L>(5);
  BOOST_CHECK_EQUAL(w5_, (T)1.0/36.0);
  const auto w6_ = weight<T, L>(6);
  BOOST_CHECK_EQUAL(w6_, (T)1.0/9.0);
  const auto w7_ = weight<T, L>(7);
  BOOST_CHECK_EQUAL(w7_, (T)1.0/36.0);
  const auto w8_ = weight<T, L>(8);
  BOOST_CHECK_EQUAL(w8_, (T)1.0/9.0);

  const auto c0_0_ = celerity<T, L>(0, 0);
  BOOST_CHECK_EQUAL(c0_0_, (T)0.0);
  const auto c0_1_ = celerity<T, L>(0, 1);
  BOOST_CHECK_EQUAL(c0_1_, (T)0.0);
  const auto c1_0_ = celerity<T, L>(1, 0);
  BOOST_CHECK_EQUAL(c1_0_, (T)-1.0);
  const auto c1_1_ = celerity<T, L>(1, 1);
  BOOST_CHECK_EQUAL(c1_1_, (T)1.0);
  const auto c2_0_ = celerity<T, L>(2, 0);
  BOOST_CHECK_EQUAL(c2_0_, (T)-1.0);
  const auto c2_1_ = celerity<T, L>(2, 1);
  BOOST_CHECK_EQUAL(c2_1_, (T)0.0);
  const auto c3_0_ = celerity<T, L>(3, 0);
  BOOST_CHECK_EQUAL(c3_0_, (T)-1.0);
  const auto c3_1_ = celerity<T, L>(3, 1);
  BOOST_CHECK_EQUAL(c3_1_, (T)-1.0);
  const auto c4_0_ = celerity<T, L>(4, 0);
  BOOST_CHECK_EQUAL(c4_0_, (T)0.0);
  const auto c4_1_ = celerity<T, L>(4, 1);
  BOOST_CHECK_EQUAL(c4_1_, (T)-1.0);
  const auto c5_0_ = celerity<T, L>(5, 0);
  BOOST_CHECK_EQUAL(c5_0_, (T)1.0);
  const auto c5_1_ = celerity<T, L>(5, 1);
  BOOST_CHECK_EQUAL(c5_1_, (T)-1.0);
  const auto c6_0_ = celerity<T, L>(6, 0);
  BOOST_CHECK_EQUAL(c6_0_, (T)1.0);
  const auto c6_1_ = celerity<T, L>(6, 1);
  BOOST_CHECK_EQUAL(c6_1_, (T)0.0);
  const auto c7_0_ = celerity<T, L>(7, 0);
  BOOST_CHECK_EQUAL(c7_0_, (T)1.0);
  const auto c7_1_ = celerity<T, L>(7, 1);
  BOOST_CHECK_EQUAL(c7_1_, (T)1.0);
  const auto c8_0_ = celerity<T, L>(8, 0);
  BOOST_CHECK_EQUAL(c8_0_, (T)0.0);
  const auto c8_1_ = celerity<T, L>(8, 1);
  BOOST_CHECK_EQUAL(c8_1_, (T)1.0);

}

BOOST_AUTO_TEST_SUITE_END()
