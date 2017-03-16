#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "input.h"
#include "commons.h"

#include <math.h>


using namespace lbm;

BOOST_AUTO_TEST_SUITE(LatticeType)

BOOST_AUTO_TEST_CASE(doubleD1Q3) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D1Q3;
  typedef double dataType;

  typedef Parameters<dataType, latticeType> P;

  const auto dimD_ = P::dimD;
  BOOST_CHECK_EQUAL(dimD_, 1);
  const auto dimQ_ = P::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 3);

  const auto lX_g_ = P::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = P::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, 1);
  const auto lZ_g_ = P::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, 1);

  const auto hX_ = P::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = P::hY;
  BOOST_CHECK_EQUAL(hY_, 0);
  const auto hZ_ = P::hZ;
  BOOST_CHECK_EQUAL(hZ_, 0);

  const auto inv_cs2_ = P::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, (T)(3.0));
  const auto cs2_ = P::cs2;
  BOOST_CHECK_EQUAL(cs2_, (T)(1.0/3.0));

  const auto w0_ = P::weight[0];
  BOOST_CHECK_EQUAL(w0_, 2.0/3.0);
  const auto w1_ = P::weight[1];
  BOOST_CHECK_EQUAL(w1_, 1.0/3.0);
  const auto w2_ = P::weight[2];
  BOOST_CHECK_EQUAL(w2_, 1.0/3.0);

  const auto c_ = P::celerity();
  const auto c_r = MathVector<MathVector<T, dimD_>, dimQ_>
    {{
        {{0.0}},
        {{-1.0}},
        {{1.0}}
      }};

  BOOST_CHECK_EQUAL(c_, c_r);


  const auto c0_0_ = P::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, (T)0.0);
  const auto c1_0_ = P::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, (T)-1.0);
  const auto c2_0_ = P::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, (T)1.0);

  const auto c0_ = P::celerity()[0];
  const auto c0_r = MathVector<T, dimD_>{{0.0}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = P::celerity()[1];
  const auto c1_r = MathVector<T, dimD_>{{-1.0}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = P::celerity()[2];
  const auto c2_r = MathVector<T, dimD_>{{1.0}};
  BOOST_CHECK_EQUAL(c2_, c2_r);

}

BOOST_AUTO_TEST_CASE(doubleD2Q9) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D2Q9;
  typedef double dataType;

  typedef Parameters<dataType, latticeType> P;

  const auto dimD_ = P::dimD;
  BOOST_CHECK_EQUAL(dimD_, 2);
  const auto dimQ_ = P::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 9);

  const auto lX_g_ = P::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = P::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, lengthY_g);
  const auto lZ_g_ = P::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, 1);

  const auto hX_ = P::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = P::hY;
  BOOST_CHECK_EQUAL(hY_, 1);
  const auto hZ_ = P::hZ;
  BOOST_CHECK_EQUAL(hZ_, 0);

  const auto inv_cs2_ = P::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, (T)(3.0));
  const auto cs2_ = P::cs2;
  BOOST_CHECK_EQUAL(cs2_, (T)(1.0/3.0));

  const auto w0_ = P::weight[0];
  BOOST_CHECK_EQUAL(w0_, (T)4.0/9.0);
  const auto w1_ = P::weight[1];
  BOOST_CHECK_EQUAL(w1_, (T)1.0/36.0);
  const auto w2_ = P::weight[2];
  BOOST_CHECK_EQUAL(w2_, (T)1.0/9.0);
  const auto w3_ = P::weight[3];
  BOOST_CHECK_EQUAL(w3_, (T)1.0/36.0);
  const auto w4_ = P::weight[4];
  BOOST_CHECK_EQUAL(w4_, (T)1.0/9.0);
  const auto w5_ = P::weight[5];
  BOOST_CHECK_EQUAL(w5_, (T)1.0/36.0);
  const auto w6_ = P::weight[6];
  BOOST_CHECK_EQUAL(w6_, (T)1.0/9.0);
  const auto w7_ = P::weight[7];
  BOOST_CHECK_EQUAL(w7_, (T)1.0/36.0);
  const auto w8_ = P::weight[8];
  BOOST_CHECK_EQUAL(w8_, (T)1.0/9.0);

  const auto c0_ = P::celerity()[0];
  const auto c0_r = MathVector<T, dimD_>{{0.0, 0.0}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = P::celerity()[1];
  const auto c1_r = MathVector<T, dimD_>{{-1.0, 1.0}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = P::celerity()[2];
  const auto c2_r = MathVector<T, dimD_>{{-1.0, 0.0}};
  BOOST_CHECK_EQUAL(c2_, c2_r);
  const auto c3_ = P::celerity()[3];
  const auto c3_r = MathVector<T, dimD_>{{-1.0, -1.0}};
  BOOST_CHECK_EQUAL(c3_, c3_r);
  const auto c4_ = P::celerity()[4];
  const auto c4_r = MathVector<T, dimD_>{{0.0, -1.0}};
  BOOST_CHECK_EQUAL(c4_, c4_r);
  const auto c5_ = P::celerity()[5];
  const auto c5_r = MathVector<T, dimD_>{{1.0, -1.0}};
  BOOST_CHECK_EQUAL(c5_, c5_r);
  const auto c6_ = P::celerity()[6];
  const auto c6_r = MathVector<T, dimD_>{{1.0, 0.0}};
  BOOST_CHECK_EQUAL(c6_, c6_r);
  const auto c7_ = P::celerity()[7];
  const auto c7_r = MathVector<T, dimD_>{{1.0, 1.0}};
  BOOST_CHECK_EQUAL(c7_, c7_r);
  const auto c8_ = P::celerity()[8];
  const auto c8_r = MathVector<T, dimD_>{{0.0, 1.0}};
  BOOST_CHECK_EQUAL(c8_, c8_r);

  const auto c0_0_ = P::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, (T)0.0);
  const auto c0_1_ = P::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_1_, (T)0.0);
  const auto c1_0_ = P::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, (T)-1.0);
  const auto c1_1_ = P::celerity()[1][1];
  BOOST_CHECK_EQUAL(c1_1_, (T)1.0);
  const auto c2_0_ = P::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, (T)-1.0);
  const auto c2_1_ = P::celerity()[2][1];
  BOOST_CHECK_EQUAL(c2_1_, (T)0.0);
  const auto c3_0_ = P::celerity()[3][0];
  BOOST_CHECK_EQUAL(c3_0_, (T)-1.0);
  const auto c3_1_ = P::celerity()[3][1];
  BOOST_CHECK_EQUAL(c3_1_, (T)-1.0);
  const auto c4_0_ = P::celerity()[4][0];
  BOOST_CHECK_EQUAL(c4_0_, (T)0.0);
  const auto c4_1_ = P::celerity()[4][1];
  BOOST_CHECK_EQUAL(c4_1_, (T)-1.0);
  const auto c5_0_ = P::celerity()[5][0];
  BOOST_CHECK_EQUAL(c5_0_, (T)1.0);
  const auto c5_1_ = P::celerity()[5][1];
  BOOST_CHECK_EQUAL(c5_1_, (T)-1.0);
  const auto c6_0_ = P::celerity()[6][0];
  BOOST_CHECK_EQUAL(c6_0_, (T)1.0);
  const auto c6_1_ = P::celerity()[6][1];
  BOOST_CHECK_EQUAL(c6_1_, (T)0.0);
  const auto c7_0_ = P::celerity()[7][0];
  BOOST_CHECK_EQUAL(c7_0_, (T)1.0);
  const auto c7_1_ = P::celerity()[7][1];
  BOOST_CHECK_EQUAL(c7_1_, (T)1.0);
  const auto c8_0_ = P::celerity()[8][0];
  BOOST_CHECK_EQUAL(c8_0_, (T)0.0);
  const auto c8_1_ = P::celerity()[8][1];
  BOOST_CHECK_EQUAL(c8_1_, (T)1.0);

}

BOOST_AUTO_TEST_SUITE_END()
