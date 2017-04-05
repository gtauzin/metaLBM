#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
namespace tt = boost::test_tools;

#define NPROCS 1
#include "lattice.h"

#include <math.h>
#include <iostream>

using namespace lbm;

BOOST_AUTO_TEST_SUITE(Lattice)

BOOST_AUTO_TEST_CASE(isotropyD1Q3) {
  constexpr LatticeType latticeType = ::lbm::LatticeType::D1Q3;
  typedef double valueType;
  typedef ::lbm::Lattice<valueType, latticeType> L;

  valueType sumWeight = 0.0;
  valueType sumWeight_r = 1.0;

  MathVector<valueType, L::dimD> sumWeightCelerity{{0.0}};
  MathVector<valueType, L::dimD> sumWeightCelerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD> sumWeight2Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD> sumWeight2Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD> sumWeight3Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD> sumWeight3Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD> sumWeight4Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD> sumWeight4Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD*L::dimD> sumWeight5Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD*L::dimD> sumWeight5Celerity_r{{0.0}};


  for(int iQ = 0; iQ < L::dimQ; ++iQ) {
    sumWeight += L::weight()[iQ];

    for(int d_a = 0; d_a < L::dimD; ++d_a) {
      int idx_a = d_a;
      sumWeightCelerity[idx_a] += L::weight()[iQ]*L::celerity()[iQ][d_a];

      for(int d_b = 0; d_b < L::dimD; ++d_b) {
        int idx_b = L::dimD * idx_a + d_b;
        sumWeight2Celerity[idx_b] += L::weight()[iQ]*L::celerity()[iQ][d_a]*L::celerity()[iQ][d_b];

        if(d_a == d_b) {
          sumWeight2Celerity_r[idx_b] = L::cs2;
        }

        for(int d_c = 0; d_c < L::dimD; ++d_c) {
          int idx_c = L::dimD * idx_b + d_c;
          sumWeight3Celerity[idx_c] += L::weight()[iQ]*L::celerity()[iQ][d_a]
            *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c];

          for(int d_d = 0; d_d < L::dimD; ++d_d) {
            int idx_d = L::dimD * idx_c + d_d;
            sumWeight4Celerity[idx_d] += L::weight()[iQ]*L::celerity()[iQ][d_a]
              *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c]*L::celerity()[iQ][d_d];

            if((d_a == d_b && d_c == d_d)
               && (d_a == d_c && d_b == d_d)
               && (d_a == d_d && d_b == d_c)) {
              sumWeight4Celerity_r[idx_d] = 3.0*L::cs2*L::cs2;
            }
            else if(((d_a == d_b && d_c == d_d) && (d_a == d_c && d_b == d_d))
                    || ((d_a == d_b && d_c == d_d) && (d_a == d_d && d_b == d_d))
                    || ((d_a == d_d && d_b == d_c) && (d_a == d_d && d_b == d_d))) {
              sumWeight4Celerity_r[idx_d] = 2.0*L::cs2*L::cs2;
            }
            else if((d_a == d_b && d_c == d_d)
                    || (d_a == d_c && d_b == d_d)
                    || (d_a == d_d && d_b == d_c)) {
              sumWeight4Celerity_r[idx_d] = 1.0*L::cs2*L::cs2;
            }

            for(int d_e = 0; d_e < L::dimD; ++d_e) {
              int idx_e = L::dimD * idx_d + d_e;
              sumWeight5Celerity[idx_e] += L::weight()[iQ]*L::celerity()[iQ][d_a]
                *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c]*L::celerity()[iQ][d_d]
                *L::celerity()[iQ][d_e];
            }
          }
        }
      }
    }
  }

  BOOST_TEST(sumWeight == sumWeight_r, tt::tolerance(1e-15));
  for(int i = 0; i < L::dimD; ++i) {
    BOOST_TEST(sumWeightCelerity[i] == sumWeightCelerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight2Celerity[i] == sumWeight2Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight3Celerity[i] == sumWeight3Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight4Celerity[i] == sumWeight4Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight5Celerity[i] == sumWeight5Celerity_r[i], tt::tolerance(1e-15));
  }
}

BOOST_AUTO_TEST_CASE(doubleD1Q3) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D1Q3;
  typedef double valueType;
  typedef ::lbm::Lattice<valueType, latticeType> L;

  const auto dimD_ = L::dimD;
  BOOST_CHECK_EQUAL(dimD_, 1);
  const auto dimQ_ = L::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 3);

  const auto lX_g_ = L::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = L::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, 1);
  const auto lZ_g_ = L::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, 1);

  const auto hX_ = L::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = L::hY;
  BOOST_CHECK_EQUAL(hY_, 0);
  const auto hZ_ = L::hZ;
  BOOST_CHECK_EQUAL(hZ_, 0);

  const auto inv_cs2_ = L::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, 3.0);
  const auto cs2_ = L::cs2;
  BOOST_CHECK_EQUAL(cs2_, 1.0/3.0);

  const auto w0_ = L::weight()[0];
  BOOST_CHECK_EQUAL(w0_, 2.0/3.0);
  const auto w1_ = L::weight()[1];
  BOOST_CHECK_EQUAL(w1_, 1.0/6.0);
  const auto w2_ = L::weight()[2];
  BOOST_CHECK_EQUAL(w2_, 1.0/6.0);

  const auto c_ = L::celerity();
  const auto c_r = MathVector<MathVector<valueType, dimD_>, dimQ_>
    {{
        {{0.0}},
        {{-1.0}},
        {{1.0}}
      }};

  BOOST_CHECK_EQUAL(c_, c_r);


  const auto c0_0_ = L::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, 0.0);
  const auto c1_0_ = L::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, -1.0);
  const auto c2_0_ = L::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, 1.0);

  const auto c0_ = L::celerity()[0];
  const auto c0_r = MathVector<valueType, dimD_>{{0.0}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = L::celerity()[1];
  const auto c1_r = MathVector<valueType, dimD_>{{-1.0}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = L::celerity()[2];
  const auto c2_r = MathVector<valueType, dimD_>{{1.0}};
  BOOST_CHECK_EQUAL(c2_, c2_r);

}

BOOST_AUTO_TEST_CASE(floatD1Q3) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D1Q3;
  typedef float valueType;
  typedef ::lbm::Lattice<valueType, latticeType> L;

  const auto dimD_ = L::dimD;
  BOOST_CHECK_EQUAL(dimD_, 1);
  const auto dimQ_ = L::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 3);

  const auto lX_g_ = L::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = L::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, 1);
  const auto lZ_g_ = L::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, 1);

  const auto hX_ = L::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = L::hY;
  BOOST_CHECK_EQUAL(hY_, 0);
  const auto hZ_ = L::hZ;
  BOOST_CHECK_EQUAL(hZ_, 0);

  const auto inv_cs2_ = L::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, 3.0f);
  const auto cs2_ = L::cs2;
  BOOST_CHECK_EQUAL(cs2_, 1.0f/3.0f);

  const auto w0_ = L::weight()[0];
  BOOST_CHECK_EQUAL(w0_, 2.0f/3.0f);
  const auto w1_ = L::weight()[1];
  BOOST_CHECK_EQUAL(w1_, 1.0f/6.0f);
  const auto w2_ = L::weight()[2];
  BOOST_CHECK_EQUAL(w2_, 1.0f/6.0f);

  const auto c_ = L::celerity();
  const auto c_r = MathVector<MathVector<valueType, dimD_>, dimQ_>
    {{
        {{0.0f}},
        {{-1.0f}},
        {{1.0f}}
      }};

  BOOST_CHECK_EQUAL(c_, c_r);


  const auto c0_0_ = L::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, 0.0f);

  const auto c1_0_ = L::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, -1.0f);

  const auto c2_0_ = L::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, 1.0f);

  const auto c0_ = L::celerity()[0];
  const auto c0_r = MathVector<valueType, dimD_>{{0.0f}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = L::celerity()[1];
  const auto c1_r = MathVector<valueType, dimD_>{{-1.0f}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = L::celerity()[2];
  const auto c2_r = MathVector<valueType, dimD_>{{1.0f}};
  BOOST_CHECK_EQUAL(c2_, c2_r);
}


BOOST_AUTO_TEST_CASE(isotropyD2Q5) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D2Q5;
  typedef double valueType;
  typedef ::lbm::Lattice<valueType, latticeType> L;

  valueType sumWeight = 0.0;
  valueType sumWeight_r = 1.0;

  MathVector<valueType, L::dimD> sumWeightCelerity{{0.0}};
  MathVector<valueType, L::dimD> sumWeightCelerity_r{{0.0}};

  for(int iQ = 0; iQ < L::dimQ; ++iQ) {
    sumWeight += L::weight()[iQ];

    for(int d_a = 0; d_a < L::dimD; ++d_a) {
      int idx_a = d_a;
      sumWeightCelerity[idx_a] += L::weight()[iQ]*L::celerity()[iQ][d_a];
    }
  }

  BOOST_TEST(sumWeight == sumWeight_r, tt::tolerance(1e-15));
  for(int i = 0; i < L::dimD; ++i) {
    BOOST_TEST(sumWeightCelerity[i] == sumWeightCelerity_r[i], tt::tolerance(1e-15));
  }
}

BOOST_AUTO_TEST_CASE(doubleD2Q5) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D2Q5;
  typedef double valueType;

  typedef ::lbm::Lattice<valueType, latticeType> L;

  const auto dimD_ = L::dimD;
  BOOST_CHECK_EQUAL(dimD_, 2);
  const auto dimQ_ = L::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 5);

  const auto lX_g_ = L::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = L::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, lengthY_g);
  const auto lZ_g_ = L::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, 1);

  const auto hX_ = L::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = L::hY;
  BOOST_CHECK_EQUAL(hY_, 1);
  const auto hZ_ = L::hZ;
  BOOST_CHECK_EQUAL(hZ_, 0);

  const auto inv_cs2_ = L::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, 3.0);
  const auto cs2_ = L::cs2;
  BOOST_CHECK_EQUAL(cs2_, (1.0/3.0));

  const auto w0_ = L::weight()[0];
  BOOST_CHECK_EQUAL(w0_, 4.0/6.0);
  const auto w1_ = L::weight()[1];
  BOOST_CHECK_EQUAL(w1_, 1.0/12.0);
  const auto w2_ = L::weight()[2];
  BOOST_CHECK_EQUAL(w2_, 1.0/12.0);
  const auto w3_ = L::weight()[3];
  BOOST_CHECK_EQUAL(w3_, 1.0/12.0);
  const auto w4_ = L::weight()[4];
  BOOST_CHECK_EQUAL(w4_, 1.0/12.0);

  const auto c0_ = L::celerity()[0];
  const auto c0_r = MathVector<valueType, dimD_>{{0.0, 0.0}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = L::celerity()[1];
  const auto c1_r = MathVector<valueType, dimD_>{{-1.0, 0.0}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = L::celerity()[2];
  const auto c2_r = MathVector<valueType, dimD_>{{0.0, -1.0}};
  BOOST_CHECK_EQUAL(c2_, c2_r);
  const auto c3_ = L::celerity()[3];
  const auto c3_r = MathVector<valueType, dimD_>{{1.0, 0.0}};
  BOOST_CHECK_EQUAL(c3_, c3_r);
  const auto c4_ = L::celerity()[4];
  const auto c4_r = MathVector<valueType, dimD_>{{0.0, 1.0}};
  BOOST_CHECK_EQUAL(c4_, c4_r);

  const auto c0_0_ = L::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, 0.0);
  const auto c0_1_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_1_, 0.0);

  const auto c1_0_ = L::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, -1.0);
  const auto c1_1_ = L::celerity()[1][1];
  BOOST_CHECK_EQUAL(c1_1_, 0.0);

  const auto c2_0_ = L::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, 0.0);
  const auto c2_1_ = L::celerity()[2][1];
  BOOST_CHECK_EQUAL(c2_1_, -1.0);

  const auto c3_0_ = L::celerity()[3][0];
  BOOST_CHECK_EQUAL(c3_0_, 1.0);
  const auto c3_1_ = L::celerity()[3][1];
  BOOST_CHECK_EQUAL(c3_1_, 0.0);

  const auto c4_0_ = L::celerity()[4][0];
  BOOST_CHECK_EQUAL(c4_0_, 0.0);
  const auto c4_1_ = L::celerity()[4][1];
  BOOST_CHECK_EQUAL(c4_1_, 1.0);
}

BOOST_AUTO_TEST_CASE(floatD2Q5) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D2Q5;
  typedef float valueType;

  typedef ::lbm::Lattice<valueType, latticeType> L;

  const auto dimD_ = L::dimD;
  BOOST_CHECK_EQUAL(dimD_, 2);
  const auto dimQ_ = L::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 5);

  const auto lX_g_ = L::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = L::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, lengthY_g);
  const auto lZ_g_ = L::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, 1);

  const auto hX_ = L::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = L::hY;
  BOOST_CHECK_EQUAL(hY_, 1);
  const auto hZ_ = L::hZ;
  BOOST_CHECK_EQUAL(hZ_, 0);

  const auto inv_cs2_ = L::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, 3.0f);
  const auto cs2_ = L::cs2;
  BOOST_CHECK_EQUAL(cs2_, (1.0f/3.0f));

  const auto w0_ = L::weight()[0];
  BOOST_CHECK_EQUAL(w0_, 4.0f/6.0f);
  const auto w1_ = L::weight()[1];
  BOOST_CHECK_EQUAL(w1_, 1.0f/12.0f);
  const auto w2_ = L::weight()[2];
  BOOST_CHECK_EQUAL(w2_, 1.0f/12.0f);
  const auto w3_ = L::weight()[3];
  BOOST_CHECK_EQUAL(w3_, 1.0f/12.0f);
  const auto w4_ = L::weight()[4];
  BOOST_CHECK_EQUAL(w4_, 1.0f/12.0f);

  const auto c0_ = L::celerity()[0];
  const auto c0_r = MathVector<valueType, dimD_>{{0.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = L::celerity()[1];
  const auto c1_r = MathVector<valueType, dimD_>{{-1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = L::celerity()[2];
  const auto c2_r = MathVector<valueType, dimD_>{{0.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c2_, c2_r);
  const auto c3_ = L::celerity()[3];
  const auto c3_r = MathVector<valueType, dimD_>{{1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c3_, c3_r);
  const auto c4_ = L::celerity()[4];
  const auto c4_r = MathVector<valueType, dimD_>{{0.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c4_, c4_r);

  const auto c0_0_ = L::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, 0.0f);
  const auto c0_1_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_1_, 0.0f);

  const auto c1_0_ = L::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, -1.0f);
  const auto c1_1_ = L::celerity()[1][1];
  BOOST_CHECK_EQUAL(c1_1_, 0.0f);

  const auto c2_0_ = L::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, 0.0f);
  const auto c2_1_ = L::celerity()[2][1];
  BOOST_CHECK_EQUAL(c2_1_, -1.0f);

  const auto c3_0_ = L::celerity()[3][0];
  BOOST_CHECK_EQUAL(c3_0_, 1.0f);
  const auto c3_1_ = L::celerity()[3][1];
  BOOST_CHECK_EQUAL(c3_1_, 0.0f);

  const auto c4_0_ = L::celerity()[4][0];
  BOOST_CHECK_EQUAL(c4_0_, 0.0f);
  const auto c4_1_ = L::celerity()[4][1];
  BOOST_CHECK_EQUAL(c4_1_, 1.0f);
}

BOOST_AUTO_TEST_CASE(isotropyD2Q9) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D2Q9;
  typedef double valueType;
  typedef ::lbm::Lattice<valueType, latticeType> L;

  valueType sumWeight = 0.0;
  valueType sumWeight_r = 1.0;

  MathVector<valueType, L::dimD> sumWeightCelerity{{0.0}};
  MathVector<valueType, L::dimD> sumWeightCelerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD> sumWeight2Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD> sumWeight2Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD> sumWeight3Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD> sumWeight3Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD> sumWeight4Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD> sumWeight4Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD*L::dimD> sumWeight5Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD*L::dimD> sumWeight5Celerity_r{{0.0}};


  for(int iQ = 0; iQ < L::dimQ; ++iQ) {
    sumWeight += L::weight()[iQ];

    for(int d_a = 0; d_a < L::dimD; ++d_a) {
      int idx_a = d_a;
      sumWeightCelerity[idx_a] += L::weight()[iQ]*L::celerity()[iQ][d_a];

      for(int d_b = 0; d_b < L::dimD; ++d_b) {
        int idx_b = L::dimD * idx_a + d_b;
        sumWeight2Celerity[idx_b] += L::weight()[iQ]*L::celerity()[iQ][d_a]*L::celerity()[iQ][d_b];

        if(d_a == d_b) {
          sumWeight2Celerity_r[idx_b] = L::cs2;
        }

        for(int d_c = 0; d_c < L::dimD; ++d_c) {
          int idx_c = L::dimD * idx_b + d_c;
          sumWeight3Celerity[idx_c] += L::weight()[iQ]*L::celerity()[iQ][d_a]
            *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c];

          for(int d_d = 0; d_d < L::dimD; ++d_d) {
            int idx_d = L::dimD * idx_c + d_d;
            sumWeight4Celerity[idx_d] += L::weight()[iQ]*L::celerity()[iQ][d_a]
              *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c]*L::celerity()[iQ][d_d];

            if((d_a == d_b && d_c == d_d)
               && (d_a == d_c && d_b == d_d)
               && (d_a == d_d && d_b == d_c)) {
              sumWeight4Celerity_r[idx_d] = 3.0*L::cs2*L::cs2;
            }
            else if(((d_a == d_b && d_c == d_d) && (d_a == d_c && d_b == d_d))
                    || ((d_a == d_b && d_c == d_d) && (d_a == d_d && d_b == d_d))
                    || ((d_a == d_d && d_b == d_c) && (d_a == d_d && d_b == d_d))) {
              sumWeight4Celerity_r[idx_d] = 2.0*L::cs2*L::cs2;
            }
            else if((d_a == d_b && d_c == d_d)
                    || (d_a == d_c && d_b == d_d)
                    || (d_a == d_d && d_b == d_c)) {
              sumWeight4Celerity_r[idx_d] = 1.0*L::cs2*L::cs2;
            }

            for(int d_e = 0; d_e < L::dimD; ++d_e) {
              int idx_e = L::dimD * idx_d + d_e;
              sumWeight5Celerity[idx_e] += L::weight()[iQ]*L::celerity()[iQ][d_a]
                *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c]*L::celerity()[iQ][d_d]
                *L::celerity()[iQ][d_e];
            }
          }
        }
      }
    }
  }

  BOOST_TEST(sumWeight == sumWeight_r, tt::tolerance(1e-15));
  for(int i = 0; i < L::dimD; ++i) {
    BOOST_TEST(sumWeightCelerity[i] == sumWeightCelerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight2Celerity[i] == sumWeight2Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight3Celerity[i] == sumWeight3Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight4Celerity[i] == sumWeight4Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight5Celerity[i] == sumWeight5Celerity_r[i], tt::tolerance(1e-15));
  }
}

BOOST_AUTO_TEST_CASE(doubleD2Q9) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D2Q9;
  typedef double valueType;

  typedef ::lbm::Lattice<valueType, latticeType> L;

  const auto dimD_ = L::dimD;
  BOOST_CHECK_EQUAL(dimD_, 2);
  const auto dimQ_ = L::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 9);

  const auto lX_g_ = L::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = L::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, lengthY_g);
  const auto lZ_g_ = L::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, 1);

  const auto hX_ = L::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = L::hY;
  BOOST_CHECK_EQUAL(hY_, 1);
  const auto hZ_ = L::hZ;
  BOOST_CHECK_EQUAL(hZ_, 0);

  const auto inv_cs2_ = L::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, (3.0));
  const auto cs2_ = L::cs2;
  BOOST_CHECK_EQUAL(cs2_, (1.0/3.0));

  const auto w0_ = L::weight()[0];
  BOOST_CHECK_EQUAL(w0_, 4.0/9.0);
  const auto w1_ = L::weight()[1];
  BOOST_CHECK_EQUAL(w1_, 1.0/36.0);
  const auto w2_ = L::weight()[2];
  BOOST_CHECK_EQUAL(w2_, 1.0/9.0);
  const auto w3_ = L::weight()[3];
  BOOST_CHECK_EQUAL(w3_, 1.0/36.0);
  const auto w4_ = L::weight()[4];
  BOOST_CHECK_EQUAL(w4_, 1.0/9.0);
  const auto w5_ = L::weight()[5];
  BOOST_CHECK_EQUAL(w5_, 1.0/36.0);
  const auto w6_ = L::weight()[6];
  BOOST_CHECK_EQUAL(w6_, 1.0/9.0);
  const auto w7_ = L::weight()[7];
  BOOST_CHECK_EQUAL(w7_, 1.0/36.0);
  const auto w8_ = L::weight()[8];
  BOOST_CHECK_EQUAL(w8_, 1.0/9.0);

  const auto c0_ = L::celerity()[0];
  const auto c0_r = MathVector<valueType, dimD_>{{0.0, 0.0}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = L::celerity()[1];
  const auto c1_r = MathVector<valueType, dimD_>{{-1.0, 1.0}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = L::celerity()[2];
  const auto c2_r = MathVector<valueType, dimD_>{{-1.0, 0.0}};
  BOOST_CHECK_EQUAL(c2_, c2_r);
  const auto c3_ = L::celerity()[3];
  const auto c3_r = MathVector<valueType, dimD_>{{-1.0, -1.0}};
  BOOST_CHECK_EQUAL(c3_, c3_r);
  const auto c4_ = L::celerity()[4];
  const auto c4_r = MathVector<valueType, dimD_>{{0.0, -1.0}};
  BOOST_CHECK_EQUAL(c4_, c4_r);
  const auto c5_ = L::celerity()[5];
  const auto c5_r = MathVector<valueType, dimD_>{{1.0, -1.0}};
  BOOST_CHECK_EQUAL(c5_, c5_r);
  const auto c6_ = L::celerity()[6];
  const auto c6_r = MathVector<valueType, dimD_>{{1.0, 0.0}};
  BOOST_CHECK_EQUAL(c6_, c6_r);
  const auto c7_ = L::celerity()[7];
  const auto c7_r = MathVector<valueType, dimD_>{{1.0, 1.0}};
  BOOST_CHECK_EQUAL(c7_, c7_r);
  const auto c8_ = L::celerity()[8];
  const auto c8_r = MathVector<valueType, dimD_>{{0.0, 1.0}};
  BOOST_CHECK_EQUAL(c8_, c8_r);

  const auto c0_0_ = L::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, 0.0);
  const auto c0_1_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_1_, 0.0);

  const auto c1_0_ = L::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, -1.0);
  const auto c1_1_ = L::celerity()[1][1];
  BOOST_CHECK_EQUAL(c1_1_, 1.0);

  const auto c2_0_ = L::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, -1.0);
  const auto c2_1_ = L::celerity()[2][1];
  BOOST_CHECK_EQUAL(c2_1_, 0.0);

  const auto c3_0_ = L::celerity()[3][0];
  BOOST_CHECK_EQUAL(c3_0_, -1.0);
  const auto c3_1_ = L::celerity()[3][1];
  BOOST_CHECK_EQUAL(c3_1_, -1.0);

  const auto c4_0_ = L::celerity()[4][0];
  BOOST_CHECK_EQUAL(c4_0_, 0.0);
  const auto c4_1_ = L::celerity()[4][1];
  BOOST_CHECK_EQUAL(c4_1_, -1.0);

  const auto c5_0_ = L::celerity()[5][0];
  BOOST_CHECK_EQUAL(c5_0_, 1.0);
  const auto c5_1_ = L::celerity()[5][1];
  BOOST_CHECK_EQUAL(c5_1_, -1.0);

  const auto c6_0_ = L::celerity()[6][0];
  BOOST_CHECK_EQUAL(c6_0_, 1.0);
  const auto c6_1_ = L::celerity()[6][1];
  BOOST_CHECK_EQUAL(c6_1_, 0.0);

  const auto c7_0_ = L::celerity()[7][0];
  BOOST_CHECK_EQUAL(c7_0_, 1.0);
  const auto c7_1_ = L::celerity()[7][1];
  BOOST_CHECK_EQUAL(c7_1_, 1.0);

  const auto c8_0_ = L::celerity()[8][0];
  BOOST_CHECK_EQUAL(c8_0_, 0.0);
  const auto c8_1_ = L::celerity()[8][1];
  BOOST_CHECK_EQUAL(c8_1_, 1.0);
}

BOOST_AUTO_TEST_CASE(floatD2Q9) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D2Q9;
  typedef float valueType;

  typedef ::lbm::Lattice<valueType, latticeType> L;

  const auto dimD_ = L::dimD;
  BOOST_CHECK_EQUAL(dimD_, 2);
  const auto dimQ_ = L::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 9);

  const auto lX_g_ = L::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = L::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, lengthY_g);
  const auto lZ_g_ = L::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, 1);

  const auto hX_ = L::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = L::hY;
  BOOST_CHECK_EQUAL(hY_, 1);
  const auto hZ_ = L::hZ;
  BOOST_CHECK_EQUAL(hZ_, 0);

  const auto inv_cs2_ = L::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, (3.0f));
  const auto cs2_ = L::cs2;
  BOOST_CHECK_EQUAL(cs2_, (1.0f/3.0f));

  const auto w0_ = L::weight()[0];
  BOOST_CHECK_EQUAL(w0_, 4.0f/9.0f);
  const auto w1_ = L::weight()[1];
  BOOST_CHECK_EQUAL(w1_, 1.0f/36.0f);
  const auto w2_ = L::weight()[2];
  BOOST_CHECK_EQUAL(w2_, 1.0f/9.0f);
  const auto w3_ = L::weight()[3];
  BOOST_CHECK_EQUAL(w3_, 1.0f/36.0f);
  const auto w4_ = L::weight()[4];
  BOOST_CHECK_EQUAL(w4_, 1.0f/9.0f);
  const auto w5_ = L::weight()[5];
  BOOST_CHECK_EQUAL(w5_, 1.0f/36.0f);
  const auto w6_ = L::weight()[6];
  BOOST_CHECK_EQUAL(w6_, 1.0f/9.0f);
  const auto w7_ = L::weight()[7];
  BOOST_CHECK_EQUAL(w7_, 1.0f/36.0f);
  const auto w8_ = L::weight()[8];
  BOOST_CHECK_EQUAL(w8_, 1.0f/9.0f);

  const auto c0_ = L::celerity()[0];
  const auto c0_r = MathVector<valueType, dimD_>{{0.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = L::celerity()[1];
  const auto c1_r = MathVector<valueType, dimD_>{{-1.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = L::celerity()[2];
  const auto c2_r = MathVector<valueType, dimD_>{{-1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c2_, c2_r);
  const auto c3_ = L::celerity()[3];
  const auto c3_r = MathVector<valueType, dimD_>{{-1.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c3_, c3_r);
  const auto c4_ = L::celerity()[4];
  const auto c4_r = MathVector<valueType, dimD_>{{0.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c4_, c4_r);
  const auto c5_ = L::celerity()[5];
  const auto c5_r = MathVector<valueType, dimD_>{{1.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c5_, c5_r);
  const auto c6_ = L::celerity()[6];
  const auto c6_r = MathVector<valueType, dimD_>{{1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c6_, c6_r);
  const auto c7_ = L::celerity()[7];
  const auto c7_r = MathVector<valueType, dimD_>{{1.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c7_, c7_r);
  const auto c8_ = L::celerity()[8];
  const auto c8_r = MathVector<valueType, dimD_>{{0.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c8_, c8_r);

  const auto c0_0_ = L::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, 0.0f);
  const auto c0_1_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_1_, 0.0f);

  const auto c1_0_ = L::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, -1.0f);
  const auto c1_1_ = L::celerity()[1][1];
  BOOST_CHECK_EQUAL(c1_1_, 1.0f);

  const auto c2_0_ = L::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, -1.0f);
  const auto c2_1_ = L::celerity()[2][1];
  BOOST_CHECK_EQUAL(c2_1_, 0.0f);

  const auto c3_0_ = L::celerity()[3][0];
  BOOST_CHECK_EQUAL(c3_0_, -1.0f);
  const auto c3_1_ = L::celerity()[3][1];
  BOOST_CHECK_EQUAL(c3_1_, -1.0f);

  const auto c4_0_ = L::celerity()[4][0];
  BOOST_CHECK_EQUAL(c4_0_, 0.0f);
  const auto c4_1_ = L::celerity()[4][1];
  BOOST_CHECK_EQUAL(c4_1_, -1.0f);

  const auto c5_0_ = L::celerity()[5][0];
  BOOST_CHECK_EQUAL(c5_0_, 1.0f);
  const auto c5_1_ = L::celerity()[5][1];
  BOOST_CHECK_EQUAL(c5_1_, -1.0f);

  const auto c6_0_ = L::celerity()[6][0];
  BOOST_CHECK_EQUAL(c6_0_, 1.0f);
  const auto c6_1_ = L::celerity()[6][1];
  BOOST_CHECK_EQUAL(c6_1_, 0.0f);

  const auto c7_0_ = L::celerity()[7][0];
  BOOST_CHECK_EQUAL(c7_0_, 1.0f);
  const auto c7_1_ = L::celerity()[7][1];
  BOOST_CHECK_EQUAL(c7_1_, 1.0f);

  const auto c8_0_ = L::celerity()[8][0];
  BOOST_CHECK_EQUAL(c8_0_, 0.0f);
  const auto c8_1_ = L::celerity()[8][1];
  BOOST_CHECK_EQUAL(c8_1_, 1.0f);
}

BOOST_AUTO_TEST_CASE(isotropyD3Q15) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D3Q15;
  typedef double valueType;
  typedef ::lbm::Lattice<valueType, latticeType> L;

  valueType sumWeight = 0.0;
  valueType sumWeight_r = 1.0;

  MathVector<valueType, L::dimD> sumWeightCelerity{{0.0}};
  MathVector<valueType, L::dimD> sumWeightCelerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD> sumWeight2Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD> sumWeight2Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD> sumWeight3Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD> sumWeight3Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD> sumWeight4Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD> sumWeight4Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD*L::dimD> sumWeight5Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD*L::dimD> sumWeight5Celerity_r{{0.0}};


  for(int iQ = 0; iQ < L::dimQ; ++iQ) {
    sumWeight += L::weight()[iQ];

    for(int d_a = 0; d_a < L::dimD; ++d_a) {
      int idx_a = d_a;
      sumWeightCelerity[idx_a] += L::weight()[iQ]*L::celerity()[iQ][d_a];

      for(int d_b = 0; d_b < L::dimD; ++d_b) {
        int idx_b = L::dimD * idx_a + d_b;
        sumWeight2Celerity[idx_b] += L::weight()[iQ]*L::celerity()[iQ][d_a]*L::celerity()[iQ][d_b];

        if(d_a == d_b) {
          sumWeight2Celerity_r[idx_b] = L::cs2;
        }

        for(int d_c = 0; d_c < L::dimD; ++d_c) {
          int idx_c = L::dimD * idx_b + d_c;
          sumWeight3Celerity[idx_c] += L::weight()[iQ]*L::celerity()[iQ][d_a]
            *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c];

          for(int d_d = 0; d_d < L::dimD; ++d_d) {
            int idx_d = L::dimD * idx_c + d_d;
            sumWeight4Celerity[idx_d] += L::weight()[iQ]*L::celerity()[iQ][d_a]
              *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c]*L::celerity()[iQ][d_d];

            if((d_a == d_b && d_c == d_d)
               && (d_a == d_c && d_b == d_d)
               && (d_a == d_d && d_b == d_c)) {
              sumWeight4Celerity_r[idx_d] = 3.0*L::cs2*L::cs2;
            }
            else if(((d_a == d_b && d_c == d_d) && (d_a == d_c && d_b == d_d))
                    || ((d_a == d_b && d_c == d_d) && (d_a == d_d && d_b == d_d))
                    || ((d_a == d_d && d_b == d_c) && (d_a == d_d && d_b == d_d))) {
              sumWeight4Celerity_r[idx_d] = 2.0*L::cs2*L::cs2;
            }
            else if((d_a == d_b && d_c == d_d)
                    || (d_a == d_c && d_b == d_d)
                    || (d_a == d_d && d_b == d_c)) {
              sumWeight4Celerity_r[idx_d] = 1.0*L::cs2*L::cs2;
            }

            for(int d_e = 0; d_e < L::dimD; ++d_e) {
              int idx_e = L::dimD * idx_d + d_e;
              sumWeight5Celerity[idx_e] += L::weight()[iQ]*L::celerity()[iQ][d_a]
                *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c]*L::celerity()[iQ][d_d]
                *L::celerity()[iQ][d_e];
            }
          }
        }
      }
    }
  }

  BOOST_TEST(sumWeight == sumWeight_r, tt::tolerance(1e-15));
  for(int i = 0; i < L::dimD; ++i) {
    BOOST_TEST(sumWeightCelerity[i] == sumWeightCelerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight2Celerity[i] == sumWeight2Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight3Celerity[i] == sumWeight3Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight4Celerity[i] == sumWeight4Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight5Celerity[i] == sumWeight5Celerity_r[i], tt::tolerance(1e-15));
  }
}


BOOST_AUTO_TEST_CASE(doubleD3Q15) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D3Q15;
  typedef double valueType;

  typedef ::lbm::Lattice<valueType, latticeType> L;

  const auto dimD_ = L::dimD;
  BOOST_CHECK_EQUAL(dimD_, 3);
  const auto dimQ_ = L::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 15);

  const auto lX_g_ = L::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = L::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, lengthY_g);
  const auto lZ_g_ = L::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, lengthZ_g);

  const auto hX_ = L::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = L::hY;
  BOOST_CHECK_EQUAL(hY_, 1);
  const auto hZ_ = L::hZ;
  BOOST_CHECK_EQUAL(hZ_, 1);

  const auto inv_cs2_ = L::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, (3.0));
  const auto cs2_ = L::cs2;
  BOOST_CHECK_EQUAL(cs2_, (1.0/3.0));

  const auto w0_ = L::weight()[0];
  BOOST_CHECK_EQUAL(w0_, 2.0/9.0);
  const auto w1_ = L::weight()[1];
  BOOST_CHECK_EQUAL(w1_, 1.0/9.0);
  const auto w2_ = L::weight()[2];
  BOOST_CHECK_EQUAL(w2_, 1.0/9.0);
  const auto w3_ = L::weight()[3];
  BOOST_CHECK_EQUAL(w3_, 1.0/9.0);
  const auto w4_ = L::weight()[4];
  BOOST_CHECK_EQUAL(w4_, 1.0/72.0);
  const auto w5_ = L::weight()[5];
  BOOST_CHECK_EQUAL(w5_, 1.0/72.0);
  const auto w6_ = L::weight()[6];
  BOOST_CHECK_EQUAL(w6_, 1.0/72.0);
  const auto w7_ = L::weight()[7];
  BOOST_CHECK_EQUAL(w7_, 1.0/72.0);
  const auto w8_ = L::weight()[8];
  BOOST_CHECK_EQUAL(w8_, 1.0/9.0);
  const auto w9_ = L::weight()[9];
  BOOST_CHECK_EQUAL(w9_, 1.0/9.0);
  const auto w10_ = L::weight()[10];
  BOOST_CHECK_EQUAL(w10_, 1.0/9.0);
  const auto w11_ = L::weight()[11];
  BOOST_CHECK_EQUAL(w11_, 1.0/72.0);
  const auto w12_ = L::weight()[12];
  BOOST_CHECK_EQUAL(w12_, 1.0/72.0);
  const auto w13_ = L::weight()[13];
  BOOST_CHECK_EQUAL(w13_, 1.0/72.0);
  const auto w14_ = L::weight()[14];
  BOOST_CHECK_EQUAL(w14_, 1.0/72.0);

  const auto c0_ = L::celerity()[0];
  const auto c0_r = MathVector<valueType, dimD_>{{0.0, 0.0, 0.0}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = L::celerity()[1];
  const auto c1_r = MathVector<valueType, dimD_>{{-1.0, 0.0, 0.0}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = L::celerity()[2];
  const auto c2_r = MathVector<valueType, dimD_>{{0.0, -1.0, 0.0}};
  BOOST_CHECK_EQUAL(c2_, c2_r);
  const auto c3_ = L::celerity()[3];
  const auto c3_r = MathVector<valueType, dimD_>{{0.0, 0.0, -1.0}};
  BOOST_CHECK_EQUAL(c3_, c3_r);
  const auto c4_ = L::celerity()[4];
  const auto c4_r = MathVector<valueType, dimD_>{{-1.0, -1.0, -1.0}};
  BOOST_CHECK_EQUAL(c4_, c4_r);
  const auto c5_ = L::celerity()[5];
  const auto c5_r = MathVector<valueType, dimD_>{{-1.0, -1.0, 1.0}};
  BOOST_CHECK_EQUAL(c5_, c5_r);
  const auto c6_ = L::celerity()[6];
  const auto c6_r = MathVector<valueType, dimD_>{{-1.0, 1.0, -1.0}};
  BOOST_CHECK_EQUAL(c6_, c6_r);
  const auto c7_ = L::celerity()[7];
  const auto c7_r = MathVector<valueType, dimD_>{{-1.0, 1.0, 1.0}};
  BOOST_CHECK_EQUAL(c7_, c7_r);
  const auto c8_ = L::celerity()[8];
  const auto c8_r = MathVector<valueType, dimD_>{{1.0, 0.0, 0.0}};
  BOOST_CHECK_EQUAL(c8_, c8_r);
  const auto c9_ = L::celerity()[9];
  const auto c9_r = MathVector<valueType, dimD_>{{0.0, 1.0, 0.0}};
  BOOST_CHECK_EQUAL(c9_, c9_r);
  const auto c10_ = L::celerity()[10];
  const auto c10_r = MathVector<valueType, dimD_>{{0.0, 0.0, 1.0}};
  BOOST_CHECK_EQUAL(c10_, c10_r);
  const auto c11_ = L::celerity()[11];
  const auto c11_r = MathVector<valueType, dimD_>{{1.0, 1.0, 1.0}};
  BOOST_CHECK_EQUAL(c11_, c11_r);
  const auto c12_ = L::celerity()[12];
  const auto c12_r = MathVector<valueType, dimD_>{{1.0, 1.0, -1.0}};
  BOOST_CHECK_EQUAL(c12_, c12_r);
  const auto c13_ = L::celerity()[13];
  const auto c13_r = MathVector<valueType, dimD_>{{1.0, -1.0, 1.0}};
  BOOST_CHECK_EQUAL(c13_, c13_r);
  const auto c14_ = L::celerity()[14];
  const auto c14_r = MathVector<valueType, dimD_>{{1.0, -1.0, -1.0}};
  BOOST_CHECK_EQUAL(c14_, c14_r);

  const auto c0_0_ = L::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, 0.0);
  const auto c0_1_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_1_, 0.0);
  const auto c0_2_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_2_, 0.0);

  const auto c1_0_ = L::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, -1.0);
  const auto c1_1_ = L::celerity()[1][1];
  BOOST_CHECK_EQUAL(c1_1_, 0.0);
  const auto c1_2_ = L::celerity()[1][2];
  BOOST_CHECK_EQUAL(c1_2_, 0.0);

  const auto c2_0_ = L::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, 0.0);
  const auto c2_1_ = L::celerity()[2][1];
  BOOST_CHECK_EQUAL(c2_1_, -1.0);
  const auto c2_2_ = L::celerity()[2][2];
  BOOST_CHECK_EQUAL(c2_2_, 0.0);

  const auto c3_0_ = L::celerity()[3][0];
  BOOST_CHECK_EQUAL(c3_0_, 0.0);
  const auto c3_1_ = L::celerity()[3][1];
  BOOST_CHECK_EQUAL(c3_1_, 0.0);
  const auto c3_2_ = L::celerity()[3][2];
  BOOST_CHECK_EQUAL(c3_2_, -1.0);

  const auto c4_0_ = L::celerity()[4][0];
  BOOST_CHECK_EQUAL(c4_0_, -1.0);
  const auto c4_1_ = L::celerity()[4][1];
  BOOST_CHECK_EQUAL(c4_1_, -1.0);
  const auto c4_2_ = L::celerity()[4][2];
  BOOST_CHECK_EQUAL(c4_2_, -1.0);

  const auto c5_0_ = L::celerity()[5][0];
  BOOST_CHECK_EQUAL(c5_0_, -1.0);
  const auto c5_1_ = L::celerity()[5][1];
  BOOST_CHECK_EQUAL(c5_1_, -1.0);
  const auto c5_2_ = L::celerity()[5][2];
  BOOST_CHECK_EQUAL(c5_2_, 1.0);

  const auto c6_0_ = L::celerity()[6][0];
  BOOST_CHECK_EQUAL(c6_0_, -1.0);
  const auto c6_1_ = L::celerity()[6][1];
  BOOST_CHECK_EQUAL(c6_1_, 1.0);
  const auto c6_2_ = L::celerity()[6][2];
  BOOST_CHECK_EQUAL(c6_2_, -1.0);

  const auto c7_0_ = L::celerity()[7][0];
  BOOST_CHECK_EQUAL(c7_0_, -1.0);
  const auto c7_1_ = L::celerity()[7][1];
  BOOST_CHECK_EQUAL(c7_1_, 1.0);
  const auto c7_2_ = L::celerity()[7][2];
  BOOST_CHECK_EQUAL(c7_2_, 1.0);

  const auto c8_0_ = L::celerity()[8][0];
  BOOST_CHECK_EQUAL(c8_0_, 1.0);
  const auto c8_1_ = L::celerity()[8][1];
  BOOST_CHECK_EQUAL(c8_1_, 0.0);
  const auto c8_2_ = L::celerity()[8][2];
  BOOST_CHECK_EQUAL(c8_2_, 0.0);

  const auto c9_0_ = L::celerity()[9][0];
  BOOST_CHECK_EQUAL(c9_0_, 0.0);
  const auto c9_1_ = L::celerity()[9][1];
  BOOST_CHECK_EQUAL(c9_1_, 1.0);
  const auto c9_2_ = L::celerity()[9][2];
  BOOST_CHECK_EQUAL(c9_2_, 0.0);

  const auto c10_0_ = L::celerity()[10][0];
  BOOST_CHECK_EQUAL(c10_0_, 0.0);
  const auto c10_1_ = L::celerity()[10][1];
  BOOST_CHECK_EQUAL(c10_1_, 0.0);
  const auto c10_2_ = L::celerity()[10][2];
  BOOST_CHECK_EQUAL(c10_2_, 1.0);

  const auto c11_0_ = L::celerity()[11][0];
  BOOST_CHECK_EQUAL(c11_0_, 1.0);
  const auto c11_1_ = L::celerity()[11][1];
  BOOST_CHECK_EQUAL(c11_1_, 1.0);
  const auto c11_2_ = L::celerity()[11][2];
  BOOST_CHECK_EQUAL(c11_2_, 1.0);

  const auto c12_0_ = L::celerity()[12][0];
  BOOST_CHECK_EQUAL(c12_0_, 1.0);
  const auto c12_1_ = L::celerity()[12][1];
  BOOST_CHECK_EQUAL(c12_1_, 1.0);
  const auto c12_2_ = L::celerity()[12][2];
  BOOST_CHECK_EQUAL(c12_2_, -1.0);

  const auto c13_0_ = L::celerity()[13][0];
  BOOST_CHECK_EQUAL(c13_0_, 1.0);
  const auto c13_1_ = L::celerity()[13][1];
  BOOST_CHECK_EQUAL(c13_1_, -1.0);
  const auto c13_2_ = L::celerity()[13][2];
  BOOST_CHECK_EQUAL(c13_2_, 1.0);

  const auto c14_0_ = L::celerity()[14][0];
  BOOST_CHECK_EQUAL(c14_0_, 1.0);
  const auto c14_1_ = L::celerity()[14][1];
  BOOST_CHECK_EQUAL(c14_1_, -1.0);
  const auto c14_2_ = L::celerity()[14][2];
  BOOST_CHECK_EQUAL(c14_2_, -1.0);
}

BOOST_AUTO_TEST_CASE(floatD3Q15) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D3Q15;
  typedef float valueType;

  typedef ::lbm::Lattice<valueType, latticeType> L;

  const auto dimD_ = L::dimD;
  BOOST_CHECK_EQUAL(dimD_, 3);
  const auto dimQ_ = L::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 15);

  const auto lX_g_ = L::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = L::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, lengthY_g);
  const auto lZ_g_ = L::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, lengthZ_g);

  const auto hX_ = L::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = L::hY;
  BOOST_CHECK_EQUAL(hY_, 1);
  const auto hZ_ = L::hZ;
  BOOST_CHECK_EQUAL(hZ_, 1);

  const auto inv_cs2_ = L::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, (3.0f));
  const auto cs2_ = L::cs2;
  BOOST_CHECK_EQUAL(cs2_, (1.0f/3.0f));

  const auto w0_ = L::weight()[0];
  BOOST_CHECK_EQUAL(w0_, 2.0f/9.0f);
  const auto w1_ = L::weight()[1];
  BOOST_CHECK_EQUAL(w1_, 1.0f/9.0f);
  const auto w2_ = L::weight()[2];
  BOOST_CHECK_EQUAL(w2_, 1.0f/9.0f);
  const auto w3_ = L::weight()[3];
  BOOST_CHECK_EQUAL(w3_, 1.0f/9.0f);
  const auto w4_ = L::weight()[4];
  BOOST_CHECK_EQUAL(w4_, 1.0f/72.0f);
  const auto w5_ = L::weight()[5];
  BOOST_CHECK_EQUAL(w5_, 1.0f/72.0f);
  const auto w6_ = L::weight()[6];
  BOOST_CHECK_EQUAL(w6_, 1.0f/72.0f);
  const auto w7_ = L::weight()[7];
  BOOST_CHECK_EQUAL(w7_, 1.0f/72.0f);
  const auto w8_ = L::weight()[8];
  BOOST_CHECK_EQUAL(w8_, 1.0f/9.0f);
  const auto w9_ = L::weight()[9];
  BOOST_CHECK_EQUAL(w9_, 1.0f/9.0f);
  const auto w10_ = L::weight()[10];
  BOOST_CHECK_EQUAL(w10_, 1.0f/9.0f);
  const auto w11_ = L::weight()[11];
  BOOST_CHECK_EQUAL(w11_, 1.0f/72.0f);
  const auto w12_ = L::weight()[12];
  BOOST_CHECK_EQUAL(w12_, 1.0f/72.0f);
  const auto w13_ = L::weight()[13];
  BOOST_CHECK_EQUAL(w13_, 1.0f/72.0f);
  const auto w14_ = L::weight()[14];
  BOOST_CHECK_EQUAL(w14_, 1.0f/72.0f);

  const auto c0_ = L::celerity()[0];
  const auto c0_r = MathVector<valueType, dimD_>{{0.0f, 0.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = L::celerity()[1];
  const auto c1_r = MathVector<valueType, dimD_>{{-1.0f, 0.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = L::celerity()[2];
  const auto c2_r = MathVector<valueType, dimD_>{{0.0f, -1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c2_, c2_r);
  const auto c3_ = L::celerity()[3];
  const auto c3_r = MathVector<valueType, dimD_>{{0.0f, 0.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c3_, c3_r);
  const auto c4_ = L::celerity()[4];
  const auto c4_r = MathVector<valueType, dimD_>{{-1.0f, -1.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c4_, c4_r);
  const auto c5_ = L::celerity()[5];
  const auto c5_r = MathVector<valueType, dimD_>{{-1.0f, -1.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c5_, c5_r);
  const auto c6_ = L::celerity()[6];
  const auto c6_r = MathVector<valueType, dimD_>{{-1.0f, 1.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c6_, c6_r);
  const auto c7_ = L::celerity()[7];
  const auto c7_r = MathVector<valueType, dimD_>{{-1.0f, 1.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c7_, c7_r);
  const auto c8_ = L::celerity()[8];
  const auto c8_r = MathVector<valueType, dimD_>{{1.0f, 0.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c8_, c8_r);
  const auto c9_ = L::celerity()[9];
  const auto c9_r = MathVector<valueType, dimD_>{{0.0f, 1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c9_, c9_r);
  const auto c10_ = L::celerity()[10];
  const auto c10_r = MathVector<valueType, dimD_>{{0.0f, 0.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c10_, c10_r);
  const auto c11_ = L::celerity()[11];
  const auto c11_r = MathVector<valueType, dimD_>{{1.0f, 1.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c11_, c11_r);
  const auto c12_ = L::celerity()[12];
  const auto c12_r = MathVector<valueType, dimD_>{{1.0f, 1.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c12_, c12_r);
  const auto c13_ = L::celerity()[13];
  const auto c13_r = MathVector<valueType, dimD_>{{1.0f, -1.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c13_, c13_r);
  const auto c14_ = L::celerity()[14];
  const auto c14_r = MathVector<valueType, dimD_>{{1.0f, -1.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c14_, c14_r);

  const auto c0_0_ = L::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, 0.0f);
  const auto c0_1_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_1_, 0.0f);
  const auto c0_2_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_2_, 0.0f);

  const auto c1_0_ = L::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, -1.0f);
  const auto c1_1_ = L::celerity()[1][1];
  BOOST_CHECK_EQUAL(c1_1_, 0.0f);
  const auto c1_2_ = L::celerity()[1][2];
  BOOST_CHECK_EQUAL(c1_2_, 0.0f);

  const auto c2_0_ = L::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, 0.0f);
  const auto c2_1_ = L::celerity()[2][1];
  BOOST_CHECK_EQUAL(c2_1_, -1.0f);
  const auto c2_2_ = L::celerity()[2][2];
  BOOST_CHECK_EQUAL(c2_2_, 0.0f);

  const auto c3_0_ = L::celerity()[3][0];
  BOOST_CHECK_EQUAL(c3_0_, 0.0f);
  const auto c3_1_ = L::celerity()[3][1];
  BOOST_CHECK_EQUAL(c3_1_, 0.0f);
  const auto c3_2_ = L::celerity()[3][2];
  BOOST_CHECK_EQUAL(c3_2_, -1.0f);

  const auto c4_0_ = L::celerity()[4][0];
  BOOST_CHECK_EQUAL(c4_0_, -1.0f);
  const auto c4_1_ = L::celerity()[4][1];
  BOOST_CHECK_EQUAL(c4_1_, -1.0f);
  const auto c4_2_ = L::celerity()[4][2];
  BOOST_CHECK_EQUAL(c4_2_, -1.0f);

  const auto c5_0_ = L::celerity()[5][0];
  BOOST_CHECK_EQUAL(c5_0_, -1.0f);
  const auto c5_1_ = L::celerity()[5][1];
  BOOST_CHECK_EQUAL(c5_1_, -1.0f);
  const auto c5_2_ = L::celerity()[5][2];
  BOOST_CHECK_EQUAL(c5_2_, 1.0f);

  const auto c6_0_ = L::celerity()[6][0];
  BOOST_CHECK_EQUAL(c6_0_, -1.0f);
  const auto c6_1_ = L::celerity()[6][1];
  BOOST_CHECK_EQUAL(c6_1_, 1.0f);
  const auto c6_2_ = L::celerity()[6][2];
  BOOST_CHECK_EQUAL(c6_2_, -1.0f);

  const auto c7_0_ = L::celerity()[7][0];
  BOOST_CHECK_EQUAL(c7_0_, -1.0f);
  const auto c7_1_ = L::celerity()[7][1];
  BOOST_CHECK_EQUAL(c7_1_, 1.0f);
  const auto c7_2_ = L::celerity()[7][2];
  BOOST_CHECK_EQUAL(c7_2_, 1.0f);

  const auto c8_0_ = L::celerity()[8][0];
  BOOST_CHECK_EQUAL(c8_0_, 1.0f);
  const auto c8_1_ = L::celerity()[8][1];
  BOOST_CHECK_EQUAL(c8_1_, 0.0f);
  const auto c8_2_ = L::celerity()[8][2];
  BOOST_CHECK_EQUAL(c8_2_, 0.0f);

  const auto c9_0_ = L::celerity()[9][0];
  BOOST_CHECK_EQUAL(c9_0_, 0.0f);
  const auto c9_1_ = L::celerity()[9][1];
  BOOST_CHECK_EQUAL(c9_1_, 1.0f);
  const auto c9_2_ = L::celerity()[9][2];
  BOOST_CHECK_EQUAL(c9_2_, 0.0f);

  const auto c10_0_ = L::celerity()[10][0];
  BOOST_CHECK_EQUAL(c10_0_, 0.0f);
  const auto c10_1_ = L::celerity()[10][1];
  BOOST_CHECK_EQUAL(c10_1_, 0.0f);
  const auto c10_2_ = L::celerity()[10][2];
  BOOST_CHECK_EQUAL(c10_2_, 1.0f);

  const auto c11_0_ = L::celerity()[11][0];
  BOOST_CHECK_EQUAL(c11_0_, 1.0f);
  const auto c11_1_ = L::celerity()[11][1];
  BOOST_CHECK_EQUAL(c11_1_, 1.0f);
  const auto c11_2_ = L::celerity()[11][2];
  BOOST_CHECK_EQUAL(c11_2_, 1.0f);

  const auto c12_0_ = L::celerity()[12][0];
  BOOST_CHECK_EQUAL(c12_0_, 1.0f);
  const auto c12_1_ = L::celerity()[12][1];
  BOOST_CHECK_EQUAL(c12_1_, 1.0f);
  const auto c12_2_ = L::celerity()[12][2];
  BOOST_CHECK_EQUAL(c12_2_, -1.0f);

  const auto c13_0_ = L::celerity()[13][0];
  BOOST_CHECK_EQUAL(c13_0_, 1.0f);
  const auto c13_1_ = L::celerity()[13][1];
  BOOST_CHECK_EQUAL(c13_1_, -1.0f);
  const auto c13_2_ = L::celerity()[13][2];
  BOOST_CHECK_EQUAL(c13_2_, 1.0f);

  const auto c14_0_ = L::celerity()[14][0];
  BOOST_CHECK_EQUAL(c14_0_, 1.0f);
  const auto c14_1_ = L::celerity()[14][1];
  BOOST_CHECK_EQUAL(c14_1_, -1.0f);
  const auto c14_2_ = L::celerity()[14][2];
  BOOST_CHECK_EQUAL(c14_2_, -1.0f);
}

BOOST_AUTO_TEST_CASE(isotropyD3Q19) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D3Q19;
  typedef double valueType;
  typedef ::lbm::Lattice<valueType, latticeType> L;

  valueType sumWeight = 0.0;
  valueType sumWeight_r = 1.0;

  MathVector<valueType, L::dimD> sumWeightCelerity{{0.0}};
  MathVector<valueType, L::dimD> sumWeightCelerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD> sumWeight2Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD> sumWeight2Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD> sumWeight3Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD> sumWeight3Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD> sumWeight4Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD> sumWeight4Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD*L::dimD> sumWeight5Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD*L::dimD> sumWeight5Celerity_r{{0.0}};


  for(int iQ = 0; iQ < L::dimQ; ++iQ) {
    sumWeight += L::weight()[iQ];

    for(int d_a = 0; d_a < L::dimD; ++d_a) {
      int idx_a = d_a;
      sumWeightCelerity[idx_a] += L::weight()[iQ]*L::celerity()[iQ][d_a];

      for(int d_b = 0; d_b < L::dimD; ++d_b) {
        int idx_b = L::dimD * idx_a + d_b;
        sumWeight2Celerity[idx_b] += L::weight()[iQ]*L::celerity()[iQ][d_a]*L::celerity()[iQ][d_b];

        if(d_a == d_b) {
          sumWeight2Celerity_r[idx_b] = L::cs2;
        }

        for(int d_c = 0; d_c < L::dimD; ++d_c) {
          int idx_c = L::dimD * idx_b + d_c;
          sumWeight3Celerity[idx_c] += L::weight()[iQ]*L::celerity()[iQ][d_a]
            *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c];

          for(int d_d = 0; d_d < L::dimD; ++d_d) {
            int idx_d = L::dimD * idx_c + d_d;
            sumWeight4Celerity[idx_d] += L::weight()[iQ]*L::celerity()[iQ][d_a]
              *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c]*L::celerity()[iQ][d_d];

            if((d_a == d_b && d_c == d_d)
               && (d_a == d_c && d_b == d_d)
               && (d_a == d_d && d_b == d_c)) {
              sumWeight4Celerity_r[idx_d] = 3.0*L::cs2*L::cs2;
            }
            else if(((d_a == d_b && d_c == d_d) && (d_a == d_c && d_b == d_d))
                    || ((d_a == d_b && d_c == d_d) && (d_a == d_d && d_b == d_d))
                    || ((d_a == d_d && d_b == d_c) && (d_a == d_d && d_b == d_d))) {
              sumWeight4Celerity_r[idx_d] = 2.0*L::cs2*L::cs2;
            }
            else if((d_a == d_b && d_c == d_d)
                    || (d_a == d_c && d_b == d_d)
                    || (d_a == d_d && d_b == d_c)) {
              sumWeight4Celerity_r[idx_d] = 1.0*L::cs2*L::cs2;
            }

            for(int d_e = 0; d_e < L::dimD; ++d_e) {
              int idx_e = L::dimD * idx_d + d_e;
              sumWeight5Celerity[idx_e] += L::weight()[iQ]*L::celerity()[iQ][d_a]
                *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c]*L::celerity()[iQ][d_d]
                *L::celerity()[iQ][d_e];
            }
          }
        }
      }
    }
  }

  BOOST_TEST(sumWeight == sumWeight_r, tt::tolerance(1e-15));
  for(int i = 0; i < L::dimD; ++i) {
    BOOST_TEST(sumWeightCelerity[i] == sumWeightCelerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight2Celerity[i] == sumWeight2Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight3Celerity[i] == sumWeight3Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight4Celerity[i] == sumWeight4Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight5Celerity[i] == sumWeight5Celerity_r[i], tt::tolerance(1e-15));
  }
}


BOOST_AUTO_TEST_CASE(doubleD3Q19) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D3Q19;
  typedef double valueType;

  typedef ::lbm::Lattice<valueType, latticeType> L;

  const auto dimD_ = L::dimD;
  BOOST_CHECK_EQUAL(dimD_, 3);
  const auto dimQ_ = L::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 19);

  const auto lX_g_ = L::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = L::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, lengthY_g);
  const auto lZ_g_ = L::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, lengthZ_g);

  const auto hX_ = L::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = L::hY;
  BOOST_CHECK_EQUAL(hY_, 1);
  const auto hZ_ = L::hZ;
  BOOST_CHECK_EQUAL(hZ_, 1);

  const auto inv_cs2_ = L::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, (3.0));
  const auto cs2_ = L::cs2;
  BOOST_CHECK_EQUAL(cs2_, (1.0/3.0));

  const auto w0_ = L::weight()[0];
  BOOST_CHECK_EQUAL(w0_, 1.0/3.0);
  const auto w1_ = L::weight()[1];
  BOOST_CHECK_EQUAL(w1_, 1.0/18.0);
  const auto w2_ = L::weight()[2];
  BOOST_CHECK_EQUAL(w2_, 1.0/18.0);
  const auto w3_ = L::weight()[3];
  BOOST_CHECK_EQUAL(w3_, 1.0/18.0);
  const auto w4_ = L::weight()[4];
  BOOST_CHECK_EQUAL(w4_, 1.0/36.0);
  const auto w5_ = L::weight()[5];
  BOOST_CHECK_EQUAL(w5_, 1.0/36.0);
  const auto w6_ = L::weight()[6];
  BOOST_CHECK_EQUAL(w6_, 1.0/36.0);
  const auto w7_ = L::weight()[7];
  BOOST_CHECK_EQUAL(w7_, 1.0/36.0);
  const auto w8_ = L::weight()[8];
  BOOST_CHECK_EQUAL(w8_, 1.0/36.0);
  const auto w9_ = L::weight()[9];
  BOOST_CHECK_EQUAL(w9_, 1.0/36.0);
  const auto w10_ = L::weight()[10];
  BOOST_CHECK_EQUAL(w10_, 1.0/18.0);
  const auto w11_ = L::weight()[11];
  BOOST_CHECK_EQUAL(w11_, 1.0/18.0);
  const auto w12_ = L::weight()[12];
  BOOST_CHECK_EQUAL(w12_, 1.0/18.0);
  const auto w13_ = L::weight()[13];
  BOOST_CHECK_EQUAL(w13_, 1.0/36.0);
  const auto w14_ = L::weight()[14];
  BOOST_CHECK_EQUAL(w14_, 1.0/36.0);
  const auto w15_ = L::weight()[15];
  BOOST_CHECK_EQUAL(w15_, 1.0/36.0);
  const auto w16_ = L::weight()[16];
  BOOST_CHECK_EQUAL(w16_, 1.0/36.0);
  const auto w17_ = L::weight()[17];
  BOOST_CHECK_EQUAL(w17_, 1.0/36.0);
  const auto w18_ = L::weight()[18];
  BOOST_CHECK_EQUAL(w18_, 1.0/36.0);

  const auto c0_ = L::celerity()[0];
  const auto c0_r = MathVector<valueType, dimD_>{{0.0, 0.0, 0.0}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = L::celerity()[1];
  const auto c1_r = MathVector<valueType, dimD_>{{-1.0, 0.0, 0.0}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = L::celerity()[2];
  const auto c2_r = MathVector<valueType, dimD_>{{0.0, -1.0, 0.0}};
  BOOST_CHECK_EQUAL(c2_, c2_r);
  const auto c3_ = L::celerity()[3];
  const auto c3_r = MathVector<valueType, dimD_>{{0.0, 0.0, -1.0}};
  BOOST_CHECK_EQUAL(c3_, c3_r);
  const auto c4_ = L::celerity()[4];
  const auto c4_r = MathVector<valueType, dimD_>{{-1.0, -1.0, 0.0}};
  BOOST_CHECK_EQUAL(c4_, c4_r);
  const auto c5_ = L::celerity()[5];
  const auto c5_r = MathVector<valueType, dimD_>{{-1.0, 1.0, 0.0}};
  BOOST_CHECK_EQUAL(c5_, c5_r);
  const auto c6_ = L::celerity()[6];
  const auto c6_r = MathVector<valueType, dimD_>{{-1.0, 0.0, -1.0}};
  BOOST_CHECK_EQUAL(c6_, c6_r);
  const auto c7_ = L::celerity()[7];
  const auto c7_r = MathVector<valueType, dimD_>{{-1.0, 0.0, 1.0}};
  BOOST_CHECK_EQUAL(c7_, c7_r);
  const auto c8_ = L::celerity()[8];
  const auto c8_r = MathVector<valueType, dimD_>{{0.0, -1.0, -1.0}};
  BOOST_CHECK_EQUAL(c8_, c8_r);
  const auto c9_ = L::celerity()[9];
  const auto c9_r = MathVector<valueType, dimD_>{{0.0, -1.0, 1.0}};
  BOOST_CHECK_EQUAL(c9_, c9_r);
  const auto c10_ = L::celerity()[10];
  const auto c10_r = MathVector<valueType, dimD_>{{1.0, 0.0, 0.0}};
  BOOST_CHECK_EQUAL(c10_, c10_r);
  const auto c11_ = L::celerity()[11];
  const auto c11_r = MathVector<valueType, dimD_>{{0.0, 1.0, 0.0}};
  BOOST_CHECK_EQUAL(c11_, c11_r);
  const auto c12_ = L::celerity()[12];
  const auto c12_r = MathVector<valueType, dimD_>{{0.0, 0.0, 1.0}};
  BOOST_CHECK_EQUAL(c12_, c12_r);
  const auto c13_ = L::celerity()[13];
  const auto c13_r = MathVector<valueType, dimD_>{{1.0, 1.0, 0.0}};
  BOOST_CHECK_EQUAL(c13_, c13_r);
  const auto c14_ = L::celerity()[14];
  const auto c14_r = MathVector<valueType, dimD_>{{1.0, -1.0, 0.0}};
  BOOST_CHECK_EQUAL(c14_, c14_r);
  const auto c15_ = L::celerity()[15];
  const auto c15_r = MathVector<valueType, dimD_>{{1.0, 0.0, 1.0}};
  BOOST_CHECK_EQUAL(c15_, c15_r);
  const auto c16_ = L::celerity()[16];
  const auto c16_r = MathVector<valueType, dimD_>{{1.0, 0.0, -1.0}};
  BOOST_CHECK_EQUAL(c16_, c16_r);
  const auto c17_ = L::celerity()[17];
  const auto c17_r = MathVector<valueType, dimD_>{{0.0, 1.0, 1.0}};
  BOOST_CHECK_EQUAL(c17_, c17_r);
  const auto c18_ = L::celerity()[18];
  const auto c18_r = MathVector<valueType, dimD_>{{0.0, 1.0, -1.0}};
  BOOST_CHECK_EQUAL(c18_, c18_r);

  const auto c0_0_ = L::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, 0.0);
  const auto c0_1_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_1_, 0.0);
  const auto c0_2_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_2_, 0.0);

  const auto c1_0_ = L::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, -1.0);
  const auto c1_1_ = L::celerity()[1][1];
  BOOST_CHECK_EQUAL(c1_1_, 0.0);
  const auto c1_2_ = L::celerity()[1][2];
  BOOST_CHECK_EQUAL(c1_2_, 0.0);

  const auto c2_0_ = L::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, 0.0);
  const auto c2_1_ = L::celerity()[2][1];
  BOOST_CHECK_EQUAL(c2_1_, -1.0);
  const auto c2_2_ = L::celerity()[2][2];
  BOOST_CHECK_EQUAL(c2_2_, 0.0);

  const auto c3_0_ = L::celerity()[3][0];
  BOOST_CHECK_EQUAL(c3_0_, 0.0);
  const auto c3_1_ = L::celerity()[3][1];
  BOOST_CHECK_EQUAL(c3_1_, 0.0);
  const auto c3_2_ = L::celerity()[3][2];
  BOOST_CHECK_EQUAL(c3_2_, -1.0);

  const auto c4_0_ = L::celerity()[4][0];
  BOOST_CHECK_EQUAL(c4_0_, -1.0);
  const auto c4_1_ = L::celerity()[4][1];
  BOOST_CHECK_EQUAL(c4_1_, -1.0);
  const auto c4_2_ = L::celerity()[4][2];
  BOOST_CHECK_EQUAL(c4_2_, 0.0);

  const auto c5_0_ = L::celerity()[5][0];
  BOOST_CHECK_EQUAL(c5_0_, -1.0);
  const auto c5_1_ = L::celerity()[5][1];
  BOOST_CHECK_EQUAL(c5_1_, 1.0);
  const auto c5_2_ = L::celerity()[5][2];
  BOOST_CHECK_EQUAL(c5_2_, 0.0);

  const auto c6_0_ = L::celerity()[6][0];
  BOOST_CHECK_EQUAL(c6_0_, -1.0);
  const auto c6_1_ = L::celerity()[6][1];
  BOOST_CHECK_EQUAL(c6_1_, 0.0);
  const auto c6_2_ = L::celerity()[6][2];
  BOOST_CHECK_EQUAL(c6_2_, -1.0);

  const auto c7_0_ = L::celerity()[7][0];
  BOOST_CHECK_EQUAL(c7_0_, -1.0);
  const auto c7_1_ = L::celerity()[7][1];
  BOOST_CHECK_EQUAL(c7_1_, 0.0);
  const auto c7_2_ = L::celerity()[7][2];
  BOOST_CHECK_EQUAL(c7_2_, 1.0);

  const auto c8_0_ = L::celerity()[8][0];
  BOOST_CHECK_EQUAL(c8_0_, 0.0);
  const auto c8_1_ = L::celerity()[8][1];
  BOOST_CHECK_EQUAL(c8_1_, -1.0);
  const auto c8_2_ = L::celerity()[8][2];
  BOOST_CHECK_EQUAL(c8_2_, -1.0);

  const auto c9_0_ = L::celerity()[9][0];
  BOOST_CHECK_EQUAL(c9_0_, 0.0);
  const auto c9_1_ = L::celerity()[9][1];
  BOOST_CHECK_EQUAL(c9_1_, -1.0);
  const auto c9_2_ = L::celerity()[9][2];
  BOOST_CHECK_EQUAL(c9_2_, 1.0);

  const auto c10_0_ = L::celerity()[10][0];
  BOOST_CHECK_EQUAL(c10_0_, 1.0);
  const auto c10_1_ = L::celerity()[10][1];
  BOOST_CHECK_EQUAL(c10_1_, 0.0);
  const auto c10_2_ = L::celerity()[10][2];
  BOOST_CHECK_EQUAL(c10_2_, 0.0);

  const auto c11_0_ = L::celerity()[11][0];
  BOOST_CHECK_EQUAL(c11_0_, 0.0);
  const auto c11_1_ = L::celerity()[11][1];
  BOOST_CHECK_EQUAL(c11_1_, 1.0);
  const auto c11_2_ = L::celerity()[11][2];
  BOOST_CHECK_EQUAL(c11_2_, 0.0);

  const auto c12_0_ = L::celerity()[12][0];
  BOOST_CHECK_EQUAL(c12_0_, 0.0);
  const auto c12_1_ = L::celerity()[12][1];
  BOOST_CHECK_EQUAL(c12_1_, 0.0);
  const auto c12_2_ = L::celerity()[12][2];
  BOOST_CHECK_EQUAL(c12_2_, 1.0);

  const auto c13_0_ = L::celerity()[13][0];
  BOOST_CHECK_EQUAL(c13_0_, 1.0);
  const auto c13_1_ = L::celerity()[13][1];
  BOOST_CHECK_EQUAL(c13_1_, 1.0);
  const auto c13_2_ = L::celerity()[13][2];
  BOOST_CHECK_EQUAL(c13_2_, 0.0);

  const auto c14_0_ = L::celerity()[14][0];
  BOOST_CHECK_EQUAL(c14_0_, 1.0);
  const auto c14_1_ = L::celerity()[14][1];
  BOOST_CHECK_EQUAL(c14_1_, -1.0);
  const auto c14_2_ = L::celerity()[14][2];
  BOOST_CHECK_EQUAL(c14_2_, 0.0);

  const auto c15_0_ = L::celerity()[15][0];
  BOOST_CHECK_EQUAL(c15_0_, 1.0);
  const auto c15_1_ = L::celerity()[15][1];
  BOOST_CHECK_EQUAL(c15_1_, 0.0);
  const auto c15_2_ = L::celerity()[15][2];
  BOOST_CHECK_EQUAL(c15_2_, 1.0);

  const auto c16_0_ = L::celerity()[16][0];
  BOOST_CHECK_EQUAL(c16_0_, 1.0);
  const auto c16_1_ = L::celerity()[16][1];
  BOOST_CHECK_EQUAL(c16_1_, 0.0);
  const auto c16_2_ = L::celerity()[16][2];
  BOOST_CHECK_EQUAL(c16_2_, -1.0);

  const auto c17_0_ = L::celerity()[17][0];
  BOOST_CHECK_EQUAL(c17_0_, 0.0);
  const auto c17_1_ = L::celerity()[17][1];
  BOOST_CHECK_EQUAL(c17_1_, 1.0);
  const auto c17_2_ = L::celerity()[17][2];
  BOOST_CHECK_EQUAL(c17_2_, 1.0);

  const auto c18_0_ = L::celerity()[18][0];
  BOOST_CHECK_EQUAL(c18_0_, 0.0);
  const auto c18_1_ = L::celerity()[18][1];
  BOOST_CHECK_EQUAL(c18_1_, 1.0);
  const auto c18_2_ = L::celerity()[18][2];
  BOOST_CHECK_EQUAL(c18_2_, -1.0);

}

BOOST_AUTO_TEST_CASE(floatD3Q19) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D3Q19;
  typedef float valueType;

  typedef ::lbm::Lattice<valueType, latticeType> L;

  const auto dimD_ = L::dimD;
  BOOST_CHECK_EQUAL(dimD_, 3);
  const auto dimQ_ = L::dimQ;
  BOOST_CHECK_EQUAL(dimQ_, 19);

  const auto lX_g_ = L::lX_g;
  BOOST_CHECK_EQUAL(lX_g_, lengthX_g);
  const auto lY_g_ = L::lY_g;
  BOOST_CHECK_EQUAL(lY_g_, lengthY_g);
  const auto lZ_g_ = L::lZ_g;
  BOOST_CHECK_EQUAL(lZ_g_, lengthZ_g);

  const auto hX_ = L::hX;
  BOOST_CHECK_EQUAL(hX_, 1);
  const auto hY_ = L::hY;
  BOOST_CHECK_EQUAL(hY_, 1);
  const auto hZ_ = L::hZ;
  BOOST_CHECK_EQUAL(hZ_, 1);

  const auto inv_cs2_ = L::inv_cs2;
  BOOST_CHECK_EQUAL(inv_cs2_, (3.0f));
  const auto cs2_ = L::cs2;
  BOOST_CHECK_EQUAL(cs2_, (1.0f/3.0f));

  const auto w0_ = L::weight()[0];
  BOOST_CHECK_EQUAL(w0_, 1.0f/3.0f);
  const auto w1_ = L::weight()[1];
  BOOST_CHECK_EQUAL(w1_, 1.0f/18.0f);
  const auto w2_ = L::weight()[2];
  BOOST_CHECK_EQUAL(w2_, 1.0f/18.0f);
  const auto w3_ = L::weight()[3];
  BOOST_CHECK_EQUAL(w3_, 1.0f/18.0f);
  const auto w4_ = L::weight()[4];
  BOOST_CHECK_EQUAL(w4_, 1.0f/36.0f);
  const auto w5_ = L::weight()[5];
  BOOST_CHECK_EQUAL(w5_, 1.0f/36.0f);
  const auto w6_ = L::weight()[6];
  BOOST_CHECK_EQUAL(w6_, 1.0f/36.0f);
  const auto w7_ = L::weight()[7];
  BOOST_CHECK_EQUAL(w7_, 1.0f/36.0f);
  const auto w8_ = L::weight()[8];
  BOOST_CHECK_EQUAL(w8_, 1.0f/36.0f);
  const auto w9_ = L::weight()[9];
  BOOST_CHECK_EQUAL(w9_, 1.0f/36.0f);
  const auto w10_ = L::weight()[10];
  BOOST_CHECK_EQUAL(w10_, 1.0f/18.0f);
  const auto w11_ = L::weight()[11];
  BOOST_CHECK_EQUAL(w11_, 1.0f/18.0f);
  const auto w12_ = L::weight()[12];
  BOOST_CHECK_EQUAL(w12_, 1.0f/18.0f);
  const auto w13_ = L::weight()[13];
  BOOST_CHECK_EQUAL(w13_, 1.0f/36.0f);
  const auto w14_ = L::weight()[14];
  BOOST_CHECK_EQUAL(w14_, 1.0f/36.0f);
  const auto w15_ = L::weight()[15];
  BOOST_CHECK_EQUAL(w15_, 1.0f/36.0f);
  const auto w16_ = L::weight()[16];
  BOOST_CHECK_EQUAL(w16_, 1.0f/36.0f);
  const auto w17_ = L::weight()[17];
  BOOST_CHECK_EQUAL(w17_, 1.0f/36.0f);
  const auto w18_ = L::weight()[18];
  BOOST_CHECK_EQUAL(w18_, 1.0f/36.0f);

  const auto c0_ = L::celerity()[0];
  const auto c0_r = MathVector<valueType, dimD_>{{0.0f, 0.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c0_, c0_r);
  const auto c1_ = L::celerity()[1];
  const auto c1_r = MathVector<valueType, dimD_>{{-1.0f, 0.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c1_, c1_r);
  const auto c2_ = L::celerity()[2];
  const auto c2_r = MathVector<valueType, dimD_>{{0.0f, -1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c2_, c2_r);
  const auto c3_ = L::celerity()[3];
  const auto c3_r = MathVector<valueType, dimD_>{{0.0f, 0.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c3_, c3_r);
  const auto c4_ = L::celerity()[4];
  const auto c4_r = MathVector<valueType, dimD_>{{-1.0f, -1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c4_, c4_r);
  const auto c5_ = L::celerity()[5];
  const auto c5_r = MathVector<valueType, dimD_>{{-1.0f, 1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c5_, c5_r);
  const auto c6_ = L::celerity()[6];
  const auto c6_r = MathVector<valueType, dimD_>{{-1.0f, 0.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c6_, c6_r);
  const auto c7_ = L::celerity()[7];
  const auto c7_r = MathVector<valueType, dimD_>{{-1.0f, 0.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c7_, c7_r);
  const auto c8_ = L::celerity()[8];
  const auto c8_r = MathVector<valueType, dimD_>{{0.0f, -1.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c8_, c8_r);
  const auto c9_ = L::celerity()[9];
  const auto c9_r = MathVector<valueType, dimD_>{{0.0f, -1.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c9_, c9_r);
  const auto c10_ = L::celerity()[10];
  const auto c10_r = MathVector<valueType, dimD_>{{1.0f, 0.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c10_, c10_r);
  const auto c11_ = L::celerity()[11];
  const auto c11_r = MathVector<valueType, dimD_>{{0.0f, 1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c11_, c11_r);
  const auto c12_ = L::celerity()[12];
  const auto c12_r = MathVector<valueType, dimD_>{{0.0f, 0.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c12_, c12_r);
  const auto c13_ = L::celerity()[13];
  const auto c13_r = MathVector<valueType, dimD_>{{1.0f, 1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c13_, c13_r);
  const auto c14_ = L::celerity()[14];
  const auto c14_r = MathVector<valueType, dimD_>{{1.0f, -1.0f, 0.0f}};
  BOOST_CHECK_EQUAL(c14_, c14_r);
  const auto c15_ = L::celerity()[15];
  const auto c15_r = MathVector<valueType, dimD_>{{1.0f, 0.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c15_, c15_r);
  const auto c16_ = L::celerity()[16];
  const auto c16_r = MathVector<valueType, dimD_>{{1.0f, 0.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c16_, c16_r);
  const auto c17_ = L::celerity()[17];
  const auto c17_r = MathVector<valueType, dimD_>{{0.0f, 1.0f, 1.0f}};
  BOOST_CHECK_EQUAL(c17_, c17_r);
  const auto c18_ = L::celerity()[18];
  const auto c18_r = MathVector<valueType, dimD_>{{0.0f, 1.0f, -1.0f}};
  BOOST_CHECK_EQUAL(c18_, c18_r);

  const auto c0_0_ = L::celerity()[0][0];
  BOOST_CHECK_EQUAL(c0_0_, 0.0f);
  const auto c0_1_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_1_, 0.0f);
  const auto c0_2_ = L::celerity()[0][1];
  BOOST_CHECK_EQUAL(c0_2_, 0.0f);

  const auto c1_0_ = L::celerity()[1][0];
  BOOST_CHECK_EQUAL(c1_0_, -1.0f);
  const auto c1_1_ = L::celerity()[1][1];
  BOOST_CHECK_EQUAL(c1_1_, 0.0f);
  const auto c1_2_ = L::celerity()[1][2];
  BOOST_CHECK_EQUAL(c1_2_, 0.0f);

  const auto c2_0_ = L::celerity()[2][0];
  BOOST_CHECK_EQUAL(c2_0_, 0.0f);
  const auto c2_1_ = L::celerity()[2][1];
  BOOST_CHECK_EQUAL(c2_1_, -1.0f);
  const auto c2_2_ = L::celerity()[2][2];
  BOOST_CHECK_EQUAL(c2_2_, 0.0f);

  const auto c3_0_ = L::celerity()[3][0];
  BOOST_CHECK_EQUAL(c3_0_, 0.0f);
  const auto c3_1_ = L::celerity()[3][1];
  BOOST_CHECK_EQUAL(c3_1_, 0.0f);
  const auto c3_2_ = L::celerity()[3][2];
  BOOST_CHECK_EQUAL(c3_2_, -1.0f);

  const auto c4_0_ = L::celerity()[4][0];
  BOOST_CHECK_EQUAL(c4_0_, -1.0f);
  const auto c4_1_ = L::celerity()[4][1];
  BOOST_CHECK_EQUAL(c4_1_, -1.0f);
  const auto c4_2_ = L::celerity()[4][2];
  BOOST_CHECK_EQUAL(c4_2_, 0.0f);

  const auto c5_0_ = L::celerity()[5][0];
  BOOST_CHECK_EQUAL(c5_0_, -1.0f);
  const auto c5_1_ = L::celerity()[5][1];
  BOOST_CHECK_EQUAL(c5_1_, 1.0f);
  const auto c5_2_ = L::celerity()[5][2];
  BOOST_CHECK_EQUAL(c5_2_, 0.0f);

  const auto c6_0_ = L::celerity()[6][0];
  BOOST_CHECK_EQUAL(c6_0_, -1.0f);
  const auto c6_1_ = L::celerity()[6][1];
  BOOST_CHECK_EQUAL(c6_1_, 0.0f);
  const auto c6_2_ = L::celerity()[6][2];
  BOOST_CHECK_EQUAL(c6_2_, -1.0f);

  const auto c7_0_ = L::celerity()[7][0];
  BOOST_CHECK_EQUAL(c7_0_, -1.0f);
  const auto c7_1_ = L::celerity()[7][1];
  BOOST_CHECK_EQUAL(c7_1_, 0.0f);
  const auto c7_2_ = L::celerity()[7][2];
  BOOST_CHECK_EQUAL(c7_2_, 1.0f);

  const auto c8_0_ = L::celerity()[8][0];
  BOOST_CHECK_EQUAL(c8_0_, 0.0f);
  const auto c8_1_ = L::celerity()[8][1];
  BOOST_CHECK_EQUAL(c8_1_, -1.0f);
  const auto c8_2_ = L::celerity()[8][2];
  BOOST_CHECK_EQUAL(c8_2_, -1.0f);

  const auto c9_0_ = L::celerity()[9][0];
  BOOST_CHECK_EQUAL(c9_0_, 0.0f);
  const auto c9_1_ = L::celerity()[9][1];
  BOOST_CHECK_EQUAL(c9_1_, -1.0f);
  const auto c9_2_ = L::celerity()[9][2];
  BOOST_CHECK_EQUAL(c9_2_, 1.0f);

  const auto c10_0_ = L::celerity()[10][0];
  BOOST_CHECK_EQUAL(c10_0_, 1.0f);
  const auto c10_1_ = L::celerity()[10][1];
  BOOST_CHECK_EQUAL(c10_1_, 0.0f);
  const auto c10_2_ = L::celerity()[10][2];
  BOOST_CHECK_EQUAL(c10_2_, 0.0f);

  const auto c11_0_ = L::celerity()[11][0];
  BOOST_CHECK_EQUAL(c11_0_, 0.0f);
  const auto c11_1_ = L::celerity()[11][1];
  BOOST_CHECK_EQUAL(c11_1_, 1.0f);
  const auto c11_2_ = L::celerity()[11][2];
  BOOST_CHECK_EQUAL(c11_2_, 0.0f);

  const auto c12_0_ = L::celerity()[12][0];
  BOOST_CHECK_EQUAL(c12_0_, 0.0f);
  const auto c12_1_ = L::celerity()[12][1];
  BOOST_CHECK_EQUAL(c12_1_, 0.0f);
  const auto c12_2_ = L::celerity()[12][2];
  BOOST_CHECK_EQUAL(c12_2_, 1.0f);

  const auto c13_0_ = L::celerity()[13][0];
  BOOST_CHECK_EQUAL(c13_0_, 1.0f);
  const auto c13_1_ = L::celerity()[13][1];
  BOOST_CHECK_EQUAL(c13_1_, 1.0f);
  const auto c13_2_ = L::celerity()[13][2];
  BOOST_CHECK_EQUAL(c13_2_, 0.0f);

  const auto c14_0_ = L::celerity()[14][0];
  BOOST_CHECK_EQUAL(c14_0_, 1.0f);
  const auto c14_1_ = L::celerity()[14][1];
  BOOST_CHECK_EQUAL(c14_1_, -1.0f);
  const auto c14_2_ = L::celerity()[14][2];
  BOOST_CHECK_EQUAL(c14_2_, 0.0f);

  const auto c15_0_ = L::celerity()[15][0];
  BOOST_CHECK_EQUAL(c15_0_, 1.0f);
  const auto c15_1_ = L::celerity()[15][1];
  BOOST_CHECK_EQUAL(c15_1_, 0.0f);
  const auto c15_2_ = L::celerity()[15][2];
  BOOST_CHECK_EQUAL(c15_2_, 1.0f);

  const auto c16_0_ = L::celerity()[16][0];
  BOOST_CHECK_EQUAL(c16_0_, 1.0f);
  const auto c16_1_ = L::celerity()[16][1];
  BOOST_CHECK_EQUAL(c16_1_, 0.0f);
  const auto c16_2_ = L::celerity()[16][2];
  BOOST_CHECK_EQUAL(c16_2_, -1.0f);

  const auto c17_0_ = L::celerity()[17][0];
  BOOST_CHECK_EQUAL(c17_0_, 0.0f);
  const auto c17_1_ = L::celerity()[17][1];
  BOOST_CHECK_EQUAL(c17_1_, 1.0f);
  const auto c17_2_ = L::celerity()[17][2];
  BOOST_CHECK_EQUAL(c17_2_, 1.0f);

  const auto c18_0_ = L::celerity()[18][0];
  BOOST_CHECK_EQUAL(c18_0_, 0.0f);
  const auto c18_1_ = L::celerity()[18][1];
  BOOST_CHECK_EQUAL(c18_1_, 1.0f);
  const auto c18_2_ = L::celerity()[18][2];
  BOOST_CHECK_EQUAL(c18_2_, -1.0f);
}

BOOST_AUTO_TEST_CASE(isotropyD3Q27) {
  constexpr ::lbm::LatticeType latticeType = ::lbm::LatticeType::D3Q27;
  typedef double valueType;
  typedef ::lbm::Lattice<valueType, latticeType> L;

  valueType sumWeight = 0.0;
  valueType sumWeight_r = 1.0;

  MathVector<valueType, L::dimD> sumWeightCelerity{{0.0}};
  MathVector<valueType, L::dimD> sumWeightCelerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD> sumWeight2Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD> sumWeight2Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD> sumWeight3Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD> sumWeight3Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD> sumWeight4Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD> sumWeight4Celerity_r{{0.0}};

  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD*L::dimD> sumWeight5Celerity{{0.0}};
  MathVector<valueType, L::dimD*L::dimD*L::dimD*L::dimD*L::dimD> sumWeight5Celerity_r{{0.0}};


  for(int iQ = 0; iQ < L::dimQ; ++iQ) {
    sumWeight += L::weight()[iQ];

    for(int d_a = 0; d_a < L::dimD; ++d_a) {
      int idx_a = d_a;
      sumWeightCelerity[idx_a] += L::weight()[iQ]*L::celerity()[iQ][d_a];

      for(int d_b = 0; d_b < L::dimD; ++d_b) {
        int idx_b = L::dimD * idx_a + d_b;
        sumWeight2Celerity[idx_b] += L::weight()[iQ]*L::celerity()[iQ][d_a]*L::celerity()[iQ][d_b];

        if(d_a == d_b) {
          sumWeight2Celerity_r[idx_b] = L::cs2;
        }

        for(int d_c = 0; d_c < L::dimD; ++d_c) {
          int idx_c = L::dimD * idx_b + d_c;
          sumWeight3Celerity[idx_c] += L::weight()[iQ]*L::celerity()[iQ][d_a]
            *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c];

          for(int d_d = 0; d_d < L::dimD; ++d_d) {
            int idx_d = L::dimD * idx_c + d_d;
            sumWeight4Celerity[idx_d] += L::weight()[iQ]*L::celerity()[iQ][d_a]
              *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c]*L::celerity()[iQ][d_d];

            if((d_a == d_b && d_c == d_d)
               && (d_a == d_c && d_b == d_d)
               && (d_a == d_d && d_b == d_c)) {
              sumWeight4Celerity_r[idx_d] = 3.0*L::cs2*L::cs2;
            }
            else if(((d_a == d_b && d_c == d_d) && (d_a == d_c && d_b == d_d))
                    || ((d_a == d_b && d_c == d_d) && (d_a == d_d && d_b == d_d))
                    || ((d_a == d_d && d_b == d_c) && (d_a == d_d && d_b == d_d))) {
              sumWeight4Celerity_r[idx_d] = 2.0*L::cs2*L::cs2;
            }
            else if((d_a == d_b && d_c == d_d)
                    || (d_a == d_c && d_b == d_d)
                    || (d_a == d_d && d_b == d_c)) {
              sumWeight4Celerity_r[idx_d] = 1.0*L::cs2*L::cs2;
            }

            for(int d_e = 0; d_e < L::dimD; ++d_e) {
              int idx_e = L::dimD * idx_d + d_e;
              sumWeight5Celerity[idx_e] += L::weight()[iQ]*L::celerity()[iQ][d_a]
                *L::celerity()[iQ][d_b]*L::celerity()[iQ][d_c]*L::celerity()[iQ][d_d]
                *L::celerity()[iQ][d_e];
            }
          }
        }
      }
    }
  }

  BOOST_TEST(sumWeight == sumWeight_r, tt::tolerance(1e-15));
  for(int i = 0; i < L::dimD; ++i) {
    BOOST_TEST(sumWeightCelerity[i] == sumWeightCelerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight2Celerity[i] == sumWeight2Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight3Celerity[i] == sumWeight3Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight4Celerity[i] == sumWeight4Celerity_r[i], tt::tolerance(1e-15));
  }
  for(int i = 0; i < L::dimD*L::dimD*L::dimD*L::dimD*L::dimD; ++i) {
    BOOST_TEST(sumWeight5Celerity[i] == sumWeight5Celerity_r[i], tt::tolerance(1e-15));
  }
}

// TODO: Add doubleD3Q27 and floatD3Q27 tests.

BOOST_AUTO_TEST_SUITE_END()
