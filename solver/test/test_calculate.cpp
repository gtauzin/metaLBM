#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
namespace tt = boost::test_tools;

#include <iostream>
#include <vector>
template<class T, std::size_t Alignment = 1>
#ifdef ALIGN_MEMORY
using vector = std::vector<T,
                           boost::alignment::aligned_allocator<T, Alignment> >;
#else
using vector = std::vector<T>;
#endif

#define NPROCS 1
#include "input.h"
#include "lattice.h"
#include "commons.h"
#include "structure.h"
#include "calculate.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(Calculate)

BOOST_AUTO_TEST_CASE(ComputeDensity) {
  vector<valueType> fPop(L::dimQ
                         *(L::lX_l+2*L::hX)
                         *(L::lY_l+2*L::hY)
                         *(L::lZ_l+2*L::hZ));
  vector<valueType> density_r(L::lX_g *L::lY_g*L::lZ_g);

  for(int idx = 0; idx < L::lX_g *L::lY_g*L::lZ_g; ++idx) {
    for(int iQ = 0; iQ < L::dimQ; ++iQ) {
      fPop[idxPop(idx, iQ)] = (valueType) (idx + iQ);
    }
    density_r[idx] = computeDensity<valueType>(fPop.data(), idx);
  }

  for(int idx = 0; idx <L::lX_g *L::lY_g*L::lZ_g; ++idx) {
    BOOST_TEST(density_r[idx] == (valueType)(L::dimQ*idx + (L::dimQ-1)*L::dimQ/2),
               tt::tolerance(1e-15));
  }
}

BOOST_AUTO_TEST_CASE(ComputeVelocity) {
  vector<valueType> fPop(L::dimQ*sizeX_l()*sizeY_l()*sizeZ_l(), (valueType) 0);
  vector<MathVector<valueType, L::dimD>> velocity_r(L::lX_g *L::lY_g*L::lZ_g,
                                                    MathVector<valueType, L::dimD>{0});

  for(int idx = 0; idx < L::lX_g *L::lY_g*L::lZ_g; ++idx) {
    for(int iQ = 0; iQ < L::dimQ; ++iQ) {
      fPop[idxPop(idx, iQ)] = (valueType) L::weight()[iQ];
    }
    velocity_r[idx] = computeVelocity<valueType>(fPop.data(), idx,
                                                 (valueType) idx+1);
  }

  for(int idx = 0; idx < L::lX_g *L::lY_g*L::lZ_g; ++idx) {
    for(int d = 0; d < L::dimD; ++d) {
      BOOST_TEST(velocity_r[idx][d] == (valueType) 0,
                 tt::tolerance(1e-15));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
