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
#include "commons.h"
#include "structure.h"
#include "forcingScheme.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestForcingScheme)

BOOST_AUTO_TEST_CASE(TestGuo) {
  auto forcingScheme = std::shared_ptr<ForcingScheme<valueType>>(new Guo<valueType>());

  vector<MathVector<valueType, L::dimD>> eqVelocityForcing_r(L::lX_g *L::lY_g*L::lZ_g);
  vector<MathVector<valueType, L::dimQ>> collisionForcing1_r(sizeX_l()*sizeY_l()*sizeZ_l());
  vector<MathVector<valueType, L::dimQ>> collisionForcing2_r(sizeX_l()*sizeY_l()*sizeZ_l());
  vector<MathVector<valueType, L::dimD>> hydroVelocityForcing_r(L::lX_g *L::lY_g*L::lZ_g);

  for(int idx = 0; idx < L::lX_g *L::lY_g*L::lZ_g; ++idx) {
    valueType density = 0;
    MathVector<valueType, L::dimD> velocity{(valueType) 0};
    for(int iQ = 0; iQ < L::dimQ; ++iQ) {
      forcingScheme->force = MathVector<valueType, L::dimD>{(valueType) 0};
      density = (valueType)(idx+1)/(valueType) (L::lX_g*L::lY_g*L::lZ_g);
      velocity = L::celerity()[iQ];
      collisionForcing1_r[idx][iQ] = forcingScheme->getCollisionForcing(iQ, density,
                                                                        velocity,
                                                                        velocity.norm2());
      forcingScheme->force = MathVector<valueType, L::dimD>{(valueType)(L::weight()[iQ]*(idx+1))/(valueType) (L::lX_g*L::lY_g*L::lZ_g)};
      density = (valueType)(idx+1)/(valueType) (L::lX_g*L::lY_g*L::lZ_g);
      velocity = MathVector<valueType, L::dimD>{(valueType) 0};
      collisionForcing2_r[idx][iQ] = forcingScheme->getCollisionForcing(iQ, density,
                                                                        velocity,
                                                                        velocity.norm2());
    }
    eqVelocityForcing_r[idx] = forcingScheme->getEqVelocityForcing(density);
    hydroVelocityForcing_r[idx] = forcingScheme->getHydroVelocityForcing(density);

  }

  for(int idx = 0; idx < L::lX_g *L::lY_g*L::lZ_g; ++idx) {
    for(int d = 0; d < L::dimD; ++d) {
      BOOST_TEST(hydroVelocityForcing_r[idx][d] == (valueType) (0.5*L::weight()[L::dimQ-1]),
                 tt::tolerance(1e-15));
      BOOST_TEST(eqVelocityForcing_r[idx][d] == (valueType) (0.5*L::weight()[L::dimQ-1]),
                 tt::tolerance(1e-15));
    }
    for(int iQ = 0; iQ < L::dimQ; ++iQ) {
      BOOST_TEST(collisionForcing1_r[idx][iQ] == (valueType) 0,
                 tt::tolerance(1e-15));
      BOOST_TEST(collisionForcing2_r[idx][iQ] == (1.0-beta)*L::weight()[iQ]*L::inv_cs2
                 *(valueType)(L::weight()[iQ]*(idx+1))/(valueType) (L::lX_g*L::lY_g*L::lZ_g)
                 *L::celerity()[iQ].dot(MathVector<valueType, L::dimD>{(valueType) 1}),
                 tt::tolerance(1e-15));
    }
  }
}

BOOST_AUTO_TEST_CASE(TestShanChen) {
  auto forcingScheme = std::shared_ptr<ForcingScheme<valueType>>(new ShanChen<valueType>());

  vector<MathVector<valueType, L::dimD>> eqVelocityForcing_r(L::lX_g *L::lY_g*L::lZ_g);
  vector<MathVector<valueType, L::dimQ>> collisionForcing1_r(sizeX_l()*sizeY_l()*sizeZ_l());
  vector<MathVector<valueType, L::dimQ>> collisionForcing2_r(sizeX_l()*sizeY_l()*sizeZ_l());
  vector<MathVector<valueType, L::dimD>> hydroVelocityForcing_r(L::lX_g *L::lY_g*L::lZ_g);

  for(int idx = 0; idx < L::lX_g *L::lY_g*L::lZ_g; ++idx) {
    valueType density = 0;
    MathVector<valueType, L::dimD> velocity{(valueType) 0};
    for(int iQ = 0; iQ < L::dimQ; ++iQ) {
      forcingScheme->force = MathVector<valueType, L::dimD>{(valueType) 0};
      density = (valueType)(idx+1)/(valueType) (L::lX_g*L::lY_g*L::lZ_g);
      velocity = L::celerity()[iQ];
      collisionForcing1_r[idx][iQ] = forcingScheme->getCollisionForcing(iQ, density,
                                                                        velocity,
                                                                        velocity.norm2());
      forcingScheme->force = MathVector<valueType, L::dimD>{(valueType)(L::weight()[iQ]*(idx+1))/(valueType) (L::lX_g*L::lY_g*L::lZ_g)};
      density = (valueType)(idx+1)/(valueType) (L::lX_g*L::lY_g*L::lZ_g);
      velocity = MathVector<valueType, L::dimD>{(valueType) 0};
      collisionForcing2_r[idx][iQ] = forcingScheme->getCollisionForcing(iQ, density,
                                                                        velocity,
                                                                        velocity.norm2());
    }
    eqVelocityForcing_r[idx] = forcingScheme->getEqVelocityForcing(density);
    hydroVelocityForcing_r[idx] = forcingScheme->getHydroVelocityForcing(density);

  }

  for(int idx = 0; idx < L::lX_g *L::lY_g*L::lZ_g; ++idx) {
    for(int d = 0; d < L::dimD; ++d) {
      BOOST_TEST(hydroVelocityForcing_r[idx][d] == (valueType) (0.5*L::weight()[L::dimQ-1]),
                 tt::tolerance(1e-15));
      BOOST_TEST(eqVelocityForcing_r[idx][d] == (valueType) (tau*L::weight()[L::dimQ-1]),
                 tt::tolerance(1e-15));
    }
    for(int iQ = 0; iQ < L::dimQ; ++iQ) {
      BOOST_TEST(collisionForcing1_r[idx][iQ] == (valueType) 0,
                 tt::tolerance(1e-15));
      BOOST_TEST(collisionForcing2_r[idx][iQ] == (valueType) 0,
                 tt::tolerance(1e-15));
    }
  }
}


BOOST_AUTO_TEST_SUITE_END()
