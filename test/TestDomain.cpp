#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
namespace tt = boost::test_tools;

#include <iostream>

#define NPROCS 2
#include "Lattice.h"
typedef double valueType;
typedef lbm::Lattice<valueType, lbm::LatticeType::D1Q3> L;
#include "Domain.h"
#include "MathVector.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestDomain)

BOOST_AUTO_TEST_CASE(TestLocal) {
  constexpr DomainType domainType = DomainType::Local;
  Domain<L::latticeType, domainType,
         PartitionningType::Generic, MemoryLayout::Generic> domain;

  BOOST_TEST(domain.volume() == (unsigned int) 16/2*1*1,
             tt::tolerance(1e-15));

  const MathVector<unsigned int, 3> iP_1{{1, 0, 0}};
  domain.updateIndex(iP_1);

  BOOST_TEST(domain.indexLocal == (unsigned int) 1,
             tt::tolerance(1e-15));

  BOOST_TEST(domain.getIndex(2) == (unsigned int) 16+1,
             tt::tolerance(1e-15));

  const MathVector<unsigned int, 3> iP_3{{3, 0, 0}};
  domain.updateIndex(iP_3);

  BOOST_TEST(domain.indexLocal == (unsigned int) 3,
             tt::tolerance(1e-15));
  BOOST_TEST(domain.getIndex(4) == (unsigned int) 32+3,
             tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_CASE(TestGlobal) {
  constexpr DomainType domainType = DomainType::Global;
  Domain<L::latticeType, domainType,
         PartitionningType::Generic, MemoryLayout::Generic> domain;

  BOOST_TEST(domain.volume() == (unsigned int) 16*1*1,
             tt::tolerance(1e-15));

  const MathVector<unsigned int, 3> iP_1{{1, 0, 0}};
  domain.updateIndex(iP_1);

  BOOST_TEST(domain.indexGlobal == (unsigned int) 1,
             tt::tolerance(1e-15));

  BOOST_TEST(domain.getIndex(2) == (unsigned int) 32+1,
             tt::tolerance(1e-15));

  const MathVector<unsigned int, 3> iP_3{{3, 0, 0}};
  domain.updateIndex(iP_3);

  BOOST_TEST(domain.indexGlobal == (unsigned int) 3,
             tt::tolerance(1e-15));
  BOOST_TEST(domain.getIndex(4) == (unsigned int) 64+3,
             tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_CASE(TestHaloAoS) {
  constexpr DomainType domainType = DomainType::Halo;
  Domain<L::latticeType, domainType,
         PartitionningType::Generic, MemoryLayout::AoS> domain;

  BOOST_TEST(domain.volume() == (unsigned int) (8+2)*1*1,
             tt::tolerance(1e-15));

  const MathVector<unsigned int, 3> iP_1{{1, 0, 0}};
  domain.updateIndex(iP_1);

  BOOST_TEST(domain.indexHalo == (unsigned int) 1,
             tt::tolerance(1e-15));
  BOOST_TEST(domain.indexLocal == (unsigned int) 0,
             tt::tolerance(1e-15));
  BOOST_TEST(domain.getIndex(2) == (unsigned int) 3+2,
             tt::tolerance(1e-15));

  const MathVector<unsigned int, 3> iP_3{{3, 0, 0}};
  domain.updateIndex(iP_3);

  BOOST_TEST(domain.indexHalo == (unsigned int) 3,
             tt::tolerance(1e-15));
    BOOST_TEST(domain.indexLocal == (unsigned int) 2,
             tt::tolerance(1e-15));
  BOOST_TEST(domain.getIndex(4) == (unsigned int) 9+4,
             tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_CASE(TestHaloSoA) {
  constexpr DomainType domainType = DomainType::Halo;
  Domain<L::latticeType, domainType,
         PartitionningType::Generic, MemoryLayout::SoA> domain;

  BOOST_TEST(domain.volume() == (unsigned int) (8+2)*1*1,
             tt::tolerance(1e-15));

  const MathVector<unsigned int, 3> iP_1{{1, 0, 0}};
  domain.updateIndex(iP_1);

  BOOST_TEST(domain.indexHalo == (unsigned int) 1,
             tt::tolerance(1e-15));
  BOOST_TEST(domain.indexLocal == (unsigned int) 0,
             tt::tolerance(1e-15));
  BOOST_TEST(domain.getIndex(2) == (unsigned int) 2*10+1,
             tt::tolerance(1e-15));

  const MathVector<unsigned int, 3> iP_3{{3, 0, 0}};
  domain.updateIndex(iP_3);

  BOOST_TEST(domain.indexHalo == (unsigned int) 3,
             tt::tolerance(1e-15));
    BOOST_TEST(domain.indexLocal == (unsigned int) 2,
             tt::tolerance(1e-15));
  BOOST_TEST(domain.getIndex(4) == (unsigned int) 4*10+3,
             tt::tolerance(1e-15));
}


BOOST_AUTO_TEST_SUITE_END()
