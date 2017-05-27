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
#include "Writer.h"
#include "MathVector.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestWriter)

BOOST_AUTO_TEST_CASE(TestVTR) {
  constexpr WriterType writerType = WriterType::Constant;
  constexpr MathVector<valueType, 3> amplitude{{(valueType) 1}};
  constexpr MathVector<int, 3> iP{{1}};

  Writer<valueType, writerType> writer(amplitude);

  BOOST_TEST(writer.write(iP)[d::X] == (valueType) 1,
             tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_SUITE_END()
