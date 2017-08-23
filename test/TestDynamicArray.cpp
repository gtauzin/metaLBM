#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
namespace tt = boost::test_tools;

#include <iostream>

#include "metaLBM/DynamicArray.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestDynamicArray)

BOOST_AUTO_TEST_CASE(TestArchitectureCPU) {
  constexpr Architecture architecture = Architecture::CPU;
  typedef double dataT;
  DynamicArray<dataT, architecture> dArray1(5, (dataT) 1);
  DynamicArray<dataT, architecture> dArray2(5, (dataT) 2);

  BOOST_CHECK(dArray1[0] == (dataT) 1);

  BOOST_CHECK(dArray2[0] == (dataT) 2);

  dArray1.swap(dArray2);

  BOOST_CHECK(dArray1[0] == (dataT) 2);

  BOOST_CHECK(dArray2[0] == (dataT) 1);

}

BOOST_AUTO_TEST_SUITE_END()
