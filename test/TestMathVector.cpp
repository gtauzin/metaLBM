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

#include "input.h"
#include "commons.h"
#include "structure.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestStructure)

BOOST_AUTO_TEST_CASE(TestUnrollFor) {
  constexpr int size = 100;
  vector<valueType> unrolledLoopVector_r(size);
  vector<valueType> forLoopVector_r(size);

  UnrolledFor<0, size>::Do([&] (int i) {
      unrolledLoopVector_r[i] = (valueType) i+1;
    });

  for(int i = 0; i < size; ++i){
    forLoopVector_r[i] = (valueType) i+1;
  }


  for(int i = 0; i < size; ++i){
    BOOST_CHECK_EQUAL(unrolledLoopVector_r[i], forLoopVector_r[i]);
  }

}

BOOST_AUTO_TEST_SUITE_END()
