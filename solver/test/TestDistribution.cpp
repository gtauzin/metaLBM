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
#include "distribution.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestDistribution)


BOOST_AUTO_TEST_SUITE_END()
