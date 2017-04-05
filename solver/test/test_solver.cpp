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
#include "solver.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestSolver)

BOOST_AUTO_TEST_CASE(TestBGK) {
  auto solver = std::shared_ptr<Solver<valueType>>(new BGK<valueType>());

  constexpr int numberTest = 1000;
  vector<MathVector<valueType, L::dimQ>> f(numberTest,
                                           MathVector<valueType, L::dimQ>{0});
  vector<MathVector<valueType, L::dimQ>> fNeq(numberTest,
                                              MathVector<valueType, L::dimQ>{0});
  vector<valueType> alpha_r(numberTest);

  for(int i = 0; i < numberTest; ++i) {
    f[i] = MathVector<valueType, L::dimQ>{(valueType) i/ (valueType) numberTest};
    fNeq[i] = MathVector<valueType, L::dimQ>{1 - 2.0* (valueType) i / (valueType) numberTest};
    valueType alphaGuess = (valueType) i;
    alpha_r[i] = solver->computeAlpha(f[i], fNeq[i], alphaGuess);
  }

    for(int i = 0; i < numberTest; ++i) {
      BOOST_TEST(alpha_r[i] == (valueType) 2,
             tt::tolerance(1e-15));

    }
}

BOOST_AUTO_TEST_SUITE_END()
