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
#include "field_tmp.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestField)

BOOST_AUTO_TEST_CASE(TestFields) {
  vector<MathVector<valueType, L::dimQ>> fPop_field((L::lX_l+2*L::hX)
                                                    *(L::lY_l+2*L::hY)
                                                    *(L::lZ_l+2*L::hZ));

  for(int idx = 0; idx < L::lX_g *L::lY_g*L::lZ_g; ++idx) {
    for(int iQ = 0; iQ < L::dimQ; ++iQ) {
      fPop_field[idx][iQ] = (valueType) (idx + iQ);
    }
  }

  std::string fPop_name = "fPop";
  auto fPop_ptr = std::shared_ptr<Field<valueType>>(new DistributionField<valueType>(fPop_name,
                                                                                     fPop_field));
  constexpr int numberFields = 1;
  auto fieldsName = std::array<std::string, numberFields>{fPop_name};
  auto fieldsArray = std::array<std::shared_ptr<Field<valueType>>, numberFields>{fPop_ptr};

  Fields<valueType, numberFields> F(fieldsName, fieldsArray);

  BOOST_TEST(F(fPop_name)->fieldName == fPop_name);

  for(int idx = 0; idx < L::lX_g *L::lY_g*L::lZ_g; ++idx) {
    for(int iQ = 0; iQ < L::dimQ; ++iQ) {
      BOOST_TEST(F(fPop_name)->field[idx][iQ] == (valueType)(L::dimQ*idx + (L::dimQ-1)*L::dimQ/2),
                 tt::tolerance(1e-15));
    }
  }
}


BOOST_AUTO_TEST_SUITE_END()
