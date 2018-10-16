#define BOOST_TEST_MODULE "C++ Unit Tests for metaLBM"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#define NPROCS 1
#include "Lattice.h"
typedef double dataType;
typedef lbm::Lattice<dataType, lbm::LatticeType::D1Q3> L;
#include "Domain.h"
#include "Field.h"
#include "Initialize.h"

using namespace lbm;

BOOST_AUTO_TEST_SUITE(TestInitializeT)

BOOST_AUTO_TEST_CASE(initGlobalAlphaT) {
  ::lbm::Field<dataType, 1, true> alphaField("alpha");

  auto alphaGlobalField = ::lbm::initGlobalAlpha<dataType>();

  alphaField.setGlobalField(alphaGlobalField);


  for(unsigned int i = 0; i < gD::volume(); ++i) {
  const auto alpha_ = alphaField.getGlobalField(i, 0);
    BOOST_CHECK_EQUAL(alpha_, (dataType) 2);
  }

  }

BOOST_AUTO_TEST_SUITE_END()
