#include <boost/test/unit_test.hpp>
#define BOOST_MATH_INSTRUMENT

#include "solver.h"
#include "input.h"
#include "commons.h"
#include "force.h"

#include <array>
#include <vector>
#include <memory>
#include <math.h>


using namespace lbm;

BOOST_AUTO_TEST_SUITE(Equilibrium)

BOOST_AUTO_TEST_CASE(Equilibrium_ConservationMacro) {
  ELBMSolver solver;

  double density = 0.99;
  double velocityX = 0.0005;
  double velocityY = -0.000001;
  double velocity2 = velocityX * velocityX + velocityY * velocityY;

  std::array<double, 9> fEq;
  double densityR = 0.0, velocityXR = 0.0, velocityYR = 0.0;
  for(int i = 0; i < 9; ++i) {
    fEq[i] =  computeEquilibrium(i, density, velocityX, velocityY, velocity2);
    densityR += fEq[i];
    velocityXR += fEq[i] * celerityX[i];
    velocityYR += fEq[i] * celerityY[i];
  }

  BOOST_CHECK_SMALL(fabs(densityR - density)/density, 1e-10);
  BOOST_CHECK_SMALL(fabs(velocityXR - velocityX)/velocityX, 1e-10);
  BOOST_CHECK_SMALL(fabs(velocityYR - velocityY)/velocityY, 1e-10);
}


BOOST_AUTO_TEST_CASE(Equilibrium_Hfunction) {
  ELBMSolver solver;

  std::array<double, 9> fEq;

  double hFunction_fEq = 0.0, hFunctionDerivative_fEq = 0.0;
  for(int i = 0; i < 9; ++i) {
    fEq[i] =  computeEquilibrium(i, 1., 0., 0., 0.);
    hFunction_fEq += fEq[i] * log(fEq[i]/weight[i]);
    hFunctionDerivative_fEq += log(fEq[i]/weight[i]) + 1;
  }

  BOOST_CHECK_SMALL(hFunctionDerivative_fEq, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(ELBM)

BOOST_AUTO_TEST_CASE(ELBM_HStepFunctor_Case0) {
  ELBMSolver solver;

  std::array<double, 9> f;
  std::array<double, 9> fNeq;

  double deviation = 0.0;
  for(int i = 0; i < 9; ++i) {
    fNeq[i] =  deviation;
    f[i] =  deviation + computeEquilibrium(i, 1., 0., 0., 0.);
  }

  HStepFunctor hStepFunctor(f, fNeq);
  double alpha = 0.0;
  std::pair<double, double> hStepFunctorPair = hStepFunctor(alpha);

  BOOST_CHECK_SMALL(std::get<0>(hStepFunctorPair), 1e-8);
}


BOOST_AUTO_TEST_CASE(ELBM_SmallDeviation) {
  ELBMSolver solver;

  std::array<double, 9> f;
  std::array<double, 9> fNeq;

  double deviation = 0.00000001;
  for(int i = 0; i <9; ++i) {
    fNeq[i] =  deviation;
    f[i] =  deviation + computeEquilibrium(i, 1., 0., 0., 0.);
  }

  double alpha = solver.computeAlpha(f, fNeq, 2.0);
  BOOST_CHECK_EQUAL(alpha, 2.);
}

BOOST_AUTO_TEST_CASE(ForceNR_dELBM_SmallDeviation) {
  ForcedNR_ELBMSolver solver;

  std::array<double, 9> f;
  std::array<double, 9> fNeq;

  double deviation = 0.00000001;
  for(int i = 0; i <9; ++i) {
    fNeq[i] =  deviation;
    f[i] =  deviation + computeEquilibrium(i, 1., 0., 0., 0.);
  }

  double alpha = solver.computeAlpha(f, fNeq, 2.0);
  BOOST_CHECK_EQUAL(alpha, 2.);
}

BOOST_AUTO_TEST_CASE(ForceBNR_dELBM_SmallDeviation) {
  ForcedBNR_ELBMSolver solver;

  std::array<double, 9> f;
  std::array<double, 9> fNeq;

  double deviation = 0.00000001;
  for(int i = 0; i <9; ++i) {
    fNeq[i] =  deviation;
    f[i] =  deviation + computeEquilibrium(i, 1., 0., 0., 0.);
  }

  double alpha = solver.computeAlpha(f, fNeq, 2.0);
  BOOST_CHECK_EQUAL(alpha, 2.);
}


BOOST_AUTO_TEST_CASE(ELBM_MaxDeviation) {
  ELBMSolver solver;

  std::array<double, 9> f;
  std::array<double, 9> fNeq;

  double deviation = 0.001;
  for(int i = 0; i <9; ++i) {
    fNeq[i] =  deviation;
    f[i] =  fNeq[i] + computeEquilibrium(i, 1., 0., 0., 0.);
  }

  double alpha = solver.computeAlpha(f, fNeq, 2.0);
  BOOST_CHECK_EQUAL(alpha, 2.);
}

BOOST_AUTO_TEST_CASE(ForcedNR_ELBM_MaxDeviation) {
  ForcedNR_ELBMSolver solver;

  std::array<double, 9> f;
  std::array<double, 9> fNeq;

  double deviation = 0.001;
  for(int i = 0; i <9; ++i) {
    fNeq[i] =  deviation;
    f[i] =  fNeq[i] + computeEquilibrium(i, 1., 0., 0., 0.);
  }

  double alpha = solver.computeAlpha(f, fNeq, 2.0);
  BOOST_CHECK_EQUAL(alpha, 2.);
}

BOOST_AUTO_TEST_CASE(ForcedBNR_ELBM_MaxDeviation) {
  ForcedBNR_ELBMSolver solver;

  std::array<double, 9> f;
  std::array<double, 9> fNeq;

  double deviation = 0.001;
  for(int i = 0; i <9; ++i) {
    fNeq[i] =  deviation;
    f[i] =  fNeq[i] + computeEquilibrium(i, 1., 0., 0., 0.);
  }

  double alpha = solver.computeAlpha(f, fNeq, 2.0);
  BOOST_CHECK_EQUAL(alpha, 2.);
}


BOOST_AUTO_TEST_CASE(ELBM_HighDeviation) {
  ELBMSolver solver;

  std::array<double, 9> f;
  std::array<double, 9> fNeq;

  double deviation = 0.01;
  for(int i = 0; i <9; ++i) {
    fNeq[i] =  deviation;
    f[i] =  fNeq[i] + computeEquilibrium(i, 1., 0., 0., 0.);
  }

  double alpha = solver.computeAlpha(f, fNeq, 2.0);
  BOOST_CHECK_PREDICATE( std::not_equal_to<double>(), (alpha)(2.0) );
}

BOOST_AUTO_TEST_CASE(ForcedNR_ELBM_HighDeviation) {
  ForcedNR_ELBMSolver solver;

  std::array<double, 9> f;
  std::array<double, 9> fNeq;

  double deviation = 0.01;
  for(int i = 0; i <9; ++i) {
    fNeq[i] =  deviation;
    f[i] =  fNeq[i] + computeEquilibrium(i, 1., 0., 0., 0.);
  }

  double alpha = solver.computeAlpha(f, fNeq, 2.0);
  BOOST_CHECK_PREDICATE( std::not_equal_to<double>(), (alpha)(2.0) );
}

BOOST_AUTO_TEST_CASE(ForcedBNR_ELBM_HighDeviation) {
  ForcedBNR_ELBMSolver solver;

  std::array<double, 9> f;
  std::array<double, 9> fNeq;

  double deviation = 0.01;
  for(int i = 0; i <9; ++i) {
    fNeq[i] =  deviation;
    f[i] =  fNeq[i] + computeEquilibrium(i, 1., 0., 0., 0.);
  }

  double alpha = solver.computeAlpha(f, fNeq, 2.0);
  BOOST_CHECK_PREDICATE( std::not_equal_to<double>(), (alpha)(2.0) );
}


BOOST_AUTO_TEST_SUITE_END()
