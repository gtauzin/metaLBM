#include <boost/test/unit_test.hpp>
#define BOOST_MATH_INSTRUMENT

#include "force.h"
#include "input.h"
#include "commons.h"

#include <array>
#include <vector>
#include <memory>
#include <math.h>


using namespace lbm;

BOOST_AUTO_TEST_SUITE(Guo)

BOOST_AUTO_TEST_CASE(VelocityForcings) {
  ForcingGuo forcing;
  double density = 1.3;

  ConstantForce force(0.001, 0.002);

  forcing.forceX = force.forceX(0, 0, 0);
  forcing.forceY = force.forceY(0, 0, 0);

  double velocityXForcing = forcing.getVelocityXForcing(density);
  double velocityYForcing = forcing.getVelocityYForcing(density);

  BOOST_CHECK_EQUAL(velocityXForcing, 0.001/(2.0*1.3));
  BOOST_CHECK_EQUAL(velocityYForcing, 0.002/(2.0*1.3));
}


BOOST_AUTO_TEST_CASE(Collision_Forcing) {
  ForcingGuo forcing;
  ConstantForce force(0.001, 0.002);

  forcing.forceX = force.forceX(0, 0, 0);
  forcing.forceY = force.forceY(0, 0, 0);

  double density = 0.99;
  double velocityX = 0.0005;
  double velocityY = -0.000001;
  double velocity2 = velocityX * velocityX * velocityY * velocityY;

  std::array<double, 9> collisionForcing;
  double sum;
  for(int i = 0; i < 9; ++i) {
    collisionForcing[i] = forcing.getCollisionForcing(i, density,
                                                             velocityX, velocityY,
                                                             velocity2);
    sum += collisionForcing[i];
  }

  BOOST_CHECK_SMALL(sum, 1e-7);
}


BOOST_AUTO_TEST_SUITE_END()
