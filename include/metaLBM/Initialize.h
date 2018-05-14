#pragma once

#include <iostream>
#include <string>

#include "Commons.h"
#include "Domain.h"
#include "DynamicArray.cuh"
#include "Equilibrium.h"
#include "Field.h"
#include "Force.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Options.h"
#include "Reader.h"

namespace lbm {

template <class T, Architecture architecture>
Field<T, 1, architecture, true> initLocalDensity(
    const unsigned int numberElements,
    const MathVector<int, 3>& rankMPI) {
  LBM_INSTRUMENT_ON("initLocalDensity<T>", 2)

  Field<T, 1, architecture, true>
      densityFieldR("density", numberElements, initDensityValue);

  switch (initDensityT) {
    case InitDensityType::Homogeneous: {
      break;
    }
    case InitDensityType::Peak: {
      if (rankMPI[d::X] == 0) {
        const T densityPeakValue = 3.0 * initDensityValue;
        Position center;

        center[d::X] =
            static_cast<unsigned int>((lSD::sLength()[d::X] - 1) * (T)0.4);
        center[d::Y] =
            static_cast<unsigned int>((lSD::sLength()[d::Y] - 1) * (T)0.3);
        center[d::Z] =
            static_cast<unsigned int>((lSD::sLength()[d::Z] - 1) * (T)0.2);

        densityFieldR.setLocalValue(center, densityPeakValue);
      }
      break;
    }
    default: { std::cout << "Wrong type of density initialization."; }
  }
  return densityFieldR;
}

template <class T, Architecture architecture>
Field<T, L::dimD, architecture, true> initLocalVelocity(
    const unsigned int numberElements) {
  LBM_INSTRUMENT_ON("initLocalVelocity<T>", 2)

  MathVector<T, L::dimD> initVelocityVectorProjected{{(T)0}};
  initVelocityVectorProjected = Project<T, T, L::dimD>::Do(initVelocityVector);

  Field<T, L::dimD, architecture, true> velocityFieldR(
      "velocity", numberElements, initVelocityVectorProjected);

  switch (initVelocityT) {
    case InitVelocityType::Homogeneous: {
      break;
    }

    default: { std::cout << "Wrong type of velocity initialization."; }
  }
  return velocityFieldR;
}

template <class T, Architecture architecture>
Field<T, L::dimD, architecture, writeForce> initLocalForce(
    const unsigned int numberElements,
    const MathVector<int, 3>& rankMPI) {
  LBM_INSTRUMENT_ON("initLocalForce<T>", 2)

  Field<T, L::dimD, architecture, writeForce>
      forceFieldR("force", numberElements, 0);

  Force<T, forceT> force(forceAmplitude, forceWaveLength, forcekMin, forcekMax);

  if (writeForce) {
    force.setLocalForceArray(forceFieldR.getLocalData(), numberElements,
                             gFD::offset(rankMPI));
  }

  return forceFieldR;
}

template <class T, Architecture architecture>
Field<T, 1, architecture, writeAlpha> initLocalAlpha(
    const unsigned int numberElements) {
  LBM_INSTRUMENT_ON("initLocalAlpha<T>", 2)

  Field<T, 1, architecture, writeAlpha>
      alphaFieldR("alpha", numberElements, (T)2);
  return alphaFieldR;
}

template <class T, Architecture architecture>
Distribution<T, architecture> initLocalDistribution(
    const Field<T, 1, architecture, true>& densityField,
    const Field<T, L::dimD, architecture, true>& velocityField,
    const MathVector<int, 3>& rankMPI) {
  LBM_INSTRUMENT_ON("initLocalDistribution<T>", 2)

  Distribution<T, architecture> distributionR(densityField.numberElements);

  if (startIteration == 0) {
    Position iP;
    for (auto iZ = lSD::sStart()[d::Z]; iZ < lSD::sEnd()[d::Z]; iZ++) {
      for (auto iY = lSD::sStart()[d::Y]; iY < lSD::sEnd()[d::Y]; iY++) {
        for (auto iX = lSD::sStart()[d::X]; iX < lSD::sEnd()[d::X]; iX++) {
          iP = Position({iX, iY, iZ});

          T density = densityField.getLocalValue(iP);
          MathVector<T, L::dimD> velocity = velocityField.getLocalVector(iP);
          T velocity2 = velocity.norm2();

          for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
            distributionR.setLocalValue(
                iP, Equilibrium_::calculate(density, velocity, velocity2, iQ),
                iQ);
          }
        }
      }
    }

  }

  else {
    DistributionReader_ distributionReader(prefix, rankMPI);
    distributionReader.openFile(startIteration);
    distributionReader.readDistribution(distributionR);
    distributionReader.closeFile();
  }

  return distributionR;
}

}  // namespace lbm
