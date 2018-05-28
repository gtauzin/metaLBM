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
                                                   const Stream<architecture>& stream) {
    LBM_INSTRUMENT_ON("initLocalDensity<T>", 2)

    Field<T, 1, architecture, true>
    densityFieldR("density", initDensityValue, stream);

    switch (initDensityT) {
    case InitDensityType::Homogeneous: {
      break;
    }
    case InitDensityType::Peak: {
      if (MPIInit::rank[d::X] == 0) {
        const T densityPeakValue = 3.0 * initDensityValue;
        Position center;

        center[d::X] =
          static_cast<unsigned int>((lSD::sLength()[d::X] - 1) * (T)0.4);
        center[d::Y] =
          static_cast<unsigned int>((lSD::sLength()[d::Y] - 1) * (T)0.3);
        center[d::Z] =
          static_cast<unsigned int>((lSD::sLength()[d::Z] - 1) * (T)0.2);

        densityFieldR.setLocalValue(center, densityPeakValue, FFTWInit::numberElements);
      }
      break;
    }
    default: { std::cout << "Wrong type of density initialization."; }
    }
    return densityFieldR;
  }

  template <class T, Architecture architecture>
  Field<T, L::dimD, architecture, true> initLocalVelocity(
    const Stream<architecture>& stream) {
    LBM_INSTRUMENT_ON("initLocalVelocity<T>", 2)

    MathVector<T, L::dimD> initVelocityVectorProjected{{(T)0}};
    initVelocityVectorProjected = Project<T, T, L::dimD>::Do(initVelocityVector);

    Field<T, L::dimD, architecture, true> velocityFieldR(
      "velocity", initVelocityVectorProjected, stream);

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
    const Stream<architecture>& stream) {
    LBM_INSTRUMENT_ON("initLocalForce<T>", 2)

    Field<T, L::dimD, architecture, writeForce>
      forceFieldR("force", 0, stream);

    Force<T, forceT> force(forceAmplitude, forceWaveLength, forcekMin, forcekMax);

    if (writeForce) {
      force.setLocalForceArray(forceFieldR.getLocalData(FFTWInit::numberElements),
                               gFD::offset(MPIInit::rank));
    }

    return forceFieldR;
  }

  template <class T, Architecture architecture>
  Field<T, 1, architecture, writeAlpha> initLocalAlpha(
    const Stream<architecture>& stream) {
    LBM_INSTRUMENT_ON("initLocalAlpha<T>", 2)

    Field<T, 1, architecture, writeAlpha>
      alphaFieldR("alpha", (T)2, stream);
    return alphaFieldR;
  }

  template <class T, Architecture architecture>
  Distribution<T, architecture> initLocalDistribution(
    Field<T, 1, architecture, true>& densityField,
    Field<T, L::dimD, architecture, true>& velocityField,
    const Stream<architecture>& stream) {
    LBM_INSTRUMENT_ON("initLocalDistribution<T>", 2)

    Distribution<T, architecture> distributionR;

    Computation<architecture, L::dimD> computationLocal(lSD::sStart(), lSD::sEnd());
    unsigned int numberElements = FFTWInit::numberElements;
    T * localDistribution = distributionR.getLocalData(FFTWInit::numberElements);
    T * localDensity = densityField.getLocalData(FFTWInit::numberElements);
    T * localVelocity = velocityField.getLocalData(FFTWInit::numberElements);

    if (startIteration == 0) {
      computationLocal.Do(stream, [=] LBM_HOST LBM_DEVICE (const Position& iP) {
        T density = densityField.getLocalValue(iP, numberElements);
        MathVector<T, L::dimD> velocity = velocityField.getLocalVector(iP, numberElements);
        T velocity2 = velocity.norm2();

        for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
          localDistribution[iQ * numberElements + lSD::getIndex(iP)] =
            Equilibrium_::calculate(density, velocity, velocity2, iQ);
        }
        });
    }

    else {
      DistributionReader_ distributionReader(prefix);
      distributionReader.openFile(startIteration);
      distributionReader.readDistribution(distributionR);
      distributionReader.closeFile();

      computationLocal.Do(stream, [=] LBM_HOST LBM_DEVICE (const Position& iP) {
          unsigned int indexLocal = lSD::getIndex(iP);

          localDensity[indexLocal] = 0;

          for(auto iD = 0; iD < L::dimD; ++iD) {
            (localVelocity + iD * numberElements)[indexLocal] = 0;
          }

          for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
            localDensity[indexLocal] +=
              localDistribution[iQ * numberElements + indexLocal];

            for(auto iD = 0; iD < L::dimD; ++iD) {
              (localVelocity + iD * numberElements)[indexLocal] += L::celerity()[iQ][iD]
                * localDistribution[iQ * numberElements + indexLocal];
            }
          }

        });
    }

    return distributionR;
  }

}  // namespace lbm
