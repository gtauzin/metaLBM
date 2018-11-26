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
  Field<T, 1, architecture, true> initDensity(const Stream<architecture>& stream) {
    LBM_INSTRUMENT_ON("initDensity<T>", 2)

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
          static_cast<unsigned int>((lSD::sLength()[d::X] - 1) * (T)0.2);
        center[d::Y] =
          static_cast<unsigned int>((lSD::sLength()[d::Y] - 1) * (T)0.5);
        center[d::Z] =
          static_cast<unsigned int>((lSD::sLength()[d::Z] - 1) * (T)0.5);

        densityFieldR.setValue(center, densityPeakValue, FFTWInit::numberElements);
      }
      break;
    }
    default: { std::cout << "Wrong type of density initialization."; }
    }
    return densityFieldR;
  }

  template <class T, Architecture architecture>
  Field<T, L::dimD, architecture, true> initVelocity(
    const Stream<architecture>& stream) {
    LBM_INSTRUMENT_ON("initVelocity<T>", 2)

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
  Field<T, L::dimD, architecture, writeForce> initForce(const Stream<architecture>& stream) {
    LBM_INSTRUMENT_ON("initForce<T>", 2)

    Field<T, L::dimD, architecture, writeForce>
      forceFieldR("force", 0, stream);

    return forceFieldR;
  }

  template <class T, Architecture architecture>
  Field<T, 1, architecture, writeAlpha> initAlpha(const Stream<architecture>& stream) {
    LBM_INSTRUMENT_ON("initAlpha<T>", 2)

    Field<T, 1, architecture, writeAlpha>
      alphaFieldR("alpha", (T)2, stream);
    return alphaFieldR;
  }

  template <class T, Architecture architecture>
  Distribution<T, architecture> initDistribution(Field<T, 1, architecture, true>& densityField,
                                                 Field<T, L::dimD, architecture, true>& velocityField,
                                                 Field<T, L::dimD, architecture, true>& strainDiagonal,
                                                 Field<T, 2 * L::dimD - 3, architecture, true>& strainSymmetric,
                                                 const Stream<architecture>& stream) {
    LBM_INSTRUMENT_ON("initDistribution<T>", 2)

    Distribution<T, architecture> distributionR;
    DistributionReader_ distributionReader(std::string(prefix), "distribution_restart");
    FieldReader_ fieldReader(std::string(prefix), "field_restart");

    Computation<architecture, L::dimD> computationLocal(lSD::sStart(), lSD::sEnd());
    unsigned int numberElements = FFTWInit::numberElements;
    T * distributionPtr = distributionR.getData(FFTWInit::numberElements);
    T * densityPtr = densityField.getData(FFTWInit::numberElements);
    T * velocityPtr = velocityField.getData(FFTWInit::numberElements);

    if (startIteration == 0) {
      computationLocal.Do(stream, [=] LBM_HOST LBM_DEVICE (const Position& iP) {
        T density = densityField.getValue(iP, numberElements);
        MathVector<T, L::dimD> velocity = velocityField.getVector(iP, numberElements);
        T velocity2 = velocity.norm2();

        for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
          distributionPtr[iQ * numberElements + lSD::getIndex(iP)] =
            Equilibrium_::calculate(density, velocity, velocity2, iQ);
        }
        });
    }

    else if(distributionReader.fileExists(startIteration)) {
      distributionReader.openFile(startIteration);
      distributionReader.readDistribution(distributionR);
      distributionReader.closeFile();

      computationLocal.Do(stream, [=] LBM_HOST LBM_DEVICE (const Position& iP) {
          unsigned int indexLocal = lSD::getIndex(iP);

          densityPtr[indexLocal] = 0;

          for(auto iD = 0; iD < L::dimD; ++iD) {
            (velocityPtr + iD * numberElements)[indexLocal] = 0;
          }

          for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
            densityPtr[indexLocal] +=
              distributionPtr[iQ * numberElements + indexLocal];

            for(auto iD = 0; iD < L::dimD; ++iD) {
              (velocityPtr + iD * numberElements)[indexLocal] += L::celerity()[iQ][iD]
                * distributionPtr[iQ * numberElements + indexLocal];
            }
          }

        });
    }

    else if(fieldReader.fileExists(startIteration)) {
      fieldReader.openFile(startIteration);
      fieldReader.readField(densityField);
      fieldReader.readField(velocityField);
      fieldReader.closeFile();

      //velocityPtr[lSD::getIndex({1,1,0})] = 10;

      const MathVector<MathVector<dataT, L::dimD>, L::dimQ> qDiagonal = calculateQDiagonal();
      const MathVector<MathVector<dataT, 2*L::dimD-3>, L::dimQ> qSymmetric = calculateQSymmetric();

      T * strainDiagonalPtr = strainDiagonal.getData(FFTWInit::numberElements);
      StrainDiagonalCalculator<double, Architecture::CPU, PartitionningType::OneD, L::dimD>
        strainDiagonalCalculator(velocityPtr, strainDiagonalPtr, globalLengthPtrdiff_t, gFD::offset(MPIInit::rank));
      strainDiagonalCalculator.executeSpace();

      T * strainSymmetricPtr = strainSymmetric.getData(FFTWInit::numberElements);
      StrainSymmetricCalculator<double, Architecture::CPU, PartitionningType::OneD, L::dimD>
        strainSymmetricCalculator(velocityPtr, strainSymmetricPtr, globalLengthPtrdiff_t, gFD::offset(MPIInit::rank));
      strainSymmetricCalculator.executeSpace();

      computationLocal.Do(stream, [=] LBM_HOST LBM_DEVICE (const Position& iP) {
          unsigned int indexLocal = lSD::getIndex(iP);

          T density = densityField.getValue(iP, numberElements);
          MathVector<T, L::dimD> velocity = velocityField.getVector(iP, numberElements);
          T velocity2 = velocity.norm2();

          for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
            T nonEquilibrium_iQ = 0;
            for (auto iD = 0; iD < L::dimD; ++iD) {
              nonEquilibrium_iQ += qDiagonal[iQ][iD] * (strainDiagonalPtr + iD * numberElements)[indexLocal];
            }

            for (auto iD = 0; iD < 2 * L::dimD - 3; ++iD) {
              nonEquilibrium_iQ += 2 * qSymmetric[iQ][iD] * (strainSymmetricPtr + iD * numberElements)[indexLocal];
            }

            nonEquilibrium_iQ *= - L::weight()[iQ] * relaxationTime * densityPtr[indexLocal] / L::cs2;

            distributionPtr[iQ * numberElements + lSD::getIndex(iP)]
              = Equilibrium_::calculate(density, velocity, velocity2, iQ) + nonEquilibrium_iQ ;
          }

        });
    }

    else {
      std::cout << "No file top restart from. Please check that filename you wish to restart from is : "
                << distributionReader.getFileName(startIteration) << " or "
                << fieldReader.getFileName(startIteration) << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return distributionR;
  }

}  // namespace lbm
