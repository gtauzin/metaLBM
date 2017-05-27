#ifndef INPUT_H
#define INPUT_H

#include <string>

#include "Options.h"
#include "StaticArray.h"

namespace lbm {
  typedef DATA_TYPE dataType;
  constexpr LatticeType latticeType = LatticeType::D1Q3;
  constexpr int lengthX_g = 16;
  constexpr int lengthY_g = 16;
  constexpr int lengthZ_g = 16;

  constexpr StaticArray<int, 3> length_g = {lengthX_g, lengthY_g, lengthZ_g};

  constexpr StaticArray<int, 3> length_g = {NPROCS, 1, 1};

  constexpr PartitionningType partitionningT = PartitionningType::OneD;

  constexpr int startIteration = 0;
  constexpr int endIteration = 10;
  constexpr int writeStep = 1;
  constexpr int backupStep = 1;
  constexpr auto prefix = "test";

  constexpr dataType tau = 0.65;
  constexpr dataType beta = 1.0/(2.0*tau);
  constexpr CollisionType collisionType = CollisionType::ForcedNR_ELBM;

  constexpr InitDensityType initDensityType = InitDensityType::Homogeneous;
  constexpr dataType initDensityValue = 1.0;
  constexpr InitVelocityType initVelocityType = InitVelocityType::Homogeneous;
  constexpr StaticArray<StaticArray<dataType, 3>, 1> initVelocityValue = { {0.0, 0.0, 0.0} };

  constexpr int numberForces = 1;
  constexpr ForcingSchemeType forcingSchemeType = ForcingSchemeType::ExactDifferenceMethod;
  constexpr StaticArray<ForceType, numberForces> forceTypeArray = { ForceType::Kolmogorov };
  constexpr StaticArray<StaticArray<dataType, 3>, numberForces> forceAmplitudeArray = { {1e-5, 0.0, 0.0}, {1e-5, 0.0, 0.0}, {1e-5, 0.0, 0.0} };

  constexpr StaticArray<StaticArray<dataType, 3>, numberForces> forceWaveLengthArray = { {4.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {4.0, 0.0, 0.0} };

  constexpr StaticArray<int, numberForces> minWavenumber = { 0 };
  constexpr StaticArray<int, numberForces> maxWavenumber = { 0 };


  constexpr int numberWriters = 2;
  constexpr StaticArray<WriterType, numberWriters> writerTypeArray = { WriterType::Backup, WriterType::VTR };

  constexpr bool writeDensity = 3;
  constexpr bool writeVelocity = 3;
  constexpr bool writeVorticity = 1;
  constexpr bool writeEntropy = 2;
  constexpr bool writeAlpha = 1;

}

#endif // INPUT_H
