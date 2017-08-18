#ifndef INPUT_H
#define INPUT_H

#include <string>

#include "Options.h"
#include "MathVector.h"

namespace lbm {
  typedef DATA_TYPE dataT;

  constexpr LatticeType latticeT = LatticeType::D2Q9;
  constexpr int lengthX_g = 32;
  constexpr int lengthY_g = 32;
  constexpr int lengthZ_g = 8;

  constexpr MathVector<unsigned int, 3> length_g = {lengthX_g, lengthY_g, lengthZ_g};
  constexpr MathVector<unsigned int, 3> length_l = {lengthX_g/NPROCS, lengthY_g, lengthZ_g};

  constexpr MathVector<unsigned int, 3> process = {NPROCS, 1, 1};

  constexpr AlgorithmType algorithmT = AlgorithmType::Pull;
  constexpr PartitionningType partitionningT = PartitionningType::OneD;
  constexpr MemoryLayout memoryL = MemoryLayout::SoA;

  constexpr int startIteration = 0;
  constexpr int endIteration = 50;
  constexpr int writeStep = 50;
  constexpr int backupStep = 50;

  constexpr dataT relaxationTime = 0.65;
  constexpr CollisionType collisionT = CollisionType::BGK;
  constexpr EquilibriumType equilibriumT = EquilibriumType::Incompressible;

  constexpr InitDensityType initDensityT = InitDensityType::Homogeneous;
  constexpr dataT initDensityValue = 1.0;
  constexpr InitVelocityType initVelocityT = InitVelocityType::Homogeneous;
  constexpr MathVector<dataT, 3> initVelocityVector = { {0.0, 0.0, 0.0} };

  constexpr ForcingSchemeType forcingSchemeT = ForcingSchemeType::ExactDifferenceMethod;
  constexpr ForceType forceT = ForceType::Kolmogorov;

  constexpr MathVector<dataT, 3> forceAmplitude = { {1e-5, 0.0, 0.0} };
  constexpr MathVector<dataT, 3> forceWaveLength = { {8.0, 0.0, 0.0} };
  constexpr int minWavenumber = 0;
  constexpr int maxWavenumber = 0;

  constexpr BoundaryType boundaryT = BoundaryType::Generic;

  constexpr WriterType writerT= WriterType::VTR;
  constexpr auto prefix = "test-zxy";

  constexpr bool writeDensity = 1;
  constexpr bool writeVelocity = 1;
  constexpr bool writeForce = 1;
  constexpr bool writeEntropy = 1;
  constexpr bool writeAlpha = 1;
}

#endif // INPUT_H
