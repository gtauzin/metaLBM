#ifndef INPUT_H
#define INPUT_H

#include <string>

#include "metaLBM/Commons.h"
#include "metaLBM/Options.h"
#include "metaLBM/MathVector.h"

namespace lbm {
  typedef double dataT;

  constexpr LatticeType latticeT = LatticeType::D2Q9;
  constexpr int lengthX_g = 2;
  constexpr int lengthY_g = 2;
  constexpr int lengthZ_g = 32;

  constexpr unsigned int startIteration = 0;
  constexpr unsigned int endIteration = 50;
  constexpr unsigned int writeStep = 1;

  constexpr unsigned int successiveWriteStep = 2;
  constexpr unsigned int backUpStep = 1;

  constexpr AlgorithmType algorithmT = AlgorithmType::Pull;
  constexpr PartitionningType partitionningT = PartitionningType::OneD;
  constexpr Implementation implementationT = Implementation::MPI;
  constexpr MemoryLayout memoryL = MemoryLayout::AoS;

  constexpr dataT relaxationTime = 0.85;
  constexpr CollisionType collisionT = CollisionType::BGK;
  constexpr EquilibriumType equilibriumT = EquilibriumType::Incompressible;

  constexpr InitDensityType initDensityT = InitDensityType::Homogeneous;
  constexpr dataT initDensityValue = 1.0;
  constexpr InitVelocityType initVelocityT = InitVelocityType::Homogeneous;
  constexpr MathVector<dataT, 3> initVelocityVector = { {0.0, 0.0, 0.0} };

  constexpr ForcingSchemeType forcingSchemeT = ForcingSchemeType::ExactDifferenceMethod;
  constexpr ForceType forceT = ForceType::Kolmogorov;

  constexpr MathVector<dataT, 3> forceAmplitude = { {0.0, 0.0, 0.0} };
  constexpr MathVector<dataT, 3> forceWaveLength = { {32.0, 32.0, 0.0} };
  constexpr int forcekMin = 1;
  constexpr int forcekMax = 2;

  constexpr BoundaryType boundaryT = BoundaryType::Generic;

  constexpr InputOutput inputOutput = InputOutput::VTR;
  constexpr InputOutputType inputOutputType = InputOutputType::Serial;
  constexpr InputOutputDataFormat inputOutputDataFormat = InputOutputDataFormat::ascii;
  constexpr auto prefix = "test_serial";

  constexpr bool writeDensity = 1;
  constexpr bool writeVelocity = 1;
  constexpr bool writeForce = 1;
  constexpr bool writeEntropy = 1;
  constexpr bool writeAlpha = 1;
}

#endif // INPUT_H
