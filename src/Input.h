#ifndef INPUT_H
#define INPUT_H

#include <string>

#include "metaLBM/Commons.h"
#include "metaLBM/Options.h"
#include "metaLBM/MathVector.h"

namespace lbm {
  typedef double dataT;
  typedef MathVector<dataT, 3> Vector;

  constexpr LatticeType latticeT = LatticeType::D3Q27;
  constexpr int globalLengthX = 16;
  constexpr int globalLengthY = 16;
  constexpr int globalLengthZ = 16;

  constexpr unsigned int startIteration = 0;
  constexpr unsigned int endIteration = 10;
  constexpr unsigned int writeStep = 1;
  constexpr unsigned int backUpStep = 10;

  constexpr unsigned int scalarAnalysisStep = 1;
  constexpr unsigned int spectralAnalysisStep = 2;

  constexpr unsigned int successiveWriteStep = 2;

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
  constexpr Vector initVelocityVector = { {0.0, 0.0, 0.0} };

  constexpr ForcingSchemeType forcingSchemeT = ForcingSchemeType::ExactDifferenceMethod;
  constexpr ForceType forceT = ForceType::ConstantShell;

  constexpr Vector forceAmplitude = { {0.00001, 0.00001, 0.00001} };
  constexpr Vector forceWaveLength = { {32.0, 32.0, 32.0} };
  constexpr int forcekMin = 1;
  constexpr int forcekMax = 2;

  constexpr BoundaryType boundaryT = BoundaryType::Generic;

  constexpr InputOutputFormat inputOutputFormatT = InputOutputFormat::ascii;
  constexpr auto prefix = "test_serial";

  constexpr bool writeForce = 1;
  constexpr bool writeEntropy = 1;
  constexpr bool writeAlpha = 0;
  constexpr bool writeVorticity = 1;

  constexpr bool analyzeTotalEnergy = 1;
  constexpr bool analyzeTotalEnstrophy = 1;

  constexpr bool analyzeEnergySpectra = 1;
  constexpr bool analyzeEnstrophySpectra = 1;
}

#endif // INPUT_H
