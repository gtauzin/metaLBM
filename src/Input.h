#ifndef INPUT_H
#define INPUT_H

#include <string>

#include "metaLBM/Commons.h"
#include "metaLBM/Options.h"
#include "metaLBM/MathVector.h"

namespace lbm {
  typedef double dataT;
  typedef MathVector<dataT, 3> Vector;

  constexpr LatticeType latticeT = LatticeType::D3Q19;
  constexpr int globalLengthX = 32;
  constexpr int globalLengthY = 32;
  constexpr int globalLengthZ = 32;

  constexpr unsigned int startIteration = 0;
  constexpr unsigned int endIteration = 500;
  constexpr unsigned int writeStep = 1001;
  constexpr unsigned int backUpStep = 1001;

  constexpr unsigned int scalarAnalysisStep = 1001;
  constexpr unsigned int spectralAnalysisStep = 1001;

  constexpr unsigned int successiveWriteStep = 2;

  constexpr AlgorithmType algorithmT = AlgorithmType::Pull;
  constexpr PartitionningType partitionningT = PartitionningType::OneD;
  constexpr Implementation implementationT = Implementation::MPI;
  constexpr MemoryLayout memoryL = MemoryLayout::SoA;

  constexpr dataT relaxationTime = 0.55;
  constexpr CollisionType collisionT = CollisionType::BGK;
  constexpr EquilibriumType equilibriumT = EquilibriumType::Incompressible;

  constexpr InitDensityType initDensityT = InitDensityType::Homogeneous;
  constexpr dataT initDensityValue = 1.0;
  constexpr InitVelocityType initVelocityT = InitVelocityType::Homogeneous;
  constexpr Vector initVelocityVector = { {0.0, 0.0, 0.0} };

  constexpr ForcingSchemeType forcingSchemeT = ForcingSchemeType::ShanChen;
  constexpr ForceType forceT = ForceType::Kolmogorov;

  constexpr Vector forceAmplitude = { {100, 100, 100} };
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
