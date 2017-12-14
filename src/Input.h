#ifndef INPUT_H
#define INPUT_H

#include <string>

#include "metaLBM/Commons.h"
#include "metaLBM/Options.h"
#include "metaLBM/MathVector.h"

namespace lbm {
  typedef DATA_TYPE dataT;

  constexpr LatticeType latticeT = LatticeType::D2Q9;
  constexpr int lengthX_g = 16;
  constexpr int lengthY_g = 16;
  constexpr int lengthZ_g = 16;

  constexpr unsigned int startIteration = 0;
  constexpr unsigned int endIteration = 10;
  constexpr unsigned int writeStep = 10;//endIteration+1;
  constexpr unsigned int successiveWriteStep = 2;
  constexpr unsigned int backupStep = endIteration+1;

  DEVICE HOST
  constexpr MathVector<unsigned int, 3> length_g() {
    return MathVector<unsigned int, 3>({lengthX_g, lengthY_g, lengthZ_g});
  }

  DEVICE HOST
    constexpr MathVector<unsigned int, 3> length_l() {
    return MathVector<unsigned int, 3>({lengthX_g/NPROCS, lengthY_g, lengthZ_g});
  }

  constexpr AlgorithmType algorithmT = AlgorithmType::Pull;
  constexpr PartitionningType partitionningT = PartitionningType::OneD;
  constexpr MemoryLayout memoryL = MemoryLayout::AoS;

  constexpr dataT relaxationTime = 0.65;
  constexpr CollisionType collisionT = CollisionType::BGK;
  constexpr EquilibriumType equilibriumT = EquilibriumType::Incompressible;

  constexpr InitDensityType initDensityT = InitDensityType::Peak;
  constexpr dataT initDensityValue = 1.0;
  constexpr InitVelocityType initVelocityT = InitVelocityType::Homogeneous;
  constexpr MathVector<dataT, 3> initVelocityVector = { {0.005, 0.002, 0.0} };

  constexpr ForcingSchemeType forcingSchemeT = ForcingSchemeType::ExactDifferenceMethod;
  constexpr ForceType forceT = ForceType::Kolmogorov;

  constexpr MathVector<dataT, 3> forceAmplitude = { {0.01, 0.01, 0.0} };
  constexpr MathVector<dataT, 3> forceWaveLength = { {8.0, 16.0, 0.0} };
  constexpr int minWavenumber = 0;
  constexpr int maxWavenumber = 0;

  constexpr BoundaryType boundaryT = BoundaryType::Generic;

  constexpr WriterType writerT = WriterType::VTR;
  constexpr auto prefix = "test_fields";
  constexpr WriterFileFormat writerFileFormatT = WriterFileFormat::ascii;

  constexpr bool writeDensity = 1;
  constexpr bool writeVelocity = 1;
  constexpr bool writeForce = 1;
  constexpr bool writeEntropy = 1;
  constexpr bool writeAlpha = 1;
}

#endif // INPUT_H
