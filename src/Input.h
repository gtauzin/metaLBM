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
  constexpr int lengthZ_g = 1;

  constexpr unsigned int startIteration = 0;
  constexpr unsigned int endIteration = 10;
  constexpr unsigned int writeStep = 1;
  constexpr unsigned int successiveWriteStep = 2;
  constexpr unsigned int backUpStep = 100;

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
  constexpr Implementation implementationT = Implementation::MPI;
  constexpr MemoryLayout memoryL = MemoryLayout::SoA;

  constexpr dataT relaxationTime = 0.65;
  constexpr CollisionType collisionT = CollisionType::BGK;
  constexpr EquilibriumType equilibriumT = EquilibriumType::Incompressible;

  constexpr InitDensityType initDensityT = InitDensityType::Homogeneous;
  constexpr dataT initDensityValue = 1.0;
  constexpr InitVelocityType initVelocityT = InitVelocityType::Homogeneous;
  constexpr MathVector<dataT, 3> initVelocityVector = { {3.0, 4.0, 5.0} };

  constexpr ForcingSchemeType forcingSchemeT = ForcingSchemeType::ExactDifferenceMethod;
  constexpr ForceType forceT = ForceType::Kolmogorov;

  constexpr MathVector<dataT, 3> forceAmplitude = { {0.000001, 0.000001, 0.0} };
  constexpr MathVector<dataT, 3> forceWaveLength = { {32.0, 32.0, 0.0} };
  constexpr int forcekMin = 1;
  constexpr int forcekMax = 2;

  constexpr BoundaryType boundaryT = BoundaryType::Generic;

  constexpr InputOutput inputOutput = InputOutput::VTR;
  constexpr InputOutputType inputOutputType = InputOutputType::Serial;
  constexpr InputOutputDataFormat inputOutputDataFormat = InputOutputDataFormat::ascii;
  constexpr auto prefix = "test_mar_1";


  constexpr bool writeDensity = 1;
  constexpr bool writeVelocity = 1;
  constexpr bool writeForce = 1;
  constexpr bool writeEntropy = 1;
  constexpr bool writeAlpha = 1;
}

#endif // INPUT_H
