#ifndef INPUT_H
#define INPUT_H

#include <array>
#include <vector>
#include <string>

/* definition to expand macro then apply to pragma message */
#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#define VAR_NAME_VALUE(var) #var "="  VALUE(var)

/* Some example here */
//#pragma message(VAR_NAME_VALUE(NPROCS))



namespace lbm {
  enum class LatticeType {D1Q3, D2Q9, D3Q27};

  enum class SolverMethod {BGK, ELBM, Approached_ELBM, ForcedNR_ELBM, ForcedBNR_ELBM, error};

  enum class InitDensityType{homogeneous, peak, error};
  enum class InitVelocityType{homogeneous, perturbated, wave, decay, error};

  enum class ForcingMethod {guo, SC, EDM, error};
  enum class ForceType {constant, sinusoidal, turbulent, turbulent_time, wavenumbers_time, error};

  enum class BoundaryType {bounceBack_halfWay, pressure_ZouHe, velocity_ZouHe, entropic, corner, error};
  enum class BoundaryPosition {top, bottom, left, right, top_left, top_right, bottom_right, bottom_left, error};
  enum class CornerPosition {top_left, top_right, bottom_right, bottom_left, error};

  enum class OutputType {vtr, backup, error};

  typedef DATA_TYPE valueType;

  constexpr LatticeType latticeType = LatticeType::D2Q9;
  constexpr int lengthX_g = 16;
  constexpr int lengthY_g = 16;
  constexpr int lengthZ_g = 1;
  constexpr int length_g = lengthX_g;
  constexpr int size_g = lengthX_g * lengthY_g * lengthZ_g;

  constexpr int lengthX_l = lengthX_g/NPROCS;
  constexpr int lengthY_l = lengthY_g;
  constexpr int lengthZ_l = lengthZ_g;
  constexpr int size_l = lengthX_l * lengthY_l * lengthZ_l;

  constexpr int startIteration = 0;
  constexpr int iterationMax = 1000;
  constexpr int writeStep = 1;
  constexpr int backupStep = 1;
  constexpr auto prefix = "test";

  constexpr double tau = 0.65;
  constexpr double beta = 1.0/(2.0*tau);
  constexpr SolverMethod solverMethod = SolverMethod::ForcedNR_ELBM;

  constexpr InitDensityType initDensityType = InitDensityType::homogeneous;
  constexpr double initDensityValue = 1.0;
  constexpr InitVelocityType initVelocityType = InitVelocityType::homogeneous;
  constexpr double initVelocityXValue = 0.0;
  constexpr double initVelocityYValue = 0.0;
  constexpr double initVelocityZValue = 0.0;

  constexpr int numberForces = 1;
  constexpr ForcingMethod forcingMethod = ForcingMethod::EDM;
  constexpr std::array<ForceType, numberForces> forceTypeArray = { ForceType::constant };
  constexpr std::array<double, numberForces> forceAmplitudeXArray = { 1e-05 };
  constexpr std::array<double, numberForces> forceAmplitudeYArray = { 1e-05 };
  constexpr std::array<double, numberForces> forceWaveLengthXArray = { 64.0 };
  constexpr std::array<double, numberForces> forceWaveLengthYArray = { 64.0 };
  constexpr int numberWavenumberPairs = 1;
  constexpr std::array<std::array<int, numberWavenumberPairs>, numberForces> forcekXArray = { {0} };
  constexpr std::array<std::array<int, numberWavenumberPairs>, numberForces> forcekYArray = { {0} };

  constexpr std::array<double, numberForces> forcePeriodTimeArray = { 0.08 };

  constexpr int numberBCs = 0;
  constexpr std::array<BoundaryType, numberBCs> boundaryTypeArray = {  };
  constexpr std::array<BoundaryPosition, numberBCs> boundaryPositionArray = {  };
  constexpr std::array<int, numberBCs> boundaryStartArray = {  };
  constexpr std::array<int, numberBCs> boundaryEndArray = {  };
  constexpr std::array<double, numberBCs> boundaryPressureArray = {  };
  constexpr std::array<double, numberBCs> boundaryVelocityXArray = {  };
  constexpr std::array<double, numberBCs> boundaryVelocityYArray = {  };

  constexpr int numberOutputs = 2;
  constexpr std::array<OutputType, numberOutputs> outputTypeArray = { OutputType::backup, OutputType::vtr };

  constexpr bool writeNextDensity = 1;
  constexpr bool writePreviousDensity = 1;
  constexpr bool writeNextVelocity = 1;
  constexpr bool writePreviousVelocity = 1;
  constexpr bool writeNextVorticity = 1;
  constexpr bool writePreviousVorticity = 1;
  constexpr bool writeNextAlpha = 1;

}

#endif // INPUT_H
