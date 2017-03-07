#ifndef INPUT_H
#define INPUT_H

#include <array>
#include <vector>
#include <string>

#undef ENABLE_LOG

namespace lbm {

  enum class SolverMethod {BGK, ELBM, error};

  enum class InitDensityType{homogeneous, peak, error};
  enum class InitVelocityType{homogeneous, perturbated, decay, error};

  enum class ForcingMethod {guo, EDM, error};
  enum class ForceType {constant, sinusoidal, turbulent, turbulent_time, error};

  enum class BoundaryType {bounceBack_halfWay, pressure_ZouHe, velocity_ZouHe, entropic, corner, error};
  enum class BoundaryPosition {top, bottom, left, right, top_left, top_right, bottom_right, bottom_left, error};
  enum class CornerPosition {top_left, top_right, bottom_right, bottom_left, error};

  enum class OutputType {vtr, backup, profile_density, profile_velocity_x, profile_velocity_y, profile_alpha, error};
  enum class ProfileConstant {x, y, error};

  constexpr int lengthX_g = 32;
  constexpr int lengthY_g = 32;
  constexpr int size_g = lengthX_g * lengthY_g;

  constexpr int lengthX_l = lengthX_g/NPROCS;
  constexpr int lengthY_l = lengthY_g;
  constexpr int size_l = lengthX_l * lengthY_l;
  constexpr int size_buf = lengthY_l * 3;
  constexpr int haloX = 1;
  constexpr int haloY = 1;

  constexpr bool restart = false;
  constexpr int iterationMax = 100;
  constexpr int writeStep = 10;
  constexpr int backupStep = 50;
  constexpr auto prefix = "test-32";

  constexpr double tau = 0.50000001;
  constexpr double beta = 1.0/(2.0*tau);
  constexpr SolverMethod solverMethod = SolverMethod::ELBM;
  constexpr ForcingMethod forcingMethod = ForcingMethod::guo;

  constexpr InitDensityType initDensityType = InitDensityType::homogeneous;
  constexpr double initDensityValue = 1.0;
  constexpr InitVelocityType initVelocityType = InitVelocityType::homogeneous;
  constexpr double initVelocityXValue = 0.0;
  constexpr double initVelocityYValue = 0.0;

  constexpr int numberForces = 1;
  constexpr std::array<ForceType, numberForces> forceTypeArray = { ForceType::turbulent_time };
  constexpr std::array<double, numberForces> forceAmplitudeXArray = { 2e-05 };
  constexpr std::array<double, numberForces> forceAmplitudeYArray = { 2e-05 };
  constexpr std::array<double, numberForces> forceWaveLengthXArray = { 0.0 };
  constexpr std::array<double, numberForces> forceWaveLengthYArray = { 0.0 };


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
  constexpr std::array<ProfileConstant, numberOutputs> outputProfileConstantArray = { ProfileConstant::x, ProfileConstant::x };
  constexpr std::array<double, numberOutputs> outputProfileConstantValueArray = { 0.5, 0.5 };

}

#endif // INPUT_H
