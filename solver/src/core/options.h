#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <array>
#include <vector>
#include <string>


namespace lbm {
  enum class LatticeType {D1Q3, D2Q5, D2Q9, D3Q15, D3Q19, D3Q27};

  enum class SolverMethod {BGK, ELBM, Approached_ELBM, ForcedNR_ELBM, ForcedBNR_ELBM, error};

  enum class InitDensityType{homogeneous, peak, error};
  enum class InitVelocityType{homogeneous, perturbated, wave, decay, error};

  enum class ForcingSchemeMethod {guo, SC, EDM, error};
  enum class ForceType {constant, sinusoidal, kolmogorov, turbulent, turbulent_time, wavenumbers_time, error};

  enum class BoundaryType {bounceBack_halfWay, pressure_ZouHe, velocity_ZouHe, entropic, corner, error};
  enum class BoundaryPosition {top, bottom, left, right, top_left, top_right, bottom_right, bottom_left, error};
  enum class CornerPosition {top_left, top_right, bottom_right, bottom_left, error};

  enum class OutputType {vtr, backup, error};
}

#endif // PARAMETERS_H
