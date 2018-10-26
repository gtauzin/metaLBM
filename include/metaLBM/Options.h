#pragma once

#include "Commons.h"

namespace lbm {

  enum d { X, Y, Z };
  std::string dName = "XYZ";

  enum p { Re, Im };

  enum class LatticeType { Generic, D1Q3, D2Q5, D2Q9, D2Q13, D2Q17, D2Q21, D2Q37,
                           D3Q15, D3Q19, D3Q27, D3Q33 };

  enum class MemoryLayout { Generic, Default, SoA, AoS };
  enum class CommunicationType { Generic, MPI, NVSHMEM_OUT, NVSHMEM_IN };
  enum class PartitionningType { Generic, OneD, TwoD, ThreeD };
  enum class Architecture { Generic, CPU, GPU, CPUPinned };
  enum class Overlapping { Off, On };

 enum class DomainType { Generic, GlobalSpace, LocalSpace, HaloSpace, BufferXSpace,
                         GlobalFourier,  LocalFourier };

 enum class AlgorithmType { Generic, Pull, Push };

 enum class InitDensityType { Homogeneous, Peak };
 enum class InitVelocityType { Homogeneous, Perturbated, Wave, Decay };

 enum class FieldType { Generic, Density, Velocity, Force, Alpha, Entropy };

 enum class EquilibriumType { Generic, Exact, TruncationMa2, TruncationMa3 };

 enum class CollisionType { GenericSRT, BGK, ELBM, Approached_ELBM, Malaspinas_ELBM,
                            Essentially1_ELBM, Essentially2_ELBM,
                            ForcedNR_ELBM, ForcedBNR_ELBM, ForcedNR_ELBM_Forcing,
                            Regularized_BGK, Regularized_ELBM, Regularized_Approached_ELBM,
                            Regularized_Malaspinas_ELBM, Regularized_Essentially1_ELBM,
                            Regularized_Essentially2_ELBM, Regularized_ForcedNR_ELBM,
                            Regularized_ForcedBNR_ELBM, Regularized_ForcedNR_ELBM_Forcing };

 enum class ForcingSchemeType { Generic, None, Guo, ShanChen, ExactDifferenceMethod };
 enum class ForceType { None, Generic, GenericTimeIndependent, GenericTimeDependent,
                        Constant, Sinusoidal, Kolmogorov, ConstantShell, EnergyRemoval,
                        Turbulent2D };

 enum class BoundaryType { Generic, None, Periodic, BounceBack_Halfway, Entropic };

 enum class InputOutput { Generic, None, DAT, HDF5, XDMF };
 enum class InputOutputFormat { Generic, ascii, binary };

}  // namespace lbm
