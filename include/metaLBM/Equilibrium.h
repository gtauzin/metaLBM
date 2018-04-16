#pragma once

#include "Commons.h"
#include "Helpers.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Options.h"

namespace lbm {

template <class T, LatticeType latticeType, EquilibriumType equilibriumType>
class Equilibrium {};

template <class T, LatticeType latticeType>
class Equilibrium<T, latticeType, EquilibriumType::Incompressible> {
 public:
  LBM_DEVICE LBM_HOST static inline T calculate(
      const T& density,
      const MathVector<T, L::dimD>& velocity,
      const T& velocity2,
      const unsigned int iQ) {
    {
        LBM_SCOREP_INSTRUMENT_OFF("Equilibrium<T, latticeType, "
                                  "EquilibriumType::Incompressible>::calculate",
                                  5)}

    T cu = L::celerity()[iQ].dot(velocity);

    T fEq_iQ = 1.0 + cu * L::inv_cs2 - 0.5 * velocity2 * L::inv_cs2 +
               0.5 * Power<T, 2>::Do(L::inv_cs2) * cu * cu -
               0.5 * Power<T, 2>::Do(L::inv_cs2) * cu * velocity2 +
               Power<T, 3>::Do(cu) * Power<T, 3>::Do(L::inv_cs2) / 6.0 +
               0.125 * velocity2 * velocity2 * Power<T, 2>::Do(L::inv_cs2) -
               0.25 * cu * cu * velocity2 * Power<T, 3>::Do(L::inv_cs2) +
               Power<T, 4>::Do(cu) * Power<T, 4>::Do(L::inv_cs2) / 24.0;

    return density * L::weight()[iQ] * fEq_iQ;
  }
};

template <class T>
class Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>
    : public Equilibrium<T,
                         LatticeType::Generic,
                         EquilibriumType::Incompressible> {
 public:
  LBM_DEVICE LBM_HOST static inline T calculate(
      const T& density,
      const MathVector<T, L::dimD>& velocity,
      const T& velocity2,
      const unsigned int iQ) {
    {LBM_SCOREP_INSTRUMENT_OFF(
        "Equilibrium<T, D1Q3, EquilibriumType::Incompressible>::calculate", 5)}

    T fEq_iQ = 1.0;

    for (auto iD = 0; iD < L::dimD; ++iD) {
      fEq_iQ *= (2.0 - sqrt(1.0 + 3.0 * velocity[iD] * velocity[iD])) *
                PowerBase((2 * velocity[iD] +
                           sqrt(1.0 + 3.0 * velocity[iD] * velocity[iD])) /
                              (1.0 - velocity[iD]),
                          L::celerity()[iQ][iD]);
    }

    return density * L::weight()[iQ] * fEq_iQ;
  }
};

template <class T>
class Equilibrium<T, LatticeType::D2Q9, EquilibriumType::Incompressible>
    : public Equilibrium<T,
                         LatticeType::D1Q3,
                         EquilibriumType::Incompressible> {
 public:
  LBM_DEVICE LBM_HOST static inline T calculate(
      const T& density,
      const MathVector<T, L::dimD>& velocity,
      const T& velocity2,
      const unsigned int iQ) {
    {LBM_SCOREP_INSTRUMENT_OFF(
        "Equilibrium<T, D2Q9, EquilibriumType::Incompressible>::calculate", 5)}

    T fEq_iQ = 1.0;

    for (auto iD = 0; iD < L::dimD; ++iD) {
      fEq_iQ *= (2.0 - sqrt(1.0 + 3.0 * velocity[iD] * velocity[iD])) *
                PowerBase((2 * velocity[iD] +
                           sqrt(1.0 + 3.0 * velocity[iD] * velocity[iD])) /
                              (1.0 - velocity[iD]),
                          uiL::celerity()[iQ][iD]);
    }

    return density * L::weight()[iQ] * fEq_iQ;
  }
};

template <class T>
class Equilibrium<T, LatticeType::D3Q27, EquilibriumType::Incompressible>
    : public Equilibrium<T,
                         LatticeType::D1Q3,
                         EquilibriumType::Incompressible> {
 public:
  LBM_DEVICE LBM_HOST static inline T calculate(
      const T& density,
      const MathVector<T, L::dimD>& velocity,
      const T& velocity2,
      const unsigned int iQ) {
    LBM_SCOREP_INSTRUMENT_OFF(
        "Equilibrium<T, D3Q27, EquilibriumType::Incompressible>::calculate", 5)

    T fEq_iQ = (T)1;
    for (auto iD = 0; iD < L::dimD; ++iD) {
      fEq_iQ *= (2.0 - sqrt(1.0 + 3.0 * velocity[iD] * velocity[iD])) *
                PowerBase((2 * velocity[iD] +
                           sqrt(1.0 + 3.0 * velocity[iD] * velocity[iD])) /
                              (1.0 - velocity[iD]),
                          L::celerity()[iQ][iD]);
    }

    return density * L::weight()[iQ] * fEq_iQ;
  }
};

typedef Equilibrium<dataT, latticeT, equilibriumT> Equilibrium_;

}  // namespace lbm
