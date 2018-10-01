#pragma once

#include "Commons.h"
#include "Helpers.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Options.h"

namespace lbm {
template <class T>
class Moment {
 public:
  LBM_DEVICE LBM_HOST LBM_INLINE
  static void calculateDensity(const T* haloDistributionPtr,
                               const Position& iP,
                               T& density) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateDensity", 5)

    density = haloDistributionPtr[hSD::getIndex(iP - uiL::celerity()[0],
                                                (unsigned int)0)];

    #pragma unroll
    for (auto iQ = 1; iQ < L::dimQ; ++iQ) {
      density +=
          haloDistributionPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)];
    }
  }

  LBM_DEVICE LBM_HOST LBM_INLINE
  static void calculateVelocity(const T* haloDistributionPtr,
                                const Position& iP,
                                const T density,
                                MathVector<T, L::dimD>& velocity) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateVelocity", 5)

    velocity =
        L::celerity()[0] * haloDistributionPtr[hSD::getIndex(
                               iP - uiL::celerity()[0], (unsigned int)0)];

    #pragma unroll
    for (auto iQ = 1; iQ < L::dimQ; ++iQ) {
      velocity +=
          L::celerity()[iQ] *
          haloDistributionPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)];
    }
    velocity /= density;
  }


  LBM_DEVICE LBM_HOST LBM_INLINE
  static void calculateObservables(const T* haloDistributionPreviousPtr, const T* haloDistributionNextPtr,
                                   const T density, const Position& iP, T& T2, T& T3, T& T4,
                                   T& T2_approx, T& T3_approx, T& T4_approx,
                                   const MathVector<MathVector<dataT, L::dimD>, L::dimQ>& qDiagonal,
                                   const MathVector<MathVector<dataT, 2*L::dimD-3>, L::dimQ>& qSymmetric,
                                   MathVector<T, L::dimD>& pi1Diagonal,
                                   MathVector<T, 2*L::dimD-3>& pi1Symmetric,
                                   T& squaredQContractedPi1, T& cubedQContractedPi1) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateObservables", 5)

    calculatePi1Diagonal(haloDistributionNextPtr, iP, pi1Diagonal);
    calculatePi1Symmetric(haloDistributionNextPtr, iP, pi1Symmetric);
    calculatePowerQContractedPi1(qDiagonal, qSymmetric, pi1Diagonal, pi1Symmetric,
                                 squaredQContractedPi1, cubedQContractedPi1);


    T non_equilibrium_0 = haloDistributionNextPtr[hSD::getIndex(iP, 0)];
    T equilibrium_0 =
      haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[0], 0)]
      - haloDistributionNextPtr[hSD::getIndex(iP, 0)];
    T equilibrium_0_approx = density * L::weight()[0];

    T2 = non_equilibrium_0 * non_equilibrium_0 / equilibrium_0;
    T3 = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0 / equilibrium_0;
    T4 = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0 / equilibrium_0 /equilibrium_0;

    T2_approx = non_equilibrium_0 * non_equilibrium_0 / equilibrium_0_approx;
    T3_approx = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0_approx / equilibrium_0_approx;
    T4_approx = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0_approx / equilibrium_0_approx /equilibrium_0_approx;


    #pragma unroll
    for (auto iQ = 1; iQ < L::dimQ; ++iQ) {
      T non_equilibrium_iQ = haloDistributionNextPtr[hSD::getIndex(iP, iQ)];
      T equilibrium_iQ =
        haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
        - haloDistributionNextPtr[hSD::getIndex(iP, iQ)];
      T equilibrium_iQ_approx = density * L::weight()[iQ];

      T2 += non_equilibrium_iQ * non_equilibrium_iQ / equilibrium_iQ;
      T3 += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ / equilibrium_iQ;
      T4 += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ / equilibrium_iQ /equilibrium_iQ;

      T2_approx += non_equilibrium_iQ * non_equilibrium_iQ / equilibrium_iQ_approx;
      T3_approx += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ_approx / equilibrium_iQ_approx;
      T4_approx += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ_approx / equilibrium_iQ_approx /equilibrium_iQ_approx;
    }
  }


  LBM_DEVICE LBM_HOST LBM_INLINE
  static void calculateObservables_forcing(const T* haloDistributionPreviousPtr, const T* haloDistributionNextPtr,
                                           const T density, const Position& iP, T& T2, T& T3, T& T4,
                                           T& T2_approx, T& T3_approx, T& T4_approx) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateVelocity", 5)

    T non_equilibrium_0 = haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[0], 0)];
    T equilibrium_0 = haloDistributionNextPtr[hSD::getIndex(iP, 0)] - non_equilibrium_0;
    T equilibrium_0_approx = density * L::weight()[0];


    T2 = non_equilibrium_0 * non_equilibrium_0 / equilibrium_0;
    T3 = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0 / equilibrium_0;
    T4 = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0 / equilibrium_0 /equilibrium_0;

    T2_approx = non_equilibrium_0 * non_equilibrium_0 / equilibrium_0_approx;
    T3_approx = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0_approx / equilibrium_0_approx;
    T4_approx = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0_approx / equilibrium_0_approx /equilibrium_0_approx;

    #pragma unroll
    for (auto iQ = 1; iQ < L::dimQ; ++iQ) {
    T non_equilibrium_iQ = haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)];
    T equilibrium_iQ = haloDistributionNextPtr[hSD::getIndex(iP, iQ)] - non_equilibrium_iQ;
    T equilibrium_iQ_approx = density * L::weight()[iQ];

      T2 += non_equilibrium_iQ * non_equilibrium_iQ / equilibrium_iQ;
      T3 += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ / equilibrium_iQ;
      T4 += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ / equilibrium_iQ /equilibrium_iQ;

      T2_approx += non_equilibrium_iQ * non_equilibrium_iQ / equilibrium_iQ_approx;
      T3_approx += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ_approx / equilibrium_iQ_approx;
      T4_approx += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ_approx / equilibrium_iQ_approx /equilibrium_iQ_approx;
    }
  }


  LBM_DEVICE LBM_HOST
  static inline T calculateEntropy(const T* haloDistributionPtr,
                                   const Position& iP,
                                   T& entropy) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateEntropy", 5)

    entropy = haloDistributionPtr[hSD::getIndex(iP - uiL::celerity()[0],
                                                (unsigned int)0)] *
              log(haloDistributionPtr[hSD::getIndex(iP - uiL::celerity()[0],
                                                    (unsigned int)0)] /
                  L::weight()[0]);

    #pragma unroll
    for (auto iQ = 1; iQ < L::dimQ; ++iQ) {
      int indexPop_iQ = hSD::getIndex(iP - uiL::celerity()[iQ], iQ);
      entropy += haloDistributionPtr[indexPop_iQ] *
                 log(haloDistributionPtr[indexPop_iQ] / L::weight()[iQ]);
    }
  }


  LBM_DEVICE LBM_HOST inline
  static void calculatePi1Diagonal(const T* haloDistributionNextPtr, const Position& iP,
                                   MathVector<T, L::dimD>& pi1Diagonal) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculatePi1Diagonal", 5)

    T non_equilibrium_0 = haloDistributionNextPtr[hSD::getIndex(iP, 0)];
    for (auto iD = 0; iD < L::dimD; ++iD) {
      pi1Diagonal[iD] = L::celerity()[0][iD] * L::celerity()[0][iD]
        * non_equilibrium_0;
    }

    T non_equilibrium_iQ;
    for (auto iQ = 1; iQ < L::dimQ; ++iQ) {
      non_equilibrium_iQ = haloDistributionNextPtr[hSD::getIndex(iP, iQ)];
      for (auto iD = 0; iD < L::dimD; ++iD) {
        pi1Diagonal[iD] += L::celerity()[iQ][iD] * L::celerity()[iQ][iD]
          * non_equilibrium_iQ;
      }
    }
  }

  LBM_DEVICE LBM_HOST inline
  static void calculatePi1Symmetric(const T* haloDistributionNextPtr, const Position& iP,
                                    MathVector<T, 2*L::dimD-3>& pi1Symmetric) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculatePi1Symmetric", 5)

    T non_equilibrium_0 = haloDistributionNextPtr[hSD::getIndex(iP, 0)];
    pi1Symmetric[d::X] = L::celerity()[0][d::X] * L::celerity()[0][d::Y]
      * non_equilibrium_0;

    for (auto iD = 1; iD < 2 * L::dimD - 3; ++iD) {
      pi1Symmetric[iD] += L::celerity()[0][iD-1] * L::celerity()[0][d::Z]
      * non_equilibrium_0;
    }

    T non_equilibrium_iQ;
    for (auto iQ = 1; iQ < L::dimQ; ++iQ) {
      non_equilibrium_iQ = haloDistributionNextPtr[hSD::getIndex(iP, iQ)];
      pi1Symmetric[d::X] += L::celerity()[iQ][d::X] * L::celerity()[iQ][d::Y]
        * non_equilibrium_iQ;

      for (auto iD = 1; iD < 2 * L::dimD - 3; ++iD) {
        pi1Symmetric[iD] += L::celerity()[iQ][iD-1] * L::celerity()[iQ][d::Z]
        * non_equilibrium_iQ;
      }
    }
  }

  LBM_DEVICE LBM_HOST inline
    static void calculatePowerQContractedPi1( const MathVector<MathVector<dataT, L::dimD>, L::dimQ>& qDiagonal,
                                              const MathVector<MathVector<dataT, 2*L::dimD-3>, L::dimQ>& qSymmetric,
                                              const MathVector<T, L::dimD>& pi1Diagonal,
                                              const MathVector<T, 2*L::dimD-3>& pi1Symmetric,
                                           T& squaredQContractedPi1, T& cubedQContractedPi1) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculatePowerQContractedPi1", 5)

    T qContractedPi1_iQ;

    squaredQContractedPi1 = (T) 0;
    cubedQContractedPi1 = (T) 0;

    for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
      qContractedPi1_iQ = 0;

      for (auto iD = 0; iD < L::dimD; ++iD) {
        qContractedPi1_iQ += qDiagonal[iQ][iD] * pi1Diagonal[iD];

      }
      for (auto iD = 0; iD < 2 * L::dimD - 3; ++iD) {
        qContractedPi1_iQ += 2 * qSymmetric[iQ][iD] * pi1Symmetric[iD];
      }
      squaredQContractedPi1 += qContractedPi1_iQ * qContractedPi1_iQ;
      cubedQContractedPi1 += qContractedPi1_iQ * qContractedPi1_iQ * qContractedPi1_iQ;
    }
  }

 };

 typedef Moment<dataT> Moment_;

}  // namespace lbm
