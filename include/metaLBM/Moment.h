#pragma once

#include "Commons.h"
#include "Helpers.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Options.h"

namespace lbm {
template <class T>
class Moment {
 private:
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
                                T& density,
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
  static void calculateT(const T* haloDistributionPreviousPtr, const T* haloDistributionNextPtr,
                         const Position& iP, T& T2, T& T3, T& T4) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateVelocity", 5)

    T non_equilibrium_0 = haloDistributionNextPtr[hSD::getIndex(iP, 0)];
    T equilibrium_0 =
      haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[0], 0)]
      - haloDistributionNextPtr[hSD::getIndex(iP, 0)];

    T2 = non_equilibrium_0 * non_equilibrium_0 / equilibrium_0;
    T3 = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0 / equilibrium_0;
    T4 = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0 / equilibrium_0 /equilibrium_0;

    #pragma unroll
    for (auto iQ = 1; iQ < L::dimQ; ++iQ) {
      T non_equilibrium_iQ = haloDistributionNextPtr[hSD::getIndex(iP, iQ)];
      T equilibrium_iQ =
        haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)]
        - haloDistributionNextPtr[hSD::getIndex(iP, iQ)];

      T2 += non_equilibrium_iQ * non_equilibrium_iQ / equilibrium_iQ;
      T3 += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ / equilibrium_iQ;
      T4 += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ / equilibrium_iQ /equilibrium_iQ;
    }
  }


  LBM_DEVICE LBM_HOST LBM_INLINE
  static void calculateT_forcing(const T* haloDistributionPreviousPtr, const T* haloDistributionNextPtr,
                                 const Position& iP, T& T2, T& T3, T& T4) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateVelocity", 5)

    T non_equilibrium_0 = haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[0], 0)];
    T equilibrium_0 = haloDistributionNextPtr[hSD::getIndex(iP, 0)] - non_equilibrium_0;


    T2 = non_equilibrium_0 * non_equilibrium_0 / equilibrium_0;
    T3 = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0 / equilibrium_0;
    T4 = non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0 * non_equilibrium_0
      / equilibrium_0 / equilibrium_0 /equilibrium_0;

    #pragma unroll
    for (auto iQ = 1; iQ < L::dimQ; ++iQ) {
    T non_equilibrium_iQ = haloDistributionPreviousPtr[hSD::getIndex(iP - uiL::celerity()[iQ], iQ)];
    T equilibrium_iQ = haloDistributionNextPtr[hSD::getIndex(iP, iQ)] - non_equilibrium_iQ;

      T2 += non_equilibrium_iQ * non_equilibrium_iQ / equilibrium_iQ;
      T3 += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ / equilibrium_iQ;
      T4 += non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ * non_equilibrium_iQ
        / equilibrium_iQ / equilibrium_iQ /equilibrium_iQ;
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
};

typedef Moment<dataT> Moment_;

}  // namespace lbm
