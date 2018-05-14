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
  LBM_DEVICE LBM_HOST LBM_INLINE static void calculateDensity(
      const T* haloDistributionPtr,
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

  LBM_DEVICE LBM_HOST LBM_INLINE static void calculateVelocity(
      const T* haloDistributionPtr,
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

  LBM_DEVICE LBM_HOST static inline T calculateEntropy(
      const T* haloDistributionPtr,
      const Position& iP,
      T& entropy) {
    LBM_INSTRUMENT_OFF("Moment<T>::calculateEntropy", 5)

    entropy = haloDistributionPtr[hSD::getIndex(iP - uiL::celerity()[0],
                                                (unsigned int)0)] *
              log(haloDistributionPtr[hSD::getIndex(iP - uiL::celerity()[0],
                                                    (unsigned int)0)] /
                  L::weight()[0]);

    for (auto iQ = 1; iQ < L::dimQ; ++iQ) {
      int indexPop_iQ = hSD::getIndex(iP - uiL::celerity()[iQ], iQ);
      entropy += haloDistributionPtr[indexPop_iQ] *
                 log(haloDistributionPtr[indexPop_iQ] / L::weight()[iQ]);
    }
  }
};

typedef Moment<dataT> Moment_;

}  // namespace lbm
