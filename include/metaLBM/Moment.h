#ifndef MOMENT_H
#define MOMENT_H

#include "Commons.h"
#include "Options.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Helpers.h"

namespace lbm {

  template <class T>
  class Moment {
  private:
  public:

    DEVICE HOST
    static inline void calculateDensity(const T * haloDistributionPtr,
                                        const Position& iP,
                                        T& density) {
      { INSTRUMENT_OFF("Moment<T>::calculateDensity",5)}

      density = haloDistributionPtr[hSD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)];

      for(auto iQ = 1; iQ < L::dimQ; ++iQ) {
        density += haloDistributionPtr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)];
      }
    }

    DEVICE HOST
    static inline void calculateVelocity(const T * haloDistributionPtr,
                                  const Position& iP,
                                  T& density,
                                  MathVector<T, L::dimD>& velocity) {
      { INSTRUMENT_OFF("Moment<T>::calculateVelocity",5) }

       velocity = L::celerity()[0]
        * haloDistributionPtr[hSD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)];

       for(auto iQ = 1; iQ < L::dimQ; ++iQ) {
         velocity += L::celerity()[iQ]
           * haloDistributionPtr[hSD::getIndex(iP-uiL::celerity()[iQ], iQ)];
       }
      velocity /= density;
    }

    DEVICE HOST
    static inline T calculateEntropy(const T * haloDistributionPtr,
                              const Position& iP,
                              T& entropy) {
      { INSTRUMENT_OFF("Moment<T>::calculateEntropy",5) }

      entropy = haloDistributionPtr[hSD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)]
        * log(haloDistributionPtr[hSD::getIndex(iP-uiL::celerity()[0], (unsigned int) 0)]
              /L::weight()[0]);

      for(auto iQ = 1; iQ < L::dimQ; ++iQ) {
        int indexPop_iQ = hSD::getIndex(iP-uiL::celerity()[iQ], iQ);
        entropy += haloDistributionPtr[indexPop_iQ]
          * log(haloDistributionPtr[indexPop_iQ]/L::weight()[iQ]);
      }
    }

  };

  typedef Moment<dataT> Moment_;

}

#endif // MOMENT_H
