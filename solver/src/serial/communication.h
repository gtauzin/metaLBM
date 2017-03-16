#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <omp.h>
#include <array>

#include "input.h"
#include "structure.h"
#include "commons.h"
#include "lattice.h"

namespace lbm {

  template<class T, LatticeType L>
    inline void push_handleHalos(Lattice<T, L>& l) {
#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iX = P::hX+1; iX < P::hX + P::lX_l-1; ++iX) {
      for (int iZ = P::hZ+1; iZ < P::hZ + P::lZ_l-1; ++iZ) {

        int srcBtoT, destBtoT;

        srcBtoT  = idxL(iX, 0, iZ);
        destBtoT = idxL(iX, P::hY+P::lY_l-1, iZ);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, iQ)] = l.f_distribution[idxPop(srcBtoT, iQ)];
          });


        srcBtoT  = idxL(iX, 2*P::hY+P::lY_l-1, iZ);
        destBtoT = idxL(iX, P::hY, iZ);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, 1)] = l.f_distribution[idxPop(srcBtoT, 1)];
          });

      }
    }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = P::hY+1; iY < P::hY+P::lY_l-1; ++iY) {
      for (int iZ = P::hZ+1; iZ < P::hZ+P::lZ_l-1; ++iZ) {

        int srcBtoT, destBtoT;

        srcBtoT  = idxL(0, iY, iZ);
        destBtoT = idxL(P::hX+P::lX_l-1, iY, iZ);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, iQ)] = l.f_distribution[idxPop(srcBtoT, iQ)];
          });

        srcBtoT  = idxL(2*P::hX+P::lX_l-1, iY, iZ);
        destBtoT = idxL(P::hX, iY, iZ);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, 1)] = l.f_distribution[idxPop(srcBtoT, 1)];
          });

      }
    }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = P::hY+1; iY < P::hY + P::lY_l-1; ++iY) {
      for (int iZ = P::hZ+1; iZ < P::hZ + P::lZ_l-1; ++iZ) {
        int srcBtoT, destBtoT;

        srcBtoT  = idxL(0, iY, iZ);
        destBtoT = idxL(P::hX+P::lX_l-1, iY, iZ);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, iQ)] = l.f_distribution[idxPop(srcBtoT, iQ)];
          });

        srcBtoT  = idxL(2*P::hX+P::lX_l-1, iY, iZ);
        destBtoT = idxL(P::hX, iY, iZ);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, 1)] = l.f_distribution[idxPop(srcBtoT, 1)];
          });

      }
    }

  }

}

#endif // COMMUNICATION_H
