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
    inline void pull_handleHalos(Lattice<T, L>& l) {
#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iX = 0; iX < 2*P::hX + P::lX_l; ++iX) {
      for (int iZ = 0; iZ < 2*P::hZ + P::lZ_l; ++iZ) {

        int srcBtoT, destBtoT;

        srcBtoT  = idxL(iX, P::hY, iZ);
        destBtoT = idxL(iX, P::hY+P::lY_l, iZ);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, iQ)] = l.f_distribution[idxPop(srcBtoT, iQ)];
          });


        srcBtoT  = idxL(iX, P::hY+P::lY_l-1, iZ);
        destBtoT = idxL(iX, 0, iZ);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, iQ)] = l.f_distribution[idxPop(srcBtoT, iQ)];
          });

      }
    }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = 0; iY < 2*P::hY+P::lY_l; ++iY) {
      for (int iZ = 0; iZ < 2*P::hZ+P::lZ_l; ++iZ) {

        int srcBtoT, destBtoT;

        srcBtoT  = idxL(P::hX, iY, iZ);
        destBtoT = idxL(P::hX + P::lX_l, iY, iZ);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, iQ)] = l.f_distribution[idxPop(srcBtoT, iQ)];
          });

        srcBtoT  = idxL(P::hX + P::lX_l-1, iY, iZ);
        destBtoT = idxL(0, iY, iZ);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, iQ)] = l.f_distribution[idxPop(srcBtoT, iQ)];
          });

      }
    }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = 0; iY < 2*P::hY + P::lY_l; ++iY) {
      for (int iX = 0; iX < 2*P::hX + P::lX_l; ++iX) {
        int srcBtoT, destBtoT;

        srcBtoT  = idxL(iX, iY, P::hZ);
        destBtoT = idxL(iX, iY, P::hZ + P::lZ_l);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, iQ)] = l.f_distribution[idxPop(srcBtoT, iQ)];
          });

        srcBtoT  = idxL(iX, iY, P::hZ + P::lZ_l-1);
        destBtoT = idxL(iX, iY, 0);

        UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
            l.f_distribution[idxPop(destBtoT, iQ)] = l.f_distribution[idxPop(srcBtoT, iQ)];
          });

      }
    }

  }

}

#endif // COMMUNICATION_H
