#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <omp.h>
#include <array>

#include "input.h"
#include "structure.h"
#include "commons.h"
#include "lattice.h"

namespace lbm {

  template<class T>
    inline void pull_handleHalos(Distribution<T>& f) {
#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iX = 0; iX < 2*L::hX + L::lX_l; ++iX) {
      for (int iZ = 0; iZ < 2*L::hZ + L::lZ_l; ++iZ) {

        int srcBtoT, destBtoT;

        srcBtoT  = idxL(iX, L::hY, iZ);
        destBtoT = idxL(iX, L::hY+L::lY_l, iZ);

        UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
            f[idxPop(destBtoT, iQ)] = f[idxPop(srcBtoT, iQ)];
          });


        srcBtoT  = idxL(iX, L::hY+L::lY_l-1, iZ);
        destBtoT = idxL(iX, 0, iZ);

        UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
            f[idxPop(destBtoT, iQ)] = f[idxPop(srcBtoT, iQ)];
          });

      }
    }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = 0; iY < 2*L::hY+L::lY_l; ++iY) {
      for (int iZ = 0; iZ < 2*L::hZ+L::lZ_l; ++iZ) {

        int srcBtoT, destBtoT;

        srcBtoT  = idxL(L::hX, iY, iZ);
        destBtoT = idxL(L::hX + L::lX_l, iY, iZ);

        UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
            f[idxPop(destBtoT, iQ)] = f[idxPop(srcBtoT, iQ)];
          });

        srcBtoT  = idxL(L::hX + L::lX_l-1, iY, iZ);
        destBtoT = idxL(0, iY, iZ);

        UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
            f[idxPop(destBtoT, iQ)] = f[idxPop(srcBtoT, iQ)];
          });

      }
    }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = 0; iY < 2*L::hY + L::lY_l; ++iY) {
      for (int iX = 0; iX < 2*L::hX + L::lX_l; ++iX) {
        int srcBtoT, destBtoT;

        srcBtoT  = idxL(iX, iY, L::hZ);
        destBtoT = idxL(iX, iY, L::hZ + L::lZ_l);

        UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
            f[idxPop(destBtoT, iQ)] = f[idxPop(srcBtoT, iQ)];
          });

        srcBtoT  = idxL(iX, iY, L::hZ + L::lZ_l-1);
        destBtoT = idxL(iX, iY, 0);

        UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
            f[idxPop(destBtoT, iQ)] = f[idxPop(srcBtoT, iQ)];
          });

      }
    }

  }

}

#endif // COMMUNICATION_H
