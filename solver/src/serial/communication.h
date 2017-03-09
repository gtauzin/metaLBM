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
    for (int iX = hX<T, L>()+1; iX < hX<T, L>()+lX_l<T, L>()-1; ++iX) {
      for (int iZ = hZ<T, L>()+1; iZ < hZ<T, L>()+lZ_l<T, L>()-1; ++iZ) {

        int srcBtoT, destBtoT;

        srcBtoT  = idxL<T, L>(iX, 0, iZ);
        destBtoT = idxL<T, L>(iX, hY<T, L>()+lY_l<T, L>()-1, iZ);

        UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
            l.f_distribution[idxPop<T, L>(destBtoT, iQ)] = l.f_distribution[idxPop<T, L>(srcBtoT, iQ)];
          });


        srcBtoT  = idxL<T, L>(iX, 2*hY<T, L>()+lY_l<T, L>()-1, iZ);
        destBtoT = idxL<T, L>(iX, hY<T, L>(), iZ);

        UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
            l.f_distribution[idxPop<T, L>(destBtoT, 1)] = l.f_distribution[idxPop<T, L>(srcBtoT, 1)];
          });

      }
    }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = hY<T, L>()+1; iY < hY<T, L>()+lY_l<T, L>()-1; ++iY) {
      for (int iZ = hZ<T, L>()+1; iZ < hZ<T, L>()+lZ_l<T, L>()-1; ++iZ) {

        int srcBtoT, destBtoT;

        srcBtoT  = idxL<T, L>(0, iY, iZ);
        destBtoT = idxL<T, L>(hX<T, L>()+lX_l<T, L>()-1, iY, iZ);

        UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
            l.f_distribution[idxPop<T, L>(destBtoT, iQ)] = l.f_distribution[idxPop<T, L>(srcBtoT, iQ)];
          });

        srcBtoT  = idxL<T, L>(2*hX<T, L>()+lX_l<T, L>()-1, iY, iZ);
        destBtoT = idxL<T, L>(hX<T, L>(), iY, iZ);

        UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
            l.f_distribution[idxPop<T, L>(destBtoT, 1)] = l.f_distribution[idxPop<T, L>(srcBtoT, 1)];
          });

      }
    }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = hY<T, L>()+1; iY < hY<T, L>()+lY_l<T, L>()-1; ++iY) {
      for (int iZ = hZ<T, L>()+1; iZ < hZ<T, L>()+lZ_l<T, L>()-1; ++iZ) {
        int srcBtoT, destBtoT;

        srcBtoT  = idxL<T, L>(0, iY, iZ);
        destBtoT = idxL<T, L>(hX<T, L>()+lX_l<T, L>()-1, iY, iZ);

        UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
            l.f_distribution[idxPop<T, L>(destBtoT, iQ)] = l.f_distribution[idxPop<T, L>(srcBtoT, iQ)];
          });

        srcBtoT  = idxL<T, L>(2*hX<T, L>()+lX_l<T, L>()-1, iY, iZ);
        destBtoT = idxL<T, L>(hX<T, L>(), iY, iZ);

        UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
            l.f_distribution[idxPop<T, L>(destBtoT, 1)] = l.f_distribution[idxPop<T, L>(srcBtoT, 1)];
          });

      }
    }

  }

}

#endif // COMMUNICATION_H
