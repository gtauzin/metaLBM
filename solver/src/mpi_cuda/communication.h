#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include "input.h"
#include "commons.h"
#include "lattice.h"

#include <omp.h>
#include <array>

namespace lbm {

 inline void push_handleHalos(Lattice& l) {
    #pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iX = haloX+1; iX < haloX+lengthX_l-1; ++iX) {
      int srcBtoT, destBtoT;

      srcBtoT  = idxL(iX, 0);
      destBtoT = idxL(iX, haloY+lengthY_l-1);
      l.f_distribution[idxPop(destBtoT, 3)] = l.f_distribution[idxPop(srcBtoT, 3)];
      l.f_distribution[idxPop(destBtoT, 4)] = l.f_distribution[idxPop(srcBtoT, 4)];
      l.f_distribution[idxPop(destBtoT, 5)] = l.f_distribution[idxPop(srcBtoT, 5)];

      srcBtoT  = idxL(iX, 2*haloY+lengthY_l-1);
      destBtoT = idxL(iX, haloY);
      l.f_distribution[idxPop(destBtoT, 1)] = l.f_distribution[idxPop(srcBtoT, 1)];
      l.f_distribution[idxPop(destBtoT, 8)] = l.f_distribution[idxPop(srcBtoT, 8)];
      l.f_distribution[idxPop(destBtoT, 7)] = l.f_distribution[idxPop(srcBtoT, 7)];
    }

    int srcBtoT, destBtoT;

    srcBtoT = idxL(haloX, 0);
    destBtoT = idxL(haloX, haloY+lengthY_l-1);
    l.f_distribution[idxPop(destBtoT, 3)] = l.f_distribution[idxPop(srcBtoT, 3)];
    l.f_distribution[idxPop(destBtoT, 4)] = l.f_distribution[idxPop(srcBtoT, 4)];

    srcBtoT = idxL(haloX+lengthX_l-1, 0);
    destBtoT = idxL(haloX+lengthX_l-1, haloY+lengthY_l-1);
    l.f_distribution[idxPop(destBtoT, 4)] = l.f_distribution[idxPop(srcBtoT, 4)];
    l.f_distribution[idxPop(destBtoT, 5)] = l.f_distribution[idxPop(srcBtoT, 5)];

    srcBtoT = idxL(haloX, 2*haloY+lengthY_l-1);
    destBtoT = idxL(haloX, haloY);
    l.f_distribution[idxPop(destBtoT, 1)] = l.f_distribution[idxPop(srcBtoT, 1)];
    l.f_distribution[idxPop(destBtoT, 8)] = l.f_distribution[idxPop(srcBtoT, 8)];

    srcBtoT = idxL(haloX+lengthX_l-1, 2*haloY+lengthY_l-1);
    destBtoT = idxL(haloX+lengthX_l-1, haloY);
    l.f_distribution[idxPop(destBtoT, 7)] = l.f_distribution[idxPop(srcBtoT, 7)];
    l.f_distribution[idxPop(destBtoT, 8)] = l.f_distribution[idxPop(srcBtoT, 8)];
  }

inline void push_packToSend_mpi(Lattice& l,
                                std::array<double, size_buf>& sendLeft_buffer,
                                std::array<double, size_buf>& sendRight_buffer) {
#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = haloY+1; iY < haloY+lengthY_l-1; ++iY) {
      int idx_left = idxL(0, iY);
      sendLeft_buffer[(iY-haloY)*3 + 0] = l.f_distribution[idxPop(idx_left, 1)];
      sendLeft_buffer[(iY-haloY)*3 + 1] = l.f_distribution[idxPop(idx_left, 2)];
      sendLeft_buffer[(iY-haloY)*3 + 2] = l.f_distribution[idxPop(idx_left, 3)];
    }

    int iY, idx_left;

    iY = haloY;
    idx_left = idxL(0, haloY);
    sendLeft_buffer[(iY-haloY)*3 + 1] = l.f_distribution[idxPop(idx_left, 2)];
    sendLeft_buffer[(iY-haloY)*3 + 2] = l.f_distribution[idxPop(idx_left, 3)];

    idx_left = idxL(0, 2*haloY+lengthY_l-1);
    sendLeft_buffer[(iY-haloY)*3 + 0] = l.f_distribution[idxPop(idx_left, 1)];

    iY = haloY+lengthY_l-1;
    idx_left = idxL(0, haloY+lengthY_l-1);
    sendLeft_buffer[(iY-haloY)*3 + 0] = l.f_distribution[idxPop(idx_left, 1)];
    sendLeft_buffer[(iY-haloY)*3 + 1] = l.f_distribution[idxPop(idx_left, 2)];

    idx_left = idxL(0, 0);
    sendLeft_buffer[(iY-haloY)*3 + 2] = l.f_distribution[idxPop(idx_left, 3)];

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = haloY+1; iY < haloY+lengthY_l-1; ++iY) {
      int idx_right = idxL(2*haloX+lengthX_l-1, iY);
      sendRight_buffer[(iY-haloY)*3 + 0] = l.f_distribution[idxPop(idx_right, 7)];
      sendRight_buffer[(iY-haloY)*3 + 1] = l.f_distribution[idxPop(idx_right, 6)];
      sendRight_buffer[(iY-haloY)*3 + 2] = l.f_distribution[idxPop(idx_right, 5)];
    }

    int idx_right;

    iY = haloY;
    idx_right = idxL(2*haloX+lengthX_l-1, haloY);
    sendRight_buffer[(iY-haloY)*3 + 1] = l.f_distribution[idxPop(idx_right, 6)];
    sendRight_buffer[(iY-haloY)*3 + 2] = l.f_distribution[idxPop(idx_right, 5)];

    idx_right = idxL(2*haloX+lengthX_l-1, 2*haloY+lengthY_l-1);
    sendRight_buffer[(iY-haloY)*3 + 0] = l.f_distribution[idxPop(idx_right, 7)];

    iY = haloY+lengthY_l-1;
    idx_right = idxL(2*haloX+lengthX_l-1, haloY+lengthY_l-1);
    sendRight_buffer[(iY-haloY)*3 + 0] = l.f_distribution[idxPop(idx_right, 7)];
    sendRight_buffer[(iY-haloY)*3 + 1] = l.f_distribution[idxPop(idx_right, 6)];

    idx_right = idxL(2*haloX+lengthX_l-1, 0);
    sendRight_buffer[(iY-haloY)*3 + 2] = l.f_distribution[idxPop(idx_right, 5)];
  }

  inline void push_unpackReceived_mpi(Lattice& l,
                                      const std::array<double, size_buf>& receiveLeft_buffer,
                                      const std::array<double, size_buf>& receiveRight_buffer) {
#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = haloY; iY < haloY+lengthY_l; ++iY) {
      int idx_right = idxL(haloX+lengthX_l-1, iY);
      l.f_distribution[idxPop(idx_right, 1)] = receiveRight_buffer[(iY-haloY)*3 + 0];
      l.f_distribution[idxPop(idx_right, 2)] = receiveRight_buffer[(iY-haloY)*3 + 1];
      l.f_distribution[idxPop(idx_right, 3)] = receiveRight_buffer[(iY-haloY)*3 + 2];
    }

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for (int iY = haloY; iY < haloY+lengthY_l; ++iY) {
      int idx_left = idxL(haloX, iY);
      l.f_distribution[idxPop(idx_left, 7)] = receiveLeft_buffer[(iY-haloY)*3 + 0];
      l.f_distribution[idxPop(idx_left, 6)] = receiveLeft_buffer[(iY-haloY)*3 + 1];
      l.f_distribution[idxPop(idx_left, 5)] = receiveLeft_buffer[(iY-haloY)*3 + 2];
    }
  }

}

#endif // COMMUNICATION_H
