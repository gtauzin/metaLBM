#ifndef COMMONS_H
#define COMMONS_H

#include <iostream>
#include <string>
#include <stdlib.h>

#include "input.h"
#include "lattice.h"
#include "structure.h"

#ifdef __CUDACC__
#define _HOST __host__
#define _DEVICE __device__
#else
#define _HOST
#define _DEVICE
#endif

#define CACHE_LINE 64

namespace lbm {

  enum d {X, Y, Z};

  constexpr int sizeX_l() { return 2*L::hX + L::lX_l; }

  constexpr int sizeY_l() { return 2*L::hY + L::lY_l; }

  constexpr int sizeZ_l() { return 2*L::hZ + L::lZ_l; }

  constexpr int s_l() { return L::lX_l * L::lY_l * L::lZ_l; }

  constexpr int s_g() { return L::lX_g * L::lY_g * L::lZ_g; }



  _HOST _DEVICE inline int idxL(const int iX, const int iY, const int iZ) {
    return sizeZ_l()*(sizeY_l()*iX + iY) + iZ;
  }

  _HOST _DEVICE inline int idxL(const MathVector<int, 3>& iP) {
    return sizeZ_l()*(sizeY_l()*iP[d::X] + iP[d::Y]) + iP[d::Z];
  }

  _HOST _DEVICE inline int idx_inL(const int iX, const int iY, const int iZ) {
    return sizeZ_l()*(sizeY_l()*(iX-L::hX) + (iY-L::hY)) + (iZ-L::hZ);
  }

  _HOST _DEVICE inline int idx_outL(const int iX, const int iY, const int iZ) {
    return sizeZ_l()*(sizeY_l()*(iX+L::hX) + (iY+L::hY)) + (iZ+L::hZ);
  }

  _HOST _DEVICE inline int idx_inF(const int iX, const int iY, const int iZ) {
    return L::lZ_l * (L::lY_l * (iX-L::hX) + (iY-L::hY)) + (iZ-L::hZ);
  }

  _HOST _DEVICE inline int idx_lF(const int iX, const int iY, const int iZ) {
    return L::lZ_l * (L::lY_l*iX + iY) + iZ;
  }

  _HOST _DEVICE inline int idx_gF(const int iX, const int iY, const int iZ) {
    const int iX_l = iX - iX/L::lX_l * L::lX_l;
    return iX/L::lX_l * s_l() + idx_lF(iX_l, iY, iZ);
  }

#if defined(_AOS)

  _HOST _DEVICE inline int idxPop(const int idx, const int iQ) {
    return idx * L::dimQ + iQ;
  }

  _HOST _DEVICE inline int idxPop(const int iX, const int iY, const int iZ,
                                  const int iQ) {
    return idxL(iX, iY, iZ) * L::dimQ + iQ;
  }

  _HOST _DEVICE constexpr int idxPop_lF(const int idx, const int iQ) {
    return idx * L::dimQ + iQ;
  }

  _HOST _DEVICE inline int idxPop_gF(const int idx, const int iQ) {
    return idx * L::dimQ + iQ;
  }

  _HOST _DEVICE inline int idxPop_gF(const int iX, const int iY, const int iZ,
                                     const int iQ) {
    return idx_gF(iX, iY, iZ) * L::dimQ + iQ;
  }

#elif defined(_SOA)

  _HOST _DEVICE inline int idxPop(const int idx, const int iQ) {
    return iQ * sizeX_l() * sizeY_l() * sizeZ_l()  + idx;
  }

  _HOST _DEVICE inline int idxPop(const int iX, const int iY, const int iZ,
                                    const int iQ) {
    return iQ * sizeX_l() * sizeY_l() * sizeZ_l() + idxL(iX, iY, iZ);
  }

  _HOST _DEVICE inline int idxPop_lF(const int idx, const int iQ) {
    return iQ * L::lX_l * L::lY_l * L::lZ_l + idx;
  }

  _HOST _DEVICE inline int idxPop_lF(const int iX, const int iY, const int iZ,
                                     const int iQ) {
    return iQ * L::lX_g * L::lY_g * L::lZ_g + idx_lF(iX, iY, iZ);
  }

  _HOST _DEVICE inline int idxPop_gF(const int idx, const int iQ) {
    return iQ * L::lX_g * L::lY_g * L::lZ_g + idx;
  }

  _HOST _DEVICE inline int idxPop_gF(const int iX, const int iY, const int iZ,
                                     const int iQ) {
    const int iX_l = iX - iX/L::lX_l * L::lX_l;
    return iX/L::lX_l * L::dimQ * s_l() + idxPop_lF(iX_l, iY, iZ, iQ);
  }


#else
# error "Error: data structure format is not defined."
#endif

}

#endif // COMMONS_H
