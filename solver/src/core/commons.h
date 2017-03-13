#ifndef COMMONS_H
#define COMMONS_H

#include <iostream>
#include <string>
#include <stdlib.h>

#include "input.h"
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

  /**
   * Parameters required to define a lattice for a parallel code.
   *
   * @tparam T data type.
   * @tparam L lattice type of the form DdQq.
   */
  template <class T, LatticeType L>
    struct Parameters;


  template <class T>
    struct Parameters<T, LatticeType::D1Q3>
    {
      static constexpr int dimD = 1;
      static constexpr int dimQ = 3;
      static constexpr int lX_g = lengthX_g;
      static constexpr int lY_g = 1;
      static constexpr int lZ_g = 1;

      static constexpr int lX_l = lengthX_l;
      static constexpr int lY_l = 1;
      static constexpr int lZ_l = 1;

      static constexpr int hX = 1;
      static constexpr int hY = 0;
      static constexpr int hZ = 0;


      static constexpr T inv_cs2 = 3.0;
      static constexpr T cs2 = 1./inv_cs2;
      static constexpr MathVector<MathVector<T, dimD>, dimQ> celerity =
        {
          MathVector<T, dimD>{{0.0}},
          MathVector<T, dimD>{{-1.0}},
          MathVector<T, dimD>{{1.0}}
        };
      static constexpr MathVector<T, dimQ> weight =
        {
          2./3., 1./3., 1./3.
        };
    };


  template <class T>
    struct Parameters<T, LatticeType::D2Q9>
    {
      static constexpr int dimD = 2;
      static constexpr int dimQ = 9;
      static constexpr int lX_g = lengthX_g;
      static constexpr int lY_g = lengthY_g;
      static constexpr int lZ_g = 1;

      static constexpr int lX_l = lengthX_l;
      static constexpr int lY_l = lengthY_l;
      static constexpr int lZ_l = 1;

      static constexpr int hX = 1;
      static constexpr int hY = 1;
      static constexpr int hZ = 0;

      static constexpr T inv_cs2 = 3.0;
      static constexpr T cs2 = 1./inv_cs2;
      static constexpr MathVector<MathVector<T, dimD>, dimQ> celerity =
        {
          MathVector<T, dimD>{{0.0, 0.0}},
          MathVector<T, dimD>{{-1.0, 1.0}},
          MathVector<T, dimD>{{-1.0, 0.0}},
          MathVector<T, dimD>{{-1.0, -1.0}},
          MathVector<T, dimD>{{0.0, -1.0}},
          MathVector<T, dimD>{{1.0, -1.0}},
          MathVector<T, dimD>{{1.0, 0.0}},
          MathVector<T, dimD>{{1.0, 1.0}},
          MathVector<T, dimD>{{0.0, 1.0}}
        };

      static constexpr MathVector<T, dimQ> weight =
        {
          4./9., 1./36., 1./9.,
          1./36., 1./9., 1./36.,
          1./9., 1./36., 1./9.
        };

    };

  template <class T>
    struct Parameters<T, LatticeType::D3Q27>
    {
      static constexpr int dimD = 3;
      static constexpr int dimQ = 27;
      static constexpr int lX_g = lengthX_g;
      static constexpr int lY_g = lengthY_g;
      static constexpr int lZ_g = lengthZ_g;

      static constexpr int lX_l = lengthX_l;
      static constexpr int lY_l = lengthY_l;
      static constexpr int lZ_l = lengthZ_l;

      static constexpr int hX = 1;
      static constexpr int hY = 1;
      static constexpr int hZ = 1;

      static constexpr T inv_cs2 = 3.0;
      static constexpr T cs2 = 1./inv_cs2;
      static constexpr MathVector<MathVector<T, dimD>, dimQ> celerity =
        {
          MathVector<T, dimD>{{0.0, 0.0, 0.0}},
          MathVector<T, dimD>{{-1.0, 0.0, 0.0}},
          MathVector<T, dimD>{{0.0, -1.0, 0.0}},
          MathVector<T, dimD>{{0.0, 0.0, -1.0}},
          MathVector<T, dimD>{{-1.0, -1.0, 0.0}},
          MathVector<T, dimD>{{-1.0, 1.0, 0.0}},
          MathVector<T, dimD>{{-1.0, 0.0, -1.0}},
          MathVector<T, dimD>{{-1.0, 0.0, 1.0}},
          MathVector<T, dimD>{{0.0, -1.0, -1.0}},
          MathVector<T, dimD>{{0.0, -1.0, 1.0}},
          MathVector<T, dimD>{{-1.0, -1.0, -1.0}},
          MathVector<T, dimD>{{-1.0, -1.0, 1.0}},
          MathVector<T, dimD>{{-1.0, 1.0, -1.0}},
          MathVector<T, dimD>{{-1.0, 1.0, 1.0}},
          MathVector<T, dimD>{{1.0, 0.0, 0.0}},
          MathVector<T, dimD>{{0.0, 1.0, 0.0}},
          MathVector<T, dimD>{{0.0, 0.0, 1.0}},
          MathVector<T, dimD>{{1.0, 1.0, 0.0}},
          MathVector<T, dimD>{{1.0, -1.0, 0.0}},
          MathVector<T, dimD>{{1.0, 0.0, 1.0}},
          MathVector<T, dimD>{{1.0, 0.0, -1.0}},
          MathVector<T, dimD>{{0.0, 1.0, 1.0}},
          MathVector<T, dimD>{{0.0, 1.0, -1.0}},
          MathVector<T, dimD>{{1.0, 1.0, 1.0}},
          MathVector<T, dimD>{{1.0, 1.0, -1.0}},
          MathVector<T, dimD>{{1.0, -1.0, 1.0}},
          MathVector<T, dimD>{{1.0, -1.0, -1.0}}
        };
      static constexpr MathVector<T, dimQ> weight =
        {
          8./27., 2./27., 2./27.,
          2./27., 1./54., 1./54.,
          1./54., 1./54., 1./54.,
          1./54., 1./216., 1./216.,
          2./216., 1./216., 2./27.,
          2./27., 2./27., 1./54.,
          1./54., 1./54., 1./54.,
          1./54., 1./54., 1./216.,
          1./216., 1./216., 1./216.
        };
    };


  template <class T, LatticeType L>
    constexpr int dimD() { return Parameters<T, L>::dimD; }

  template <class T, LatticeType L>
    constexpr int dimQ() { return Parameters<T, L>::dimQ; }


  template <class T, LatticeType L>
    constexpr int lX_g() { return Parameters<T, L>::lX_g; }

  template <class T, LatticeType L>
    constexpr int lY_g() { return Parameters<T, L>::lY_g; }

  template <class T, LatticeType L>
    constexpr int lZ_g() { return Parameters<T, L>::lZ_g; }


  template <class T, LatticeType L>
    constexpr int lX_l() { return Parameters<T, L>::lX_l; }

  template <class T, LatticeType L>
    constexpr int lY_l() { return Parameters<T, L>::lY_l; }

  template <class T, LatticeType L>
    constexpr int lZ_l() { return Parameters<T, L>::lZ_l; }


  template <class T, LatticeType L>
    constexpr int hX() { return Parameters<T, L>::hX; }

  template <class T, LatticeType L>
    constexpr int hY() { return Parameters<T, L>::hY; }

  template <class T, LatticeType L>
    constexpr int hZ() { return Parameters<T, L>::hZ; }


  template <class T, LatticeType L>
    constexpr int sizeX_l() { return 2*hX<T, L>() + lX_l<T, L>(); }

  template <class T, LatticeType L>
    constexpr int sizeY_l() { return 2*hY<T, L>() + lY_l<T, L>(); }

  template <class T, LatticeType L>
    constexpr int sizeZ_l() { return 2*hZ<T, L>() + lZ_l<T, L>(); }

  template <class T, LatticeType L>
    constexpr int s_l() { return lX_l<T, L>()*lY_l<T, L>()*lZ_l<T, L>(); }

  template <class T, LatticeType L>
    constexpr int s_g() { return lX_g<T, L>()*lY_g<T, L>()*lZ_g<T, L>(); }


  template <class T, LatticeType L>
    constexpr T inv_cs2() { return Parameters<T, L>::inv_cs2; }

  template <class T, LatticeType L>
    constexpr T cs2() { return 1.0/Parameters<T, L>::inv_cs2; }

  template <class T, LatticeType L>
    constexpr T celerity(int d, int iQ) { return Parameters<T, L>::celerity[d][iQ]; }

  template <class T, LatticeType L>
    constexpr MathVector<T, dimD<T, L>()> celerity(int iQ) { return Parameters<T, L>::celerity[iQ]; }

  template <class T, LatticeType L>
    constexpr T weight(int iQ) { return Parameters<T, L>::weight[iQ]; }


  template<class T, LatticeType L>
    _HOST _DEVICE inline int idxL(const int iX, const int iY, const int iZ) {
    return sizeZ_l<T, L>()*(sizeY_l<T, L>()*iX + iY) + iZ;
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idx_inL(const int iX, const int iY, const int iZ) {
    return sizeZ_l<T, L>()*(sizeY_l<T, L>*(iX-hX<T, L>()) + (iY-hY<T, L>())) + iZ;
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idx_outL(const int iX, const int iY, const int iZ) {
    return sizeZ_l<T, L>()*(sizeY_l<T, L>()*(iX+hX<T, L>()) + (iY+hY<T, L>())) + (iZ+hZ<T, L>());
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idx_inF(const int iX, const int iY, const int iZ) {
    return lZ_l<T, L>()*(lY_l<T, L>()*(iX-hX<T, L>()) + (iY-hY<T, L>()) + (iZ-hZ<T, L>()));
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idx_lF(const int iX, const int iY, const int iZ) {
    return lZ_l<T, L>()*(lY_l<T, L>()*iX + iY) + iZ;
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idx_gF(const int iX, const int iY, const int iZ) {
    const int iX_l = iX - iX/lX_l<T, L>() * lX_l<T, L>();
    return iX/lX_l<T, L>() * s_l<T, L>() + idx_lF<T, L>(iX_l, iY, iZ);
  }

#if defined(_AOS)

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idxPop(const int idx, const int i) {
    return idx*dimQ<T, L>() + i;
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idxPop(const int iX, const int iY, const int iZ, const int i) {
    return idxL(iX, iY, iZ)*dimQ<T, L>() + i;
  }

  template<class T, LatticeType L>
    _HOST _DEVICE constexpr int idxPop_lF(const int idx, const int i) {
    return idx*dimQ<T, L>() + i;
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idxPop_gF(const int idx, const int i) {
    return idx*dimQ<T, L>() + i;
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idxPop_gF(const int iX, const int iY, const int iZ,
                                       const int i) {
    return idx_gF(iX, iY, iZ)*dimQ<T, L>() + i;
  }

#elif defined(_SOA)

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idxPop(const int idx, const int i) {
    return i*sizeX_l<T, L>()*sizeY_l<T, L>() + idx;
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idxPop(const int iX, const int iY, const int iZ,
                                    const int i) {
    return i*sizeX_l<T, L>()*sizeY_l<T, L>() + idxL<T, L>(iX, iY, iZ);
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idxPop_lF(const int idx, const int i) {
    return i*lX_l<T, L>()*lY_l<T, L>() + idx;
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idxPop_lF(const int iX, const int iY, const int iZ,
                                       const int i) {
    return i*lX_g<T, L>()*lY_g<T, L>()*lZ_g<T, L>() + idx_lF<T, L>(iX, iY, iZ);
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idxPop_gF(const int idx, const int i) {
    return i*lX_g<T, L>()*lY_g<T, L>()*lZ_g<T, L>() + idx;
  }

  template<class T, LatticeType L>
    _HOST _DEVICE inline int idxPop_gF(const int iX, const int iY, const int iZ,
                                       const int iQ) {
    const int iX_l = iX - iX/lX_l<T, L>() * lX_l<T, L>();
    return iX/lX_l<T, L>() * dimQ<T, L>()*s_l<T, L>() + idxPop_lF<T, L>(iX_l, iY, iZ, iQ);
  }


#else
# error "Error: data structure format is not defined."
#endif


}

#endif // COMMONS_H
