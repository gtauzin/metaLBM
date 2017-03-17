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
    struct Parameters
    {
      static constexpr int dimD = 1;
      static constexpr int dimQ = 1;
      static constexpr int lX_g = 1;
      static constexpr int lY_g = 1;
      static constexpr int lZ_g = 1;

      static constexpr int lX_l = 1;
      static constexpr int lY_l = 1;
      static constexpr int lZ_l = 1;

      static constexpr int hX = 0;
      static constexpr int hY = 0;
      static constexpr int hZ = 0;

      static constexpr T inv_cs2 = (T)1;
      static constexpr T cs2 = (T)1/inv_cs2;
      static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity()
        {
          return
            {
              MathVector<T, dimD>{{(T)0}}
            };
        };
      static inline constexpr MathVector<T, dimQ> weight()
      {
        return
          {
            (T)1
          };
      }
    };


  template <class T>
    struct Parameters<T, LatticeType::D1Q3>
    {
      static constexpr int dimD = 1;
      static constexpr int dimQ = 3;
      static constexpr int lX_g = lengthX_g;
      static constexpr int lY_g = 1;
      static constexpr int lZ_g = 1;

      static constexpr int lX_l = lX_g/NPROCS;
      static constexpr int lY_l = lY_g;
      static constexpr int lZ_l = lZ_g;

      static constexpr int hX = 1;
      static constexpr int hY = 0;
      static constexpr int hZ = 0;

      static constexpr T inv_cs2 = (T)3;
      static constexpr T cs2 = (T)1/inv_cs2;
      static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity()
        {
          return
            {
              MathVector<T, dimD>{{(T)0}},
              MathVector<T, dimD>{{(T)-1}},
              MathVector<T, dimD>{{(T)1}}
            };
        }
      static inline constexpr MathVector<T, dimQ> weight()
        {
          return
            {
              (T)2/(T)3, (T)1/(T)3, (T)1/(T)3
            };
        }
    };


    template <class T>
    struct Parameters<T, LatticeType::D2Q5>
    {
      static constexpr int dimD = 2;
      static constexpr int dimQ = 5;
      static constexpr int lX_g = lengthX_g;
      static constexpr int lY_g = lengthY_g;
      static constexpr int lZ_g = 1;

      static constexpr int lX_l = lX_g/NPROCS;
      static constexpr int lY_l = lY_g;
      static constexpr int lZ_l = lZ_g;

      static constexpr int hX = 1;
      static constexpr int hY = 1;
      static constexpr int hZ = 0;

      static constexpr T inv_cs2 = (T)3;
      static constexpr T cs2 = (T)1/inv_cs2;
      static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity()
        {
          return
            {
              MathVector<T, dimD>{{(T)0, (T)0}},
              MathVector<T, dimD>{{(T)-1, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)-1}},
              MathVector<T, dimD>{{(T)1, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)1}}
            };
        }
      static inline constexpr MathVector<T, dimQ> weight()
        {
          return
            {
              (T)4/(T)6,
              (T)1/(T)12, (T)1/(T)12,
              (T)1/(T)12, (T)1/(T)12
            };
        }

    };

  template <class T>
    struct Parameters<T, LatticeType::D2Q9>
    {
      static constexpr int dimD = 2;
      static constexpr int dimQ = 9;
      static constexpr int lX_g = lengthX_g;
      static constexpr int lY_g = lengthY_g;
      static constexpr int lZ_g = 1;

      static constexpr int lX_l = lX_g/NPROCS;
      static constexpr int lY_l = lY_g;
      static constexpr int lZ_l = lZ_g;

      static constexpr int hX = 1;
      static constexpr int hY = 1;
      static constexpr int hZ = 0;

      static constexpr T inv_cs2 = (T)3;
      static constexpr T cs2 = (T)1/inv_cs2;
      static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity()
        {
          return
            {
              MathVector<T, dimD>{{(T)0, (T)0}},
              MathVector<T, dimD>{{(T)-1, (T)1}},
              MathVector<T, dimD>{{(T)-1, (T)0}},
              MathVector<T, dimD>{{(T)-1, (T)-1}},
              MathVector<T, dimD>{{(T)0, (T)-1}},
              MathVector<T, dimD>{{(T)1, (T)-1}},
              MathVector<T, dimD>{{(T)1, (T)0}},
              MathVector<T, dimD>{{(T)1, (T)1}},
              MathVector<T, dimD>{{(T)0, (T)1}}
            };
        }
      static inline constexpr MathVector<T, dimQ> weight()
        {
          return
            {
              (T)4/(T)9, (T)1/(T)36, (T)1/(T)9,
              (T)1/(T)36, (T)1/(T)9, (T)1/(T)36,
              (T)1/(T)9, (T)1/(T)36, (T)1/(T)9
            };
        }

    };

  template <class T>
    struct Parameters<T, LatticeType::D3Q15>
    {
      static constexpr int dimD = 3;
      static constexpr int dimQ = 15;
      static constexpr int lX_g = lengthX_g;
      static constexpr int lY_g = lengthY_g;
      static constexpr int lZ_g = lengthZ_g;

      static constexpr int lX_l = lX_g/NPROCS;
      static constexpr int lY_l = lY_g;
      static constexpr int lZ_l = lZ_g;

      static constexpr int hX = 1;
      static constexpr int hY = 1;
      static constexpr int hZ = 1;

      static constexpr T inv_cs2 = (T)3;
      static constexpr T cs2 = (T)1/inv_cs2;
      static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity()
        {
          return
            {
              MathVector<T, dimD>{{(T)0, (T)0, (T)0}},
              MathVector<T, dimD>{{(T)-1, (T)0, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)-1, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)0, (T)-1}},
              MathVector<T, dimD>{{(T)-1, (T)-1, (T)-1}},
              MathVector<T, dimD>{{(T)-1, (T)-1, (T)1}},
              MathVector<T, dimD>{{(T)-1, (T)1, (T)-1}},
              MathVector<T, dimD>{{(T)-1, (T)1, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)0,(T)0}},
              MathVector<T, dimD>{{(T)0, (T)1, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)0, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)1, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)1, (T)-1}},
              MathVector<T, dimD>{{(T)1, (T)-1, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)-1, (T)-1}}
            };
        }
      static inline constexpr MathVector<T, dimQ> weight()
        {
          return
            {
              (T)2/(T)9,
              (T)1/(T)9, (T)1/(T)9, (T)1/(T)9,
              (T)1/(T)72, (T)1/(T)72, (T)1/(T)72, (T)1/(T)72,
              (T)1/(T)9, (T)1/(T)9, (T)1/(T)9,
              (T)1/(T)72, (T)1/(T)72, (T)1/(T)72, (T)1/(T)72
            };
        }
    };


  template <class T>
    struct Parameters<T, LatticeType::D3Q19>
    {
      static constexpr int dimD = 3;
      static constexpr int dimQ = 19;
      static constexpr int lX_g = lengthX_g;
      static constexpr int lY_g = lengthY_g;
      static constexpr int lZ_g = lengthZ_g;

      static constexpr int lX_l = lX_g/NPROCS;
      static constexpr int lY_l = lY_g;
      static constexpr int lZ_l = lZ_g;

      static constexpr int hX = 1;
      static constexpr int hY = 1;
      static constexpr int hZ = 1;

      static constexpr T inv_cs2 = (T)3;
      static constexpr T cs2 = (T)1/inv_cs2;
      static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity()
        {
          return
            {
              MathVector<T, dimD>{{(T)0, (T)0, (T)0}},
              MathVector<T, dimD>{{(T)-1, (T)0, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)-1, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)0, (T)-1}},
              MathVector<T, dimD>{{(T)-1, (T)-1, (T)0}},
              MathVector<T, dimD>{{(T)-1, (T)1, (T)0}},
              MathVector<T, dimD>{{(T)-1, (T)0, (T)-1}},
              MathVector<T, dimD>{{(T)-1, (T)0, (T)1}},
              MathVector<T, dimD>{{(T)0, (T)-1,(T)-1}},
              MathVector<T, dimD>{{(T)0, (T)-1, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)0, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)1, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)0, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)1, (T)0}},
              MathVector<T, dimD>{{(T)1, (T)-1, (T)0}},
              MathVector<T, dimD>{{(T)1, (T)0, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)0, (T)-1}},
              MathVector<T, dimD>{{(T)0, (T)1, (T)1}},
              MathVector<T, dimD>{{(T)0, (T)1, (T)-1}}

            };
        }
      static inline constexpr MathVector<T, dimQ> weight()
        {
          return
            {
              (T)1/(T)3,
              (T)1/(T)18, (T)1/(T)18, (T)1/(T)18,
              (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
              (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
              (T)1/(T)18, (T)1/(T)18, (T)1/(T)18,
              (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
              (T)1/(T)36, (T)1/(T)36, (T)1/(T)36
            };
        }
    };


   template <class T>
    struct Parameters<T, LatticeType::D3Q27>
    {
      static constexpr int dimD = 3;
      static constexpr int dimQ = 27;
      static constexpr int lX_g = lengthX_g;
      static constexpr int lY_g = lengthY_g;
      static constexpr int lZ_g = lengthZ_g;

      static constexpr int lX_l = lX_g/NPROCS;
      static constexpr int lY_l = lY_g;
      static constexpr int lZ_l = lZ_g;

      static constexpr int hX = 1;
      static constexpr int hY = 1;
      static constexpr int hZ = 1;

      static constexpr T inv_cs2 = (T)3;
      static constexpr T cs2 = (T)1/inv_cs2;
      static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity()
        {
          return
            {
              MathVector<T, dimD>{{(T)0, (T)0, (T)0}},
              MathVector<T, dimD>{{(T)-1, (T)0, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)-1, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)0, (T)-1}},
              MathVector<T, dimD>{{(T)-1, (T)-1, (T)0}},
              MathVector<T, dimD>{{(T)-1, (T)1, (T)0}},
              MathVector<T, dimD>{{(T)-1, (T)0, (T)-1}},
              MathVector<T, dimD>{{(T)-1, (T)0, (T)1}},
              MathVector<T, dimD>{{(T)0, (T)-1, (T)-1}},
              MathVector<T, dimD>{{(T)0, (T)-1, (T)1}},
              MathVector<T, dimD>{{(T)-1, (T)-1, (T)-1}},
              MathVector<T, dimD>{{(T)-1, (T)-1, (T)1}},
              MathVector<T, dimD>{{(T)-1, (T)1, (T)-1}},
              MathVector<T, dimD>{{(T)-1, (T)1, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)0, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)1, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)0, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)1, (T)0}},
              MathVector<T, dimD>{{(T)1, (T)-1, (T)0}},
              MathVector<T, dimD>{{(T)1, (T)0, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)0, (T)-1}},
              MathVector<T, dimD>{{(T)0, (T)1, (T)1}},
              MathVector<T, dimD>{{(T)0, (T)1, (T)-1}},
              MathVector<T, dimD>{{(T)1, (T)1, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)1, (T)-1}},
              MathVector<T, dimD>{{(T)1, (T)-1, (T)1}},
              MathVector<T, dimD>{{(T)1, (T)-1, (T)-1}}
            };
        }
      static inline constexpr MathVector<T, dimQ> weight()
        {
          return
            {
              (T)8/(T)27, (T)2/(T)27, (T)2/(T)27,
              (T)2/(T)27, (T)1/(T)54, (T)1/(T)54,
              (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
              (T)1/(T)54, (T)1/(T)216, (T)1/(T)216,
              (T)1/(T)216, (T)1/(T)216, (T)2/(T)27,
              (T)2/(T)27, (T)2/(T)27, (T)1/(T)54,
              (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
              (T)1/(T)54, (T)1/(T)54, (T)1/(T)216,
              (T)1/(T)216, (T)1/(T)216, (T)1/(T)216
            };
        }
    };

  typedef Parameters<valueType, latticeType> P;


  constexpr int sizeX_l() { return 2*P::hX + P::lX_l; }

  constexpr int sizeY_l() { return 2*P::hY + P::lY_l; }

  constexpr int sizeZ_l() { return 2*P::hZ + P::lZ_l; }

  constexpr int s_l() { return P::lX_l * P::lY_l * P::lZ_l; }

  constexpr int s_g() { return P::lX_g * P::lY_g * P::lZ_g; }



  _HOST _DEVICE inline int idxL(const int iX, const int iY, const int iZ) {
    return sizeZ_l()*(sizeY_l()*iX + iY) + iZ;
  }

  _HOST _DEVICE inline int idxL(const MathVector<int, 3>& iP) {
    return sizeZ_l()*(sizeY_l()*iP[d::X] + iP[d::Y]) + iP[d::Z];
  }

  _HOST _DEVICE inline int idx_inL(const int iX, const int iY, const int iZ) {
    return sizeZ_l()*(sizeY_l()*(iX-P::hX) + (iY-P::hY)) + (iZ-P::hZ);
  }

  _HOST _DEVICE inline int idx_outL(const int iX, const int iY, const int iZ) {
    return sizeZ_l()*(sizeY_l()*(iX+P::hX) + (iY+P::hY)) + (iZ+P::hZ);
  }

  _HOST _DEVICE inline int idx_inF(const int iX, const int iY, const int iZ) {
    return P::lZ_l * (P::lY_l * (iX-P::hX) + (iY-P::hY)) + (iZ-P::hZ);
  }

  _HOST _DEVICE inline int idx_lF(const int iX, const int iY, const int iZ) {
    return P::lZ_l * (P::lY_l*iX + iY) + iZ;
  }

  _HOST _DEVICE inline int idx_gF(const int iX, const int iY, const int iZ) {
    const int iX_l = iX - iX/P::lX_l * P::lX_l;
    return iX/P::lX_l * s_l() + idx_lF(iX_l, iY, iZ);
  }

#if defined(_AOS)

  _HOST _DEVICE inline int idxPop(const int idx, const int iQ) {
    return idx * P::dimQ + iQ;
  }

  _HOST _DEVICE inline int idxPop(const int iX, const int iY, const int iZ,
                                  const int iQ) {
    return idxL(iX, iY, iZ) * P::dimQ + iQ;
  }

  _HOST _DEVICE constexpr int idxPop_lF(const int idx, const int iQ) {
    return idx * P::dimQ + iQ;
  }

  _HOST _DEVICE inline int idxPop_gF(const int idx, const int iQ) {
    return idx * P::dimQ + iQ;
  }

  _HOST _DEVICE inline int idxPop_gF(const int iX, const int iY, const int iZ,
                                     const int iQ) {
    return idx_gF(iX, iY, iZ) * P::dimQ + iQ;
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
    return iQ * P::lX_l * P::lY_l * P::lZ_l + idx;
  }

  _HOST _DEVICE inline int idxPop_lF(const int iX, const int iY, const int iZ,
                                     const int iQ) {
    return iQ * P::lX_g * P::lY_g * P::lZ_g + idx_lF(iX, iY, iZ);
  }

  _HOST _DEVICE inline int idxPop_gF(const int idx, const int iQ) {
    return iQ * P::lX_g * P::lY_g * P::lZ_g + idx;
  }

  _HOST _DEVICE inline int idxPop_gF(const int iX, const int iY, const int iZ,
                                     const int iQ) {
    const int iX_l = iX - iX/P::lX_l * P::lX_l;
    return iX/P::lX_l * P::dimQ * s_l() + idxPop_lF(iX_l, iY, iZ, iQ);
  }


#else
# error "Error: data structure format is not defined."
#endif


}

#endif // COMMONS_H
