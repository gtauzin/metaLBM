#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <string>
#include <stdlib.h>

#include "MathVector.h"
#include "Commons.h"
#include "Options.h"

namespace lbm {

  /**
   * Lattice required to define a lattice for a parallel code.
   *
   * @tparam T data type.
   * @tparam LatticeT lattice type of the form DdQq.
   */

  template <class T, LatticeType LatticeT>
  struct Lattice{};

  template <class T>
  struct Lattice<T, LatticeType::D1Q3>
  {
    static constexpr LatticeType Type = LatticeType::D1Q3;
    static constexpr int dimD = 1;
    static constexpr int dimQ = 3;
    static constexpr int faceQ = 1;

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 0, 0});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Bottom[0] = {};
    static constexpr unsigned int iQ_Top[0] = {};
    static constexpr unsigned int iQ_Front[0] = {};
    static constexpr unsigned int iQ_Back[0] = {};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
      return
        {
          MathVector<T, dimD>{{(T) 0}},
          MathVector<T, dimD>{{(T)-1}},
          MathVector<T, dimD>{{(T) 1}}
        };
    }

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        {
          (T)2/(T)3, (T)1/(T)6, (T)1/(T)6
        };
    }
  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D1Q3>::iQ_Bottom[0];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D1Q3>::iQ_Top[0];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D1Q3>::iQ_Front[0];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D1Q3>::iQ_Back[0];


  template <class T>
  struct Lattice<T, LatticeType::D2Q5>
  {
    static constexpr LatticeType Type = LatticeType::D2Q5;
    static constexpr int dimD = 2;
    static constexpr int dimQ = 5;
    static constexpr int faceQ = 1;

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 1, 0});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Bottom[faceQ] = {3};
    static constexpr unsigned int iQ_Top[faceQ] = {4};
    static constexpr unsigned int iQ_Front[0] = {};
    static constexpr unsigned int iQ_Back[0] = {};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
      return
        {
          MathVector<T, dimD>{{(T) 0, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T) 0}},
          MathVector<T, dimD>{{(T) 1, (T) 0}},
          MathVector<T, dimD>{{(T) 0, (T)-1}},
          MathVector<T, dimD>{{(T) 0, (T) 1}}
        };
    }

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        {
          (T)4/(T)6,
          (T)1/(T)12, (T)1/(T)12,
          (T)1/(T)12, (T)1/(T)12
        };
    }
  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q5>::iQ_Bottom[Lattice<T, LatticeType::D2Q5>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q5>::iQ_Top[Lattice<T, LatticeType::D2Q5>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q5>::iQ_Front[0];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q5>::iQ_Back[0];


  template <class T>
  struct Lattice<T, LatticeType::D2Q9> {
    static constexpr LatticeType Type = LatticeType::D2Q9;

    static constexpr int dimD = 2;
    static constexpr int dimQ = 9;
    static constexpr int faceQ = 3;

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 1, 0});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Bottom[faceQ] = {3, 4, 7};
    static constexpr unsigned int iQ_Top[faceQ] = {1, 6, 8};
    static constexpr unsigned int iQ_Front[0] = {};
    static constexpr unsigned int iQ_Back[0] = {};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
      return
        {
          MathVector<T, dimD>{{(T) 0, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T) 1}},
          MathVector<T, dimD>{{(T)-1, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T)-1}},
          MathVector<T, dimD>{{(T) 1, (T)-1}},
          MathVector<T, dimD>{{(T) 1, (T) 0}},
          MathVector<T, dimD>{{(T) 1, (T) 1}},
          MathVector<T, dimD>{{(T) 0, (T)-1}},
          MathVector<T, dimD>{{(T) 0, (T) 1}}
        };
    }

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        {
          (T)4/(T)9, (T)1/(T)36, (T)1/(T)9,
          (T)1/(T)36, (T)1/(T)36, (T)1/(T)9,
          (T)1/(T)36, (T)1/(T)9, (T)1/(T)9
        };
    }

  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q9>::iQ_Bottom[Lattice<T, LatticeType::D2Q9>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q9>::iQ_Top[Lattice<T, LatticeType::D2Q9>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q9>::iQ_Front[0];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q9>::iQ_Back[0];


  template <class T>
  struct Lattice<T, LatticeType::D3Q15>
  {
    static constexpr LatticeType Type = LatticeType::D3Q15;

    static constexpr int dimD = 3;
    static constexpr int dimQ = 15;
    static constexpr int faceQ = 5;

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 1, 1});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Bottom[faceQ] = {2, 3, 9, 10, 11};
    static constexpr unsigned int iQ_Top[faceQ] = {4, 5, 7, 8, 13};
    static constexpr unsigned int iQ_Front[faceQ] = {2, 4, 8, 10, 12};
    static constexpr unsigned int iQ_Back[faceQ] = {3, 5, 7, 9, 14};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
      return
        {
          MathVector<T, dimD>{{(T) 0, (T) 0, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T) 0, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T)-1, (T)-1}},
          MathVector<T, dimD>{{(T)-1, (T)-1, (T) 1}},
          MathVector<T, dimD>{{(T)-1, (T) 1, (T)-1}},
          MathVector<T, dimD>{{(T)-1, (T) 1, (T) 1}},
          MathVector<T, dimD>{{(T) 1, (T) 0, (T) 0}},
          MathVector<T, dimD>{{(T) 1, (T) 1, (T) 1}},
          MathVector<T, dimD>{{(T) 1, (T) 1, (T)-1}},
          MathVector<T, dimD>{{(T) 1, (T)-1, (T) 1}},
          MathVector<T, dimD>{{(T) 1, (T)-1, (T)-1}},
          MathVector<T, dimD>{{(T) 0, (T)-1, (T) 0}},
          MathVector<T, dimD>{{(T) 0, (T) 0, (T)-1}},
          MathVector<T, dimD>{{(T) 0, (T) 1, (T) 0}},
          MathVector<T, dimD>{{(T) 0, (T) 0, (T) 1}},
        };
    }

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        {
          (T)2/(T)9,
          (T)1/(T)9, (T)1/(T)72, (T)1/(T)72,
          (T)1/(T)72, (T)1/(T)72, (T)1/(T)9, (T)1/(T)72,
          (T)1/(T)72, (T)1/(T)72, (T)1/(T)72,
          (T)1/(T)9, (T)1/(T)9, (T)1/(T)9, (T)1/(T)9
        };
    }
  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q15>::iQ_Bottom[Lattice<T, LatticeType::D3Q15>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q15>::iQ_Top[Lattice<T, LatticeType::D3Q15>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q15>::iQ_Front[Lattice<T, LatticeType::D3Q15>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q15>::iQ_Back[Lattice<T, LatticeType::D3Q15>::faceQ];


  template <class T>
  struct Lattice<T, LatticeType::D3Q19>
  {
    static constexpr LatticeType Type = LatticeType::D3Q19;

    static constexpr int dimD = 3;
    static constexpr int dimQ = 19;
    static constexpr int faceQ = 5;

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 1, 1});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Bottom[faceQ] = {2, 8, 11, 13, 14};
    static constexpr unsigned int iQ_Top[faceQ] = {3, 7, 15, 17, 18};
    static constexpr unsigned int iQ_Front[faceQ] = {4, 10, 12, 13, 18};
    static constexpr unsigned int iQ_Back[faceQ] = {5, 9, 14, 16, 17};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
      return
        {
          MathVector<T, dimD>{{(T) 0, (T) 0, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T) 0, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T)-1, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T) 1, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T) 0, (T)-1}},
          MathVector<T, dimD>{{(T)-1, (T) 0, (T) 1}},
          MathVector<T, dimD>{{(T) 1, (T) 0, (T) 0}},
          MathVector<T, dimD>{{(T) 1, (T) 1, (T) 0}},
          MathVector<T, dimD>{{(T) 1, (T)-1, (T) 0}},
          MathVector<T, dimD>{{(T) 1, (T) 0, (T) 1}},
          MathVector<T, dimD>{{(T) 1, (T) 0, (T)-1}},
          MathVector<T, dimD>{{(T) 0, (T)-1, (T) 0}},
          MathVector<T, dimD>{{(T) 0, (T) 0, (T)-1}},
          MathVector<T, dimD>{{(T) 0, (T)-1, (T)-1}},
          MathVector<T, dimD>{{(T) 0, (T)-1, (T) 1}},
          MathVector<T, dimD>{{(T) 0, (T) 1, (T) 0}},
          MathVector<T, dimD>{{(T) 0, (T) 0, (T) 1}},
          MathVector<T, dimD>{{(T) 0, (T) 1, (T) 1}},
          MathVector<T, dimD>{{(T) 0, (T) 1, (T) -1}}
        };
    }

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        {
          (T)1/(T)3,
          (T)1/(T)18, (T)1/(T)36, (T)1/(T)36,
          (T)1/(T)36, (T)1/(T)36, (T)1/(T)18,
          (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
          (T)1/(T)36, (T)1/(T)18, (T)1/(T)18,
          (T)1/(T)36, (T)1/(T)36, (T)1/(T)18,
          (T)1/(T)18, (T)1/(T)36, (T)1/(T)36
        };
    }
  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q19>::iQ_Bottom[Lattice<T, LatticeType::D3Q19>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q19>::iQ_Top[Lattice<T, LatticeType::D3Q19>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q19>::iQ_Front[Lattice<T, LatticeType::D3Q19>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q19>::iQ_Back[Lattice<T, LatticeType::D3Q19>::faceQ];


  template <class T>
  struct Lattice<T, LatticeType::D3Q27>
  {
    static constexpr LatticeType Type = LatticeType::D3Q27;

    static constexpr int dimD = 3;
    static constexpr int dimQ = 27;
    static constexpr int faceQ = 9;

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 1, 1});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Bottom[faceQ] = {2, 6, 7, 12, 17, 18, 19, 21, 22};
    static constexpr unsigned int iQ_Top[faceQ] = {3, 8, 9, 11, 15, 16, 23, 25, 26};
    static constexpr unsigned int iQ_Front[faceQ] = {4, 6, 8, 14, 16, 18, 20, 21, 26};
    static constexpr unsigned int iQ_Back[faceQ] = {5, 7, 9, 13, 15, 17, 22, 24, 25};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
      return
        {
          MathVector<T, dimD>{{(T) 0, (T) 0, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T) 0, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T)-1, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T) 1, (T) 0}},
          MathVector<T, dimD>{{(T)-1, (T) 0, (T)-1}},
          MathVector<T, dimD>{{(T)-1, (T) 0, (T) 1}},
          MathVector<T, dimD>{{(T)-1, (T)-1, (T)-1}},
          MathVector<T, dimD>{{(T)-1, (T)-1, (T) 1}},
          MathVector<T, dimD>{{(T)-1, (T) 1, (T)-1}},
          MathVector<T, dimD>{{(T)-1, (T) 1, (T) 1}},
          MathVector<T, dimD>{{(T) 1, (T) 0, (T) 0}},
          MathVector<T, dimD>{{(T) 1, (T) 1, (T) 0}},
          MathVector<T, dimD>{{(T) 1, (T)-1, (T) 0}},
          MathVector<T, dimD>{{(T) 1, (T) 0, (T) 1}},
          MathVector<T, dimD>{{(T) 1, (T) 0, (T)-1}},
          MathVector<T, dimD>{{(T) 1, (T) 1, (T) 1}},
          MathVector<T, dimD>{{(T) 1, (T) 1, (T)-1}},
          MathVector<T, dimD>{{(T) 1, (T)-1, (T) 1}},
          MathVector<T, dimD>{{(T) 1, (T)-1, (T)-1}},
          MathVector<T, dimD>{{(T) 0, (T)-1, (T) 0}},
          MathVector<T, dimD>{{(T) 0, (T) 0, (T)-1}},
          MathVector<T, dimD>{{(T) 0, (T)-1, (T)-1}},
          MathVector<T, dimD>{{(T) 0, (T)-1, (T) 1}},
          MathVector<T, dimD>{{(T) 0, (T) 1, (T) 0}},
          MathVector<T, dimD>{{(T) 0, (T) 0, (T) 1}},
          MathVector<T, dimD>{{(T) 0, (T) 1, (T) 1}},
          MathVector<T, dimD>{{(T) 0, (T) 1, (T)-1}},
        };
    }

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        {
          (T)8/(T)27, (T)2/(T)27, (T)1/(T)54,
          (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
          (T)1/(T)216, (T)1/(T)216, (T)1/(T)216,
          (T)1/(T)216, (T)2/(T)27, (T)1/(T)54,
          (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
          (T)1/(T)216, (T)1/(T)216, (T)1/(T)216,
          (T)1/(T)216, (T)2/(T)27, (T)2/(T)27,
          (T)1/(T)54, (T)1/(T)54, (T)2/(T)27,
          (T)2/(T)27, (T)1/(T)54, (T)1/(T)54
        };
    }
  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q27>::iQ_Bottom[Lattice<T, LatticeType::D3Q27>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q27>::iQ_Top[Lattice<T, LatticeType::D3Q27>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q27>::iQ_Front[Lattice<T, LatticeType::D3Q27>::faceQ];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q27>::iQ_Back[Lattice<T, LatticeType::D3Q27>::faceQ];

  typedef Lattice<dataT, latticeT> L;
  typedef Lattice<unsigned int, latticeT> uiL;

}

#endif // LATTICE_H
