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

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 0, 0});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Left[1] = {1};
    static constexpr unsigned int iQ_Right[1] = {2};
    static constexpr unsigned int iQ_Bottom[0] = {};
    static constexpr unsigned int iQ_Top[0] = {};
    static constexpr unsigned int iQ_Front[0] = {};
    static constexpr unsigned int iQ_Back[0] = {};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
      return
        {
          MathVector<T, dimD>{{(T)0}},
          MathVector<T, dimD>{{(T)-1}},
          MathVector<T, dimD>{{(T)1}}
        };
    }

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        MathVector<T, dimQ>({
          (T)2/(T)3, (T)1/(T)6, (T)1/(T)6
            });
    }
  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D1Q3>::iQ_Left[1];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D1Q3>::iQ_Right[1];
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

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 1, 0});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Left[1] = {3};
    static constexpr unsigned int iQ_Right[1] = {1};
    static constexpr unsigned int iQ_Bottom[1] = {2};
    static constexpr unsigned int iQ_Top[1] = {4};
    static constexpr unsigned int iQ_Front[0] = {};
    static constexpr unsigned int iQ_Back[0] = {};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
      return
        {
          MathVector<T, dimD>{{(T)0, (T)0}},
            MathVector<T, dimD>{{(T)-1, (T)0}},
              MathVector<T, dimD>{{(T)0, (T)-1}},
                MathVector<T, dimD>{{(T)1, (T)0}},
                  MathVector<T, dimD>{{(T)0, (T)1}}
        };
    }

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        MathVector<T, dimQ>({
          (T)4/(T)6,
            (T)1/(T)12, (T)1/(T)12,
            (T)1/(T)12, (T)1/(T)12
            });
    }
  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q5>::iQ_Left[1];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q5>::iQ_Right[1];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q5>::iQ_Bottom[1];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q5>::iQ_Top[1];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q5>::iQ_Front[0];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q5>::iQ_Back[0];


  template <class T>
  struct Lattice<T, LatticeType::D2Q9> {
    static constexpr LatticeType Type = LatticeType::D2Q9;

    static constexpr int dimD = 2;
    static constexpr int dimQ = 9;

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 1, 0});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Left[3] = {1, 2, 3};
    static constexpr unsigned int iQ_Right[3] = {5, 6, 7};
    static constexpr unsigned int iQ_Bottom[3] = {3, 4, 5};
    static constexpr unsigned int iQ_Top[3] = {1, 7, 8};
    static constexpr unsigned int iQ_Front[0] = {};
    static constexpr unsigned int iQ_Back[0] = {};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
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

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        {
          (T)4/(T)9, (T)1/(T)36, (T)1/(T)9,
          (T)1/(T)36, (T)1/(T)9, (T)1/(T)36,
          (T)1/(T)9, (T)1/(T)36, (T)1/(T)9
        };
    }

  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q9>::iQ_Left[3];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q9>::iQ_Right[3];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q9>::iQ_Bottom[3];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q9>::iQ_Top[3];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q9>::iQ_Front[0];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D2Q9>::iQ_Back[0];


  template <class T>
  struct Lattice<T, LatticeType::D2Q37>
  {
    static constexpr LatticeType Type = LatticeType::D2Q37;

    static constexpr int dimD = 2;
    static constexpr int dimQ = 37;

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({3, 3, 0});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
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
          MathVector<T, dimD>{{(T)0, (T)1}},
          MathVector<T, dimD>{{(T)0, (T)0}},
          MathVector<T, dimD>{{(T)-1, (T)1}},
          MathVector<T, dimD>{{(T)-1, (T)0}},
          MathVector<T, dimD>{{(T)-1, (T)-1}},
          MathVector<T, dimD>{{(T)0, (T)-1}},
          MathVector<T, dimD>{{(T)1, (T)-1}},
          MathVector<T, dimD>{{(T)1, (T)0}},
          MathVector<T, dimD>{{(T)1, (T)1}},
          MathVector<T, dimD>{{(T)0, (T)1}},
          MathVector<T, dimD>{{(T)0, (T)0}},
          MathVector<T, dimD>{{(T)-1, (T)1}},
          MathVector<T, dimD>{{(T)-1, (T)0}},
          MathVector<T, dimD>{{(T)-1, (T)-1}},
          MathVector<T, dimD>{{(T)0, (T)-1}},
          MathVector<T, dimD>{{(T)1, (T)-1}},
          MathVector<T, dimD>{{(T)1, (T)0}},
          MathVector<T, dimD>{{(T)1, (T)1}},
          MathVector<T, dimD>{{(T)0, (T)1}},
          MathVector<T, dimD>{{(T)0, (T)0}},
          MathVector<T, dimD>{{(T)-1, (T)1}},
          MathVector<T, dimD>{{(T)-1, (T)0}},
          MathVector<T, dimD>{{(T)-1, (T)-1}},
          MathVector<T, dimD>{{(T)0, (T)-1}},
          MathVector<T, dimD>{{(T)1, (T)-1}},
          MathVector<T, dimD>{{(T)1, (T)0}},
          MathVector<T, dimD>{{(T)1, (T)1}},
          MathVector<T, dimD>{{(T)0, (T)1}},
          MathVector<T, dimD>{{(T)0, (T)1}}
        };
    }

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        MathVector<T, dimQ>({
          (T)4/(T)9, (T)1/(T)36, (T)1/(T)9,
            (T)1/(T)36, (T)1/(T)9, (T)1/(T)36,
            (T)1/(T)9, (T)1/(T)36, (T)1/(T)9,
          (T)4/(T)9, (T)1/(T)36, (T)1/(T)9,
            (T)1/(T)36, (T)1/(T)9, (T)1/(T)36,
            (T)1/(T)9, (T)1/(T)36, (T)1/(T)9,
          (T)4/(T)9, (T)1/(T)36, (T)1/(T)9,
            (T)1/(T)36, (T)1/(T)9, (T)1/(T)36,
            (T)1/(T)9, (T)1/(T)36, (T)1/(T)9,
          (T)4/(T)9, (T)1/(T)36, (T)1/(T)9,
            (T)1/(T)36, (T)1/(T)9, (T)1/(T)36,
            (T)1/(T)9, (T)1/(T)36, (T)1/(T)9,
            (T)1/(T)9
            });
    }

  };

  template <class T>
  struct Lattice<T, LatticeType::D3Q15>
  {
    static constexpr LatticeType Type = LatticeType::D3Q15;

    static constexpr int dimD = 3;
    static constexpr int dimQ = 15;

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 1, 1});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Left[5] = {1, 4, 5, 6, 7};
    static constexpr unsigned int iQ_Right[5] = {8, 11, 12, 13, 14};
    static constexpr unsigned int iQ_Bottom[5] = {2, 4, 5, 13, 14};
    static constexpr unsigned int iQ_Top[5] = {6, 7, 9, 11, 12};
    static constexpr unsigned int iQ_Front[5] = {3, 4, 6, 9, 10};
    static constexpr unsigned int iQ_Back[5] = {5, 7, 10, 11, 13};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
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

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        MathVector<T, dimQ>({
          (T)2/(T)9,
            (T)1/(T)9, (T)1/(T)9, (T)1/(T)9,
            (T)1/(T)72, (T)1/(T)72, (T)1/(T)72, (T)1/(T)72,
            (T)1/(T)9, (T)1/(T)9, (T)1/(T)9,
            (T)1/(T)72, (T)1/(T)72, (T)1/(T)72, (T)1/(T)72
            });
    }
  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q15>::iQ_Left[5];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q15>::iQ_Right[5];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q15>::iQ_Bottom[5];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q15>::iQ_Top[5];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q15>::iQ_Front[5];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q15>::iQ_Back[5];


  template <class T>
  struct Lattice<T, LatticeType::D3Q19>
  {
    static constexpr LatticeType Type = LatticeType::D3Q19;

    static constexpr int dimD = 3;
    static constexpr int dimQ = 19;

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 1, 1});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Left[5] = {1, 4, 5, 6, 7};
    static constexpr unsigned int iQ_Right[5] = {10, 13, 14, 15, 16};
    static constexpr unsigned int iQ_Bottom[5] = {2, 4, 8, 9, 14};
    static constexpr unsigned int iQ_Top[5] = {5, 11, 13, 17, 18};
    static constexpr unsigned int iQ_Front[5] = {3, 6, 8, 16, 18};
    static constexpr unsigned int iQ_Back[5] = {7, 9, 12, 15, 17};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
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

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        MathVector<T, dimQ>({
          (T)1/(T)3,
            (T)1/(T)18, (T)1/(T)18, (T)1/(T)18,
            (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
            (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
            (T)1/(T)18, (T)1/(T)18, (T)1/(T)18,
            (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
            (T)1/(T)36, (T)1/(T)36, (T)1/(T)36
            });
    }
  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q19>::iQ_Left[5];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q19>::iQ_Right[5];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q19>::iQ_Bottom[5];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q19>::iQ_Top[5];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q19>::iQ_Front[5];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q19>::iQ_Back[5];


  template <class T>
  struct Lattice<T, LatticeType::D3Q27>
  {
    static constexpr LatticeType Type = LatticeType::D3Q27;

    static constexpr int dimD = 3;
    static constexpr int dimQ = 27;

    DEVICE HOST
    static inline constexpr Position halo() {
      return Position({1, 1, 1});
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    static constexpr unsigned int iQ_Left[9] = {1, 4, 5, 6, 7, 10, 11, 12, 13};
    static constexpr unsigned int iQ_Right[9] = {14, 17, 18, 19, 20, 23, 24, 25, 26};
    static constexpr unsigned int iQ_Bottom[9] = {2, 4, 8, 9, 10, 11, 18, 25, 26};
    static constexpr unsigned int iQ_Top[9] = {5, 12, 13, 15, 17, 21, 22, 23, 24};
    static constexpr unsigned int iQ_Front[9] = {3, 6, 8, 10, 12, 20, 22, 24, 26};
    static constexpr unsigned int iQ_Back[9] = {7, 9, 11, 13, 16, 19, 21, 23, 25};

    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
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

    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        MathVector<T, dimQ>({
          (T)8/(T)27, (T)2/(T)27, (T)2/(T)27,
            (T)2/(T)27, (T)1/(T)54, (T)1/(T)54,
            (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
            (T)1/(T)54, (T)1/(T)216, (T)1/(T)216,
            (T)1/(T)216, (T)1/(T)216, (T)2/(T)27,
            (T)2/(T)27, (T)2/(T)27, (T)1/(T)54,
            (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
            (T)1/(T)54, (T)1/(T)54, (T)1/(T)216,
            (T)1/(T)216, (T)1/(T)216, (T)1/(T)216
            });
    }
  };

  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q27>::iQ_Left[9];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q27>::iQ_Right[9];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q27>::iQ_Bottom[9];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q27>::iQ_Top[9];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q27>::iQ_Front[9];
  template <class T> constexpr unsigned int Lattice<T, LatticeType::D3Q27>::iQ_Back[9];

  typedef Lattice<dataT, latticeT> L;
  typedef Lattice<unsigned int, latticeT> uiL;

}

#endif // LATTICE_H
