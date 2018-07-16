#pragma once

#include <stdlib.h>
#include <iostream>
#include <string>

#include "Commons.h"
#include "MathVector.h"
#include "Options.h"

namespace lbm {

/**
 * Parameters required to define a lattice stencil.
 *
 * @tparam T data type.
 * @tparam LatticeT lattice type of the form DdQq.
 */
template <class T, LatticeType LatticeT>
struct Lattice {};

template <class T>
struct Lattice<T, LatticeType::D1Q3> {
  static constexpr LatticeType Type = LatticeType::D1Q3;
  static constexpr int dimD = 1;
  static constexpr int dimQ = 3;
  static constexpr int faceQ = 1;

  LBM_DEVICE LBM_HOST
  static inline constexpr Position halo() {
    return Position({1, 0, 0});
  }

  static constexpr T inv_cs2 = (T)3;
  static constexpr T cs2 = (T)1 / inv_cs2;

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Bottom() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Top() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Front() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Back() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
    return {MathVector<T, dimD>{{(T)0}}, MathVector<T, dimD>{{(T)-1}},
            MathVector<T, dimD>{{(T)1}}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, dimQ> weight() {
    return {(T)2 / (T)3, (T)1 / (T)6, (T)1 / (T)6};
  }
};

template <class T>
struct Lattice<T, LatticeType::D2Q5> {
  static constexpr LatticeType Type = LatticeType::D2Q5;
  static constexpr int dimD = 2;
  static constexpr int dimQ = 5;
  static constexpr int faceQ = 1;

  LBM_DEVICE LBM_HOST
  static inline constexpr Position halo() {
    return Position({1, 1, 0});
  }

  static constexpr T inv_cs2 = (T)3;
  static constexpr T cs2 = (T)1 / inv_cs2;

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ>
  iQ_Bottom() {
    return {{3}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ>
  iQ_Top() {
    return {{4}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Front() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Back() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
    return {
        MathVector<T, dimD>{{(T)0, (T)0}},
        MathVector<T, dimD>{{(T)-1, (T)0}},
        MathVector<T, dimD>{{(T)1, (T)0}},
        MathVector<T, dimD>{{(T)0, (T)-1}},
        MathVector<T, dimD>{{(T)0, (T)1}}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, dimQ> weight() {
    return {(T)4 / (T)6, (T)1 / (T)12, (T)1 / (T)12, (T)1 / (T)12,
            (T)1 / (T)12};
  }
};

template <class T>
struct Lattice<T, LatticeType::D2Q9> {
  static constexpr LatticeType Type = LatticeType::D2Q9;

  static constexpr int dimD = 2;
  static constexpr int dimQ = 9;
  static constexpr int faceQ = 3;

  LBM_DEVICE LBM_HOST
  static inline constexpr Position halo() {
    return Position({1, 1, 0});
  }

  static constexpr T inv_cs2 = (T)3;
  static constexpr T cs2 = (T)1 / inv_cs2;

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, faceQ> iQ_Bottom() {
    return {{3, 4, 7}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, faceQ> iQ_Top() {
    return {{1, 6, 8}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Front() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Back() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
    return {
        MathVector<T, dimD>{{(T) 0, (T) 0}},
        MathVector<T, dimD>{{(T)-1, (T) 1}},
        MathVector<T, dimD>{{(T)-1, (T) 0}},
        MathVector<T, dimD>{{(T)-1, (T)-1}},
        MathVector<T, dimD>{{(T) 1, (T)-1}},
        MathVector<T, dimD>{{(T) 1, (T) 0}},
        MathVector<T, dimD>{{(T) 1, (T) 1}},
        MathVector<T, dimD>{{(T) 0, (T)-1}},
        MathVector<T, dimD>{{(T) 0, (T) 1}}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, dimQ> weight() {
    return {(T)4 / (T)9, (T)1 / (T)36, (T)1 / (T)9, (T)1 / (T)36, (T)1 / (T)36,
            (T)1 / (T)9, (T)1 / (T)36, (T)1 / (T)9, (T)1 / (T)9};
  }
};


 template <class T>
struct Lattice<T, LatticeType::D2Q13> {
  static constexpr LatticeType Type = LatticeType::D2Q13;

  static constexpr int dimD = 2;
  static constexpr int dimQ = 13;
  static constexpr int faceQ = 4;

  LBM_DEVICE LBM_HOST static inline constexpr Position halo() {
    return Position({2, 2, 0});
  }

  static constexpr T cs2 = (T) 0.4;
  static constexpr T inv_cs2 = (T)1 / (T)cs2;

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Bottom() {
    return {{2, 6, 9, 11}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Top() {
    return {{3, 7, 10, 12}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Front() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Back() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
    return {
        MathVector<T, dimD>{{(T) 0, (T) 0}},
        MathVector<T, dimD>{{(T)-1, (T) 0}},
        MathVector<T, dimD>{{(T)-1, (T)-1}},
        MathVector<T, dimD>{{(T)-1, (T) 1}},
        MathVector<T, dimD>{{(T)-2, (T) 0}},
        MathVector<T, dimD>{{(T) 1, (T) 0}},
        MathVector<T, dimD>{{(T) 1, (T)-1}},
        MathVector<T, dimD>{{(T) 1, (T) 1}},
        MathVector<T, dimD>{{(T) 2, (T) 0}},
        MathVector<T, dimD>{{(T) 0, (T)-1}},
        MathVector<T, dimD>{{(T) 0, (T) 1}},
        MathVector<T, dimD>{{(T) 0, (T)-2}},
        MathVector<T, dimD>{{(T) 0, (T) 2}},
    };
  }

  static constexpr T w0 = (T) 5 / (T) 24;
  static constexpr T w1 = (T) 4 / (T) 15;
  static constexpr T w2 = (T) 1 / (T) 10;
  static constexpr T w4 = (T) 1 / (T) 120;

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, dimQ> weight() {
    return {w0, w1, w2, w2, w4, w1, w2, w2, w4, w1, w1, w4, w4};
  }
};


template <class T>
struct Lattice<T, LatticeType::D2Q17> {
  static constexpr LatticeType Type = LatticeType::D2Q17;

  static constexpr int dimD = 2;
  static constexpr int dimQ = 17;
  static constexpr int faceQ = 7;

  LBM_DEVICE LBM_HOST
  static inline constexpr Position halo() {
    return Position({3, 3, 0});
  }

  static constexpr T inv_cs2 = (T)2 / (T)3;
  static constexpr T cs2 = (T)1 / inv_cs2;

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Bottom() {
    return {{1, 3, 6, 8, 10, 13, 15}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Top() {
    return {{2, 4, 7, 9, 11, 14, 16}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Front() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, 0> iQ_Back() {
    return {{}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
    return {
        MathVector<T, dimD>{{(T) 0, (T) 0}},
        MathVector<T, dimD>{{(T)-1, (T)-1}},
        MathVector<T, dimD>{{(T)-1, (T) 1}},
        MathVector<T, dimD>{{(T)-2, (T)-2}},
        MathVector<T, dimD>{{(T)-2, (T) 2}},
        MathVector<T, dimD>{{(T)-3, (T) 0}},
        MathVector<T, dimD>{{(T)-3, (T)-3}},
        MathVector<T, dimD>{{(T)-3, (T) 3}},
        MathVector<T, dimD>{{(T) 1, (T)-1}},
        MathVector<T, dimD>{{(T) 1, (T) 1}},
        MathVector<T, dimD>{{(T) 2, (T)-2}},
        MathVector<T, dimD>{{(T) 2, (T) 2}},
        MathVector<T, dimD>{{(T) 3, (T) 0}},
        MathVector<T, dimD>{{(T) 3, (T)-3}},
        MathVector<T, dimD>{{(T) 3, (T) 3}},
        MathVector<T, dimD>{{(T) 0, (T)-3}},
        MathVector<T, dimD>{{(T) 0, (T) 3}},
    };
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, dimQ> weight() {
    return {(T) 0.121527777777777777777778, (T) 0.175781250000000000000000,
            (T) 0.175781250000000000000000, (T) 0.014062500000000000000000,
            (T) 0.014062500000000000000000, (T) 0.027777777777777777777778,
            (T) 0.001996527777777777777778, (T) 0.001996527777777777777778,
            (T) 0.175781250000000000000000, (T) 0.175781250000000000000000,
            (T) 0.014062500000000000000000, (T) 0.014062500000000000000000,
            (T) 0.027777777777777777777778, (T) 0.001996527777777777777778,
            (T) 0.001996527777777777777778, (T) 0.027777777777777777777778,
            (T) 0.027777777777777777777778};
  }
};


template <class T>
struct Lattice<T, LatticeType::D3Q15> {
  static constexpr LatticeType Type = LatticeType::D3Q15;

  static constexpr int dimD = 3;
  static constexpr int dimQ = 15;
  static constexpr int faceQ = 5;

  LBM_DEVICE LBM_HOST
  static inline constexpr Position halo() {
    return Position({1, 1, 1});
  }

  static constexpr T inv_cs2 = (T)3;
  static constexpr T cs2 = (T)1 / inv_cs2;

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, faceQ> iQ_Bottom() {
    return {{2, 3, 9, 10, 11}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, faceQ> iQ_Top() {
    return {{4, 5, 7, 8, 13}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, faceQ> iQ_Front() {
    return {{2, 4, 8, 10, 12}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, faceQ> iQ_Back() {
    return {{3, 5, 7, 9, 14}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
    return {
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

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, dimQ> weight() {
    return {(T)2 / (T)9,  (T)1 / (T)9, (T)1 / (T)72, (T)1 / (T)72, (T)1 / (T)72,
            (T)1 / (T)72, (T)1 / (T)9, (T)1 / (T)72, (T)1 / (T)72, (T)1 / (T)72,
            (T)1 / (T)72, (T)1 / (T)9, (T)1 / (T)9,  (T)1 / (T)9,  (T)1 / (T)9};
  }
};

template <class T>
struct Lattice<T, LatticeType::D3Q19> {
  static constexpr LatticeType Type = LatticeType::D3Q19;

  static constexpr int dimD = 3;
  static constexpr int dimQ = 19;
  static constexpr int faceQ = 5;

  LBM_DEVICE LBM_HOST
  static inline constexpr Position halo() {
    return Position({1, 1, 1});
  }

  static constexpr T inv_cs2 = (T)3;
  static constexpr T cs2 = (T)1 / inv_cs2;

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Bottom() {
    return {{2, 8, 11, 13, 14}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Top() {
    return {{3, 7, 15, 17, 18}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Front() {
    return {{4, 10, 12, 13, 18}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Back() {
    return {{5, 9, 14, 16, 17}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
    return {MathVector<T, dimD>{{(T) 0, (T) 0, (T) 0}},
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
            MathVector<T, dimD>{{(T) 0, (T) 1, (T)-1}}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, dimQ> weight() {
    return {(T)1 / (T)3,  (T)1 / (T)18, (T)1 / (T)36, (T)1 / (T)36,
            (T)1 / (T)36, (T)1 / (T)36, (T)1 / (T)18, (T)1 / (T)36,
            (T)1 / (T)36, (T)1 / (T)36, (T)1 / (T)36, (T)1 / (T)18,
            (T)1 / (T)18, (T)1 / (T)36, (T)1 / (T)36, (T)1 / (T)18,
            (T)1 / (T)18, (T)1 / (T)36, (T)1 / (T)36};
  }
};

template <class T>
struct Lattice<T, LatticeType::D3Q27> {
  static constexpr LatticeType Type = LatticeType::D3Q27;

  static constexpr int dimD = 3;
  static constexpr int dimQ = 27;
  static constexpr int faceQ = 9;

  LBM_DEVICE LBM_HOST
  static inline constexpr Position halo() {
    return Position({1, 1, 1});
  }

  static constexpr T inv_cs2 = (T)3;
  static constexpr T cs2 = (T)1 / inv_cs2;

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Bottom() {
    return {{2, 6, 7, 12, 17, 18, 19, 21, 22}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Top() {
    return {{3, 8, 9, 11, 15, 16, 23, 25, 26}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Front() {
    return {{4, 6, 8, 14, 16, 18, 20, 21, 26}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Back() {
    return {{5, 7, 9, 13, 15, 17, 22, 24, 25}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
    return {
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

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, dimQ> weight() {
    return {(T)8 / (T)27,  (T)2 / (T)27,  (T)1 / (T)54,  (T)1 / (T)54,
            (T)1 / (T)54,  (T)1 / (T)54,  (T)1 / (T)216, (T)1 / (T)216,
            (T)1 / (T)216, (T)1 / (T)216, (T)2 / (T)27,  (T)1 / (T)54,
            (T)1 / (T)54,  (T)1 / (T)54,  (T)1 / (T)54,  (T)1 / (T)216,
            (T)1 / (T)216, (T)1 / (T)216, (T)1 / (T)216, (T)2 / (T)27,
            (T)2 / (T)27,  (T)1 / (T)54,  (T)1 / (T)54,  (T)2 / (T)27,
            (T)2 / (T)27,  (T)1 / (T)54,  (T)1 / (T)54};
  }
};


template <class T>
struct Lattice<T, LatticeType::D3Q33> {
  static constexpr LatticeType Type = LatticeType::D3Q33;

  static constexpr int dimD = 3;
  static constexpr int dimQ = 33;
  static constexpr int faceQ = 10;

  LBM_DEVICE LBM_HOST
  static inline constexpr Position halo() {
    return Position({2, 2, 2});
  }

  static constexpr T cs2 = (T) 0.4156023517935171;
  static constexpr T inv_cs2 = (T)1 / cs2;

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Bottom() {
    return {{2, 6, 7, 13, 18, 19, 21, 23, 24, 30}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Top() {
    return {{3, 8, 9, 12, 16, 17, 25, 27, 28, 29}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Front() {
    return {{4, 6, 8, 15, 17, 19, 22, 23, 28, 32}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<unsigned int, faceQ> iQ_Back() {
    return {{5, 7, 9, 14, 16, 18, 24, 26, 27, 31}};
  }

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
    return {
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
        MathVector<T, dimD>{{(T)-2, (T) 0, (T) 0}},
        MathVector<T, dimD>{{(T) 1, (T) 0, (T) 0}},
        MathVector<T, dimD>{{(T) 1, (T) 1, (T) 0}},
        MathVector<T, dimD>{{(T) 1, (T)-1, (T) 0}},
        MathVector<T, dimD>{{(T) 1, (T) 0, (T) 1}},
        MathVector<T, dimD>{{(T) 1, (T) 0, (T)-1}},
        MathVector<T, dimD>{{(T) 1, (T) 1, (T) 1}},
        MathVector<T, dimD>{{(T) 1, (T) 1, (T)-1}},
        MathVector<T, dimD>{{(T) 1, (T)-1, (T) 1}},
        MathVector<T, dimD>{{(T) 1, (T)-1, (T)-1}},
        MathVector<T, dimD>{{(T) 2, (T) 0, (T) 0}},
        MathVector<T, dimD>{{(T) 0, (T)-1, (T) 0}},
        MathVector<T, dimD>{{(T) 0, (T) 0, (T)-1}},
        MathVector<T, dimD>{{(T) 0, (T)-1, (T)-1}},
        MathVector<T, dimD>{{(T) 0, (T)-1, (T) 1}},
        MathVector<T, dimD>{{(T) 0, (T) 1, (T) 0}},
        MathVector<T, dimD>{{(T) 0, (T) 0, (T) 1}},
        MathVector<T, dimD>{{(T) 0, (T) 1, (T) 1}},
        MathVector<T, dimD>{{(T) 0, (T) 1, (T)-1}},
        MathVector<T, dimD>{{(T) 0, (T) 2, (T) 0}},
        MathVector<T, dimD>{{(T) 0, (T)-2, (T) 0}},
        MathVector<T, dimD>{{(T) 0, (T) 0, (T) 2}},
        MathVector<T, dimD>{{(T) 0, (T) 0, (T)-2}},
      };
  }

  static constexpr T w0 = (T) 0.177627658370520295649084;
  static constexpr T w1 = (T) 0.103315974899246818673111;
  static constexpr T w2 = (T) 0.000513472406731114352456;
  static constexpr T w3 = (T) 0.021333928148672240120078;
  static constexpr T w4 = (T) 0.004273899693974583187026;

  LBM_DEVICE LBM_HOST
  static inline constexpr MathVector<T, dimQ> weight() {
    return {w0, w1, w2, w2, w2, w2, w3, w3, w3, w3, w4, w1, w2, w2, w2, w2,
            w3, w3, w3, w3, w4, w1, w1, w2, w2, w1, w1, w2, w2, w4, w4, w4, w4};
  }
 };


typedef Lattice<dataT, latticeT> L;
typedef Lattice<unsigned int, latticeT> uiL;

}  // namespace lbm
