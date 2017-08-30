#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <string>
#include <stdlib.h>

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"

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

    #pragma omp declare simd
    DEVICE HOST
    static inline constexpr MathVector<unsigned int, 3> halo() {
      return {1, 0, 0};
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    #pragma omp declare simd
    DEVICE HOST
    static inline constexpr MathVector<MathVector<T, dimD>, dimQ> celerity() {
      return
        {
          MathVector<T, dimD>{{(T)0}},
            MathVector<T, dimD>{{(T)-1}},
              MathVector<T, dimD>{{(T)1}}
        };
    }

    #pragma omp declare simd
    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
      return
        {
          (T)2/(T)3, (T)1/(T)6, (T)1/(T)6
            };
    }
  };


  template <class T>
  struct Lattice<T, LatticeType::D2Q5>
  {
    static constexpr LatticeType Type = LatticeType::D2Q5;
    static constexpr int dimD = 2;
    static constexpr int dimQ = 5;

    #pragma omp declare simd
    DEVICE HOST
    static inline constexpr MathVector<unsigned int, 3> halo() {
      return {1, 1, 0};
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    #pragma omp declare simd
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

    #pragma omp declare simd
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

  template <class T>
  struct Lattice<T, LatticeType::D2Q9>
  {
    static constexpr LatticeType Type = LatticeType::D2Q9;

    static constexpr int dimD = 2;
    static constexpr int dimQ = 9;

    #pragma omp declare simd
    DEVICE HOST
    static inline constexpr MathVector<unsigned int, 3> halo() {
      return {1, 1, 0};
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    #pragma omp declare simd
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

    #pragma omp declare simd
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

  template <class T>
  struct Lattice<T, LatticeType::D3Q15>
  {
    static constexpr LatticeType Type = LatticeType::D3Q15;

    static constexpr int dimD = 3;
    static constexpr int dimQ = 15;

    #pragma omp declare simd
    DEVICE HOST
    static inline constexpr MathVector<unsigned int, 3> halo() {
      return {1, 1, 1};
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    #pragma omp declare simd
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

    #pragma omp declare simd
    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
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
  struct Lattice<T, LatticeType::D3Q19>
  {
    static constexpr LatticeType Type = LatticeType::D3Q19;

    static constexpr int dimD = 3;
    static constexpr int dimQ = 19;

    #pragma omp declare simd
    DEVICE HOST
    static inline constexpr MathVector<unsigned int, 3> halo() {
      return {1, 1, 1};
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    #pragma omp declare simd
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

    #pragma omp declare simd
    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
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
  struct Lattice<T, LatticeType::D3Q27>
  {
    static constexpr LatticeType Type = LatticeType::D3Q27;

    static constexpr int dimD = 3;
    static constexpr int dimQ = 27;

    #pragma omp declare simd
    DEVICE HOST
    static inline constexpr MathVector<unsigned int, 3> halo() {
      return {1, 1, 1};
    }

    static constexpr T inv_cs2 = (T)3;
    static constexpr T cs2 = (T)1/inv_cs2;

    #pragma omp declare simd
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

    #pragma omp declare simd
    DEVICE HOST
    static inline constexpr MathVector<T, dimQ> weight() {
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

  typedef Lattice<dataT, latticeT> L;
  typedef Lattice<unsigned int, latticeT> uiL;

}

#endif // LATTICE_H
