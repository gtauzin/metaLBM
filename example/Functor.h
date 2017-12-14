#ifndef FUNCTOR_H
#define FUNCTOR_H
#include <cstdio>

#include "metaLBM/MathVector.h"

namespace lbm {


  struct Functor {
    DynamicArray<double, Architecture::GPU> dArray;
    double * arrayPtr;
    unsigned int numberElements;

  Functor(DynamicArray<double, Architecture::GPU>& dArray_in,
          unsigned int numberElements_in)
  : dArray(dArray_in)
  , arrayPtr(dArray_in.data())
  , numberElements(numberElements_in)
    {}

    __host__ __device__
    void operator()(const MathVector<unsigned int, 3>& iP) {
      printf("Before assignment: iP[0] = %d, dArray[iP[0]] = %f\n", iP[0], dArray[iP[0]]);
      dArray[iP[0]] += (double) iP[0];
      printf("After assignment: iP[0] = %d, dArray[iP[0]] = %f\n", iP[0], dArray[iP[0]]);

      printf("Before assignment: iP[0] = %d, arrayPtr[iP[0]] = %f\n", iP[0], arrayPtr[iP[0]]);
      arrayPtr[iP[0]] += (double) iP[0];
      printf("After assignment: iP[0] = %d, arrayPtr[iP[0]] = %f\n", iP[0], arrayPtr[iP[0]]);

    }

    void compute() {
      Computation<Architecture::GPU, 1>().Do(MathVector<unsigned int, 3>{{0}},
                                             MathVector<unsigned int, 3>{{numberElements}},
                                             MathVector<unsigned int, 3>{{numberElements}},
                                             *this);

      Computation<Architecture::GPU, 1>().Do(MathVector<unsigned int, 3>{{0}},
                                             MathVector<unsigned int, 3>{{numberElements}},
                                             MathVector<unsigned int, 3>{{numberElements}},
                                             *this);



                                             /* [=, *this] DEVICE (MathVector<unsigned int, 3>& iP){ */
                                             /*   printf("iP[0] = %d, value = %f\n", */
                                             /*          iP[0], dArray[iP[0]]); */
                                             /*   dArray[0] = 2; */

                                             /*   printf("iP[0] = %d, value = %f\n", */
                                             /*          iP[0], dArray[iP[0]]); */

                                             /* }); */
    }

  };
}

#endif // FUNCTOR_H
