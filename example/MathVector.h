#ifndef MATHVECTOR_H
#define MATHVECTOR_H

namespace lbm {

  template<class U, unsigned int NumberComponents>
  class MathVector {
  public:
    U sArray[NumberComponents];

    __device__ __host__
    U& operator[] (int i) {
      return sArray[i];
    }

    __device__ __host__
    const U& operator[] (int i) const {
      return sArray[i];
    }

  };
}

#endif // MATHVECTOR_H
