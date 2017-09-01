#ifndef STATICARRAY_H
#define STATICARRAY_H


#include "Commons.h"
#include "Helpers.h"

namespace lbm {

  template<class U, unsigned int Size>
  class StaticArray {
  public:
    U sArray[Size];

    #pragma omp declare simd
    DEVICE HOST
    StaticArray<U, Size>& operator=(const StaticArray<U, Size> other){
      for(unsigned int i = 0; i < Size; ++i) {
        sArray[i] = other[i];
      }

      return *this;
    }

    #pragma omp declare simd
    DEVICE HOST
    StaticArray<U, Size>& operator=(const U other[Size]){
      for(unsigned int i = 0; i < Size; ++i) {
          sArray[i] = other[i];
      }

      return *this;
    }

    #pragma omp declare simd
    DEVICE HOST
    U& operator[] (int i) {
      return sArray[i];
    }

    #pragma omp declare simd
    DEVICE HOST
    const U& operator[] (int i) const {
      return sArray[i];
    }

    #pragma omp declare simd
    DEVICE HOST
    static constexpr unsigned int size() {
      return Size;
    }
  };

  template<class U, unsigned int Size>
  std::ostream& operator<<(std::ostream& os, const StaticArray<U, Size>& sArray){
    os << "[";
    for(unsigned int i = 0; i < Size-1; ++i) {
      os << sArray[i] << " ";
    }
    os << sArray[Size-1] << "]";
    return os;
  }

  template<class U, unsigned int Size>
  #pragma omp declare simd
  DEVICE HOST
  bool operator==(StaticArray<U, Size> const &lhs,
                  StaticArray<U, Size> const &rhs) {
    for(unsigned int i = 0; i < Size; ++i){
      if (!(lhs[i] == rhs[i])) {
        return false;
      }
    }
    return true;
  }

}

#endif // STATICARRAY_H
