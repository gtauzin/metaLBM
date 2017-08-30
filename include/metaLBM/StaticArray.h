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
      UnrolledFor<0, Size>::Do([&] HOST DEVICE (int i) {
          sArray[i] = other[i];
        });

      return *this;
    }

    #pragma omp declare simd
    DEVICE HOST
    StaticArray<U, Size>& operator=(const U other[Size]){
      UnrolledFor<0, Size>::Do([&] HOST DEVICE (int i) {
          sArray[i] = other[i];
        });

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
  #pragma omp declare simd
  DEVICE HOST
  std::ostream& operator<<(std::ostream& os, const StaticArray<U, Size>& sArray){
    os << "[";
    UnrolledFor<0, Size-1>::Do([&] HOST DEVICE (int i) {
        os << sArray[i] << " ";
      });
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
