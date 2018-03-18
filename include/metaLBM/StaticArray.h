#ifndef STATICARRAY_H
#define STATICARRAY_H


#include "Commons.h"
#include "Helpers.h"

namespace lbm {

  template<class U, unsigned int Size>
  class StaticArray {
  public:
    U sArray[Size];

    DEVICE HOST
    StaticArray<U, Size>& operator=(const StaticArray<U, Size> other){
      for(auto i = 0; i < Size; ++i) {
        sArray[i] = other[i];
      }

      return *this;
    }

    DEVICE HOST
    StaticArray<U, Size>& operator=(const U other[Size]){
      for(auto i = 0; i < Size; ++i) {
          sArray[i] = other[i];
      }

      return *this;
    }

    DEVICE HOST
    U& operator[] (int i) {
      return sArray[i];
    }

    DEVICE HOST
    const U& operator[] (int i) const {
      return sArray[i];
    }

    DEVICE HOST
    U * data() {
      return sArray;
    }

    DEVICE HOST
    const U * data() const {
      return sArray;
    }


    DEVICE HOST
    static constexpr unsigned int size() {
      return Size;
    }
  };

  template<class U, unsigned int Size>
  std::ostream& operator<<(std::ostream& os, const StaticArray<U, Size>& sArray){
    os << "[";
    for(auto i = 0; i < Size-1; ++i) {
      os << sArray[i] << " ";
    }
    os << sArray[Size-1] << "]";
    return os;
  }

  template<class U, unsigned int Size>
  DEVICE HOST
  bool operator==(StaticArray<U, Size> const &lhs,
                  StaticArray<U, Size> const &rhs) {
    for(auto i = 0; i < Size; ++i){
      if (!(lhs[i] == rhs[i])) {
        return false;
      }
    }
    return true;
  }

}

#endif // STATICARRAY_H
