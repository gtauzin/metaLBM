#ifndef STATICARRAY_H
#define STATICARRAY_H

#include <Helpers.h>

namespace lbm {

  template<class U, unsigned int Size>
  class StaticArray {
  public:
    U sArray[Size];

    StaticArray<U, Size>& operator=(const StaticArray<U, Size> other){
      UnrolledFor<0, Size>::Do([&] (int i) {
          sArray[i] = other[i];
        });

      return *this;
    }

    StaticArray<U, Size>& operator=(const U other[Size]){
      UnrolledFor<0, Size>::Do([&] (int i) {
          sArray[i] = other[i];
        });

      return *this;
    }

    U& operator[] (int i) {
      return sArray[i];
    }

    //constexpr??
    const U& operator[] (int i) const {
      return sArray[i];
    }

    /* constexpr U* data() const { */
    /*   return this.sArray; */
    /* } */


    constexpr unsigned int size() {
      return Size;
    }


  };

  template<class U, unsigned int Size>
  std::ostream& operator<<(std::ostream& os, const StaticArray<U, Size>& sArray){
    os << "[";
    UnrolledFor<0, Size-1>::Do([&] (int i) {
        os << sArray[i] << " ";
      });
    os << sArray[Size-1] << "]";
    return os;
  }

  template<class U, unsigned int Size>
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
