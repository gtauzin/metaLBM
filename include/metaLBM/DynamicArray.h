#ifndef DYNAMICARRAY_H
#define DYNAMICARRAY_H

#include <cstring>
#include <stdlib.h>

#include "Commons.h"
#include "Options.h"

namespace lbm {

  template<class U, Architecture architecture>
  class DynamicArray {};


  template<class U>
  class DynamicArray<U, Architecture::Generic> {
  protected:
    unsigned int numberElements;
    U * RESTRICT dArrayPtr;

    DynamicArray(const unsigned int numberElements_in = 0,
                 U * dArrayPtr_in = NULL)
      : numberElements(numberElements_in)
      , dArrayPtr(dArrayPtr_in)
    {}

  public:
    DEVICE HOST
    U& operator[] (int i) {
      return dArrayPtr[i];
    }

    DEVICE HOST
    const U& operator[] (int i) const {
      return dArrayPtr[i];
    }

    DEVICE HOST
    U * RESTRICT data() {
      return dArrayPtr;
    }

    DEVICE HOST
    const U * RESTRICT data() const {
      return dArrayPtr;
    }

    DEVICE HOST
    void swap(DynamicArray& other) {
      DynamicArray temp = *this;
      *this = other;
      other = temp;
      //std::swap(*this, other);
    }

    DEVICE HOST
    unsigned int size() {
      return numberElements;
    }

    DEVICE HOST
    unsigned int size() const {
      return numberElements;
    }

  };


  template<class U>
  class DynamicArray<U, Architecture::CPU>
    : public DynamicArray<U, Architecture::Generic> {
  private:
    using DynamicArray<U, Architecture::Generic>::numberElements;
    using DynamicArray<U, Architecture::Generic>::dArrayPtr;

  public:
    using DynamicArray<U, Architecture::Generic>::operator[];
    using DynamicArray<U, Architecture::Generic>::data;
    using DynamicArray<U, Architecture::Generic>::size;
    using DynamicArray<U, Architecture::Generic>::swap;

    DynamicArray()
      : DynamicArray<U, Architecture::Generic>()
    {}

    DynamicArray(const unsigned int numberElements_in,
                 const U& value_in = (U) 0)
      : DynamicArray<U, Architecture::Generic>(numberElements_in)
    {
      dArrayPtr = (U*)malloc(numberElements*sizeof(U));
      for(unsigned int i = 0; i < numberElements; ++i) dArrayPtr[i] = value_in;
    }

  DynamicArray(const DynamicArray<U, Architecture::CPU>& dArray_in)
      : DynamicArray<U, Architecture::Generic>(dArray_in.size())
    {
      dArrayPtr = (U*)malloc(dArray_in.size()*sizeof(U));
      copyFrom(dArray_in);
    }

    ~DynamicArray(){
      if(dArrayPtr) {
        free(dArrayPtr);
        dArrayPtr = NULL;
      }
    }

    void copyFrom(const DynamicArray<U, Architecture::CPU>& other) {
      memcpy(dArrayPtr, other.data(), other.size()*sizeof(U));
    }

    void copyTo(DynamicArray<U, Architecture::CPU>& other) const {
      memcpy(other.data(), dArrayPtr, numberElements*sizeof(U));
    }

    void clear() {
      numberElements = 0;
      dArrayPtr = (U*)realloc(dArrayPtr, numberElements*sizeof(U));
    }
  };

  template<class U, Architecture architecture>
  bool operator==(DynamicArray<U, architecture> const &lhs,
                  DynamicArray<U, architecture> const &rhs) {
    if(lhs.size() == rhs.size) {
      for(unsigned int i = 0; i < lhs.size(); ++i){
        if (!(lhs[i] == rhs[i])) {
          return false;
        }
      }
      return true;
    }
    else {
      return false;
    }
  }

}

#endif // DYNAMICARRAY_H
