#ifndef DYNAMICARRAY_H
#define DYNAMICARRAY_H

#include <cstring>
#include <stdlib.h>

#include "Commons.h"
#include "Options.h"

#ifdef USE_FFTW
  #define MALLOC fftw_malloc
  #define FREE fftw_free
#else
  #define MALLOC malloc
  #define FREE free
#endif

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

    HOST DEVICE
    void operator()(const MathVector<unsigned int, 3>& iP, const U value) {
      dArrayPtr[iP[0]] = value;
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
    Computation<Architecture::CPU, 1> computation;

  public:
    using DynamicArray<U, Architecture::Generic>::operator[];
    using DynamicArray<U, Architecture::Generic>::operator();
    using DynamicArray<U, Architecture::Generic>::data;
    using DynamicArray<U, Architecture::Generic>::size;

    DynamicArray()
      : DynamicArray<U, Architecture::Generic>()
      , computation(MathVector<unsigned int, 3>{{0}},
                    MathVector<unsigned int, 3>{{numberElements}})
    {}

    DynamicArray(const unsigned int numberElements_in,
                 const U& value_in = (U) 0)
      : DynamicArray<U, Architecture::Generic>(numberElements_in)
      , computation(MathVector<unsigned int, 3>({0, 0, 0}),
                    MathVector<unsigned int, 3>({numberElements, 0, 0}))

    {
      dArrayPtr = (U*)MALLOC(numberElements_in*sizeof(U));

      computation.Do(*this, value_in);
    }

  DynamicArray(const DynamicArray<U, Architecture::CPU>& dArray_in)
    : DynamicArray<U, Architecture::Generic>(dArray_in.size())
    , computation(MathVector<unsigned int, 3>{{0}},
                  MathVector<unsigned int, 3>{{numberElements}})

    {
      dArrayPtr = (U*)malloc(dArray_in.size()*sizeof(U));
      copyFrom(dArray_in);
    }

    ~DynamicArray() {
      if(dArrayPtr) {
        FREE(dArrayPtr);
        dArrayPtr = NULL;
      }
    }

    void copyFrom(const DynamicArray<U, Architecture::CPU>& other) {
      memcpy(dArrayPtr, other.data(), other.size()*sizeof(U));
    }

    void copyTo(DynamicArray<U, Architecture::CPU>& other) const {
      memcpy(other.data(), dArrayPtr, numberElements*sizeof(U));
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
