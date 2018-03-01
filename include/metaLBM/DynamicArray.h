#ifndef DYNAMICARRAY_H
#define DYNAMICARRAY_H

#include <cstring>
#include <stdlib.h>

#include "Commons.h"
#include "Options.h"

#ifdef USE_FFTW
  #define MALLOC_CPU fftw_malloc
  #define FREE_CPU fftw_free
#else
  #define MALLOC_CPU malloc
  #define FREE_CPU free
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
    void operator()(const MathVector<unsigned int, 3>& iP,
                    U * RESTRICT array,  const U value) {
      array[iP[0]] = value;
      std::cout << "Allocated ptr base: " << dArrayPtr << std::endl;
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
    using Base = DynamicArray<U, Architecture::Generic>;

  protected:
    using Base::numberElements;
    using Base::dArrayPtr;
    Computation<Architecture::CPU, 1> computation;

  public:
    using Base::operator[];
    using Base::operator();
    using Base::data;
    using Base::size;

    DynamicArray()
      : Base()
      , computation(MathVector<unsigned int, 3>{{0}},
                    MathVector<unsigned int, 3>{{0}})
    {}

    DynamicArray(const unsigned int numberElements_in,
                 const U& value_in = (U) 0)
      : Base(numberElements_in)
      , computation(MathVector<unsigned int, 3>({0, 0, 0}),
                    MathVector<unsigned int, 3>({numberElements, 0, 0}))

    {
      dArrayPtr = (U*)MALLOC_CPU(numberElements_in*sizeof(U));
      std::cout << "Allocated ptr:      " << dArrayPtr << std::endl;
      computation.Do(*this, dArrayPtr, value_in);
      std::cout << "Allocated ptr from: " << Base::dArrayPtr << std::endl;
      std::cout << "Allocated ptr data: " << this->data() << std::endl;
    }

  DynamicArray(const DynamicArray<U, Architecture::CPU>& dArray_in)
    : Base(dArray_in.size())
    , computation(MathVector<unsigned int, 3>{{0}},
                  MathVector<unsigned int, 3>{{numberElements}})

    {
      dArrayPtr = (U*)MALLOC_CPU(dArray_in.size()*sizeof(U));
      copyFrom(dArray_in);
    }

    ~DynamicArray() {
      if(dArrayPtr) {
        FREE_CPU(dArrayPtr);
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
