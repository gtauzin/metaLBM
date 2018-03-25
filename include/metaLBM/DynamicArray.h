#ifndef DYNAMICARRAY_H
#define DYNAMICARRAY_H

#include <cstring>
#include <stdlib.h>

#include "Commons.h"
#include "Options.h"
#include "StaticArray.h"

#ifdef USE_FFTW
  #include <fftw3-mpi.h>
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
    U * dArrayPtr;

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
    U * data(const unsigned offset = 0) {
      return dArrayPtr + offset;
    }

    DEVICE HOST
    const U * data(const unsigned offset = 0) const {
      return dArrayPtr + offset;
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

  public:
    using Base::operator[];
    using Base::data;
    using Base::size;

    DynamicArray()
      : Base()
    {}

    DynamicArray(const unsigned int numberElements_in)
      : Base(numberElements_in)
    {
      dArrayPtr = (U*)MALLOC_CPU(numberElements_in*sizeof(U));
    }

  DynamicArray(const DynamicArray<U, Architecture::CPU>& dArray_in)
    : Base(dArray_in.size())
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
      for(auto i = 0; i < lhs.size(); ++i){
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


  template<class U, Architecture architecture,
           unsigned int NumberComponents>
  class MultiDynamicArray {};


  template<class U, unsigned int NumberComponents>
  class MultiDynamicArray<U, Architecture::CPU, NumberComponents>
    : public DynamicArray<U, Architecture::CPU> {
  private:
    using Base = DynamicArray<U, Architecture::CPU>;

  protected:
    using Base::dArrayPtr;
    DynamicArray<U *, Architecture::CPU> sMultiArray;
    unsigned int numberElements;

  public:
    MultiDynamicArray(const unsigned int numberElements_in)
      : Base(NumberComponents*numberElements_in)
      , sMultiArray(NumberComponents)
      , numberElements(numberElements_in)
    {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        sMultiArray[iC] = Base::data(iC*numberElements);
      }
    }

    void copyFrom(const MultiDynamicArray& other) {
      Base::copyFrom(other);
    }

    void copyTo(MultiDynamicArray& other) const {
      Base::copyTo(other);
    }

    DEVICE HOST
    inline unsigned int getNumberElements() {
      return numberElements;
    }

    DEVICE HOST
    U * operator[] (int iC) {
      return sMultiArray[iC];
    }

    DEVICE HOST
    const U * operator[] (int iC) const {
      return sMultiArray[iC];
    }

    DEVICE HOST
    U ** multiData() {
      return sMultiArray.data();
    }

  };


  /* template<class U, unsigned int NumberComponents> */
  /* class MultiDynamicArray<U, Architecture::CPU, NumberComponents> */
  /*   : public DynamicArray<U, Architecture::CPU> { */
  /* private: */
  /*   using Base = DynamicArray<U, Architecture::CPU>; */

  /* protected: */
  /*   using Base::dArrayPtr; */
  /*   U * sMultiArrayPtr[NumberComponents]; */
  /*   unsigned int numberElements; */

  /* public: */
  /*   MultiDynamicArray(const unsigned int numberElements_in) */
  /*     : Base(NumberComponents*numberElements_in) */
  /*     , numberElements(numberElements_in) */
  /*   { */
  /*     for(auto iC = 0; iC < NumberComponents; ++iC) { */
  /*       sMultiArrayPtr[iC] = Base::data(iC*numberElements); */
  /*     } */
  /*   } */

  /*   MultiDynamicArray(const MultiDynamicArray& multiArray_in) */
  /*     : Base(multiArray_in) */
  /*     , numberElements(multiArray_in.numberElements) */
  /*   { */
  /*     for(auto iC = 0; iC < NumberComponents; ++iC) { */
  /*       sMultiArrayPtr[iC] = Base::data(iC*numberElements); */
  /*     } */
  /*   } */

  /*   void copyFrom(const MultiDynamicArray& other) { */
  /*     Base::copyFrom(other); */
  /*   } */

  /*   void copyTo(MultiDynamicArray& other) const { */
  /*     Base::copyTo(other); */
  /*   } */

  /*   DEVICE HOST */
  /*   inline unsigned int getNumberElements() { */
  /*     return numberElements; */
  /*   } */

  /*   DEVICE HOST */
  /*   U * operator[] (int iC) { */
  /*     return sMultiArrayPtr[iC]; */
  /*   } */

  /*   DEVICE HOST */
  /*   const U * operator[] (int iC) const { */
  /*     return sMultiArrayPtr[iC]; */
  /*   } */

  /*   DEVICE HOST */
  /*   U ** multiData() { */
  /*     return sMultiArrayPtr; */
  /*   } */

  /* }; */


}

#endif // DYNAMICARRAY_H
