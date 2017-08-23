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
  private:
    U * RESTRICT dArray;
    unsigned int numberElements;

  public:
    enum exception {MEMFAIL, INVALIDSIZE};

    DynamicArray()
      : numberElements(0)
    {
      dArray = (U*)malloc(numberElements*sizeof(U));

      if (dArray == NULL)
        throw MEMFAIL;
    }

    DynamicArray(const unsigned int numberElements_in,
                 const U& value_in = (U) 0)
      : numberElements(numberElements_in)
    {
      dArray = (U*)malloc(numberElements*sizeof(U));

      if (dArray == NULL)
        throw MEMFAIL;

      for(unsigned int i = 0; i < numberElements; ++i) dArray[i] = value_in;
    }


  DynamicArray(const DynamicArray& dArray_in)
      : numberElements(dArray_in.size())
    {
      dArray = (U*)malloc(dArray_in.size()*sizeof(U));

      if (dArray == NULL)
        throw MEMFAIL;

      memcpy(dArray, dArray_in.data(), dArray_in.size()*sizeof(U));
    }

    ~DynamicArray(){
      if(dArray) {
        free(dArray);
        dArray = NULL;
      }
    }

    U& operator[] (int i) {
      return dArray[i];
    }

    const U& operator[] (int i) const {
      return dArray[i];
    }

    DynamicArray<U, Architecture::Generic>& operator=(const DynamicArray& other) {
      if (this == &other)
        return *this;

      if (other.size() == 0)
        clear();

      resize(other.size());
      memcpy(dArray, other.data(), other.size()*sizeof(U));
      return *this;
    }

    U * RESTRICT data() {
      return dArray;
    }

    const U * RESTRICT data() const {
      return dArray;
    }


    const unsigned int size() {
      return numberElements;
    }

    const unsigned int size() const {
      return numberElements;
    }

    void resize(unsigned int numberElements_in) {
      numberElements = numberElements_in;

      if (numberElements != 0)
        {
          dArray = (U*)realloc(dArray, numberElements*sizeof(U));

          if (dArray == NULL)
            throw MEMFAIL;
        }
      else {
        clear();
      }
    }

    void swap(DynamicArray<U, Architecture::Generic>& other) {
      std::swap(*this, other);
      /* U * RESTRICT temporary = other.data(); */
      /* other.data() = data(); */
      /* dArray = temporary; */
    }

    void copy(DynamicArray other) {
      memcpy(dArray, other.data(), other.size()*sizeof(U));
    }

    void clear() {
      numberElements = 0;
      dArray = (U*)realloc(dArray, numberElements*sizeof(U));
    }
  };


  template<class U>
  class DynamicArray<U, Architecture::CPU>
    :public DynamicArray<U, Architecture::Generic>{
  private:
    U * RESTRICT dArray;
    unsigned int numberElements;

  public:
    enum exception {MEMFAIL, INVALIDSIZE};

    DynamicArray()
      : numberElements(0)
    {
      dArray = (U*)malloc(numberElements*sizeof(U));

      if (dArray == NULL)
        throw MEMFAIL;
    }

    DynamicArray(const unsigned int numberElements_in,
                 const U& value_in = (U) 0)
      : numberElements(numberElements_in)
    {
      dArray = (U*)malloc(numberElements*sizeof(U));

      if (dArray == NULL)
        throw MEMFAIL;

      for(unsigned int i = 0; i < numberElements; ++i) dArray[i] = value_in;
    }


  DynamicArray(const DynamicArray& dArray_in)
      : numberElements(dArray_in.size())
    {
      dArray = (U*)malloc(dArray_in.size()*sizeof(U));

      if (dArray == NULL)
        throw MEMFAIL;

      memcpy(dArray, dArray_in.data(), dArray_in.size()*sizeof(U));
    }

    ~DynamicArray(){
      if(dArray) {
        free(dArray);
        dArray = NULL;
      }
    }

    U& operator[] (int i) {
      return dArray[i];
    }

    const U& operator[] (int i) const {
      return dArray[i];
    }

    DynamicArray<U, Architecture::CPU>& operator=(const DynamicArray& other) {
      if (this == &other)
        return *this;

      if (other.size() == 0)
        clear();

      resize(other.size());
      memcpy(dArray, other.data(), other.size()*sizeof(U));
      return *this;
    }

    U * RESTRICT data() {
      return dArray;
    }

    const U * RESTRICT data() const {
      return dArray;
    }


    const unsigned int size() {
      return numberElements;
    }

    const unsigned int size() const {
      return numberElements;
    }

    void resize(unsigned int numberElements_in) {
      numberElements = numberElements_in;

      if (numberElements != 0)
        {
          dArray = (U*)realloc(dArray, numberElements*sizeof(U));

          if (dArray == NULL)
            throw MEMFAIL;
        }
      else {
        clear();
      }
    }

    using DynamicArray<U, Architecture::Generic>::swap;

    void copy(DynamicArray other) {
      memcpy(dArray, other.data(), other.size()*sizeof(U));
    }

    void clear() {
      numberElements = 0;
      dArray = (U*)realloc(dArray, numberElements*sizeof(U));
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
