#ifndef DYNAMICARRAY_H
#define DYNAMICARRAY_H

namespace lbm {

  template<class U>
  class DynamicArray {
  protected:
    U* dArray;
    const unsigned int size;

  public:
    enum exception {MEMFAIL, INVALIDSIZE};

    DynamicArray()
      : size(0)
    {
      dArray = (U*)malloc(size*sizeof(U));

      if (dArray == NULL)
        throw MEMFAIL;
    }

    DynamicArray(const DynamicArray<U>& dynamicArray_in)
    {
      dArray = (U*)malloc(dynamicArray_in.size()*sizeof(U));

      if (dArray == NULL)
        throw MEMFAIL;

      memcpy(dArray, dynamicArray.data(), dynamicArray_in.size()*sizeof(U));
      size = dynamicArray.size();
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

    DynamicArray<U>& operator=(const DynamicArray<U>& other) {
      if (this == &other)
        return *this;

      if (other.size() == 0)
        clear();

      resize(other.size());
      memcpy(dArray, other.data(), other.size()*sizeof(U));
      return *this;
    }

    U* data() {
      return this.dArray;
    }

    const U* data() const {
      return this.dArray;
    }


    unsigned int size() {
      return size;
    }

    void resize(unsigned int size_in) {
      size = size_in;

      if (size != 0)
        {
          dArray = (U*)realloc(dArray, size*sizeof(U));

          if (dArray == NULL)
            throw MEMFAIL;
        }
      else {
        clear();
      }
    }

    void clear() {
      size = 0;
      dArray = (U*)realloc(dArray, size*sizeof(U));
    }
  };

  template<class U>
  bool operator==(DynamicArray<U> const &lhs,
                  DynamicArray<U> const &rhs) {
    for(unsigned int i = 0; i < Size; ++i){
      if (!(lhs[i] == rhs[i])) {
        return false;
      }
    }
    return true;
  }

}

#endif // DYNAMICARRAY_H
