#ifndef DYNAMICARRAY_CUH
#define DYNAMICARRAY_CUH

#include <cstring>
#include <stdlib.h>

#include "DynamicArray.h"
#include "Commons.h"
#include "Options.h"

namespace lbm {

  template<class U>
  class DynamicArray<U, Architecture::GPU>
    :public DynamicArray<U, Architecture::Generic> {
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
      CUDA_CALL( cudaMalloc(&dArrayPtr, numberElements*sizeof(U)); )
      for(unsigned int i = 0; i < numberElements; ++i) dArrayPtr[i] = value_in;
    }

    DynamicArray(const DynamicArray<U, Architecture::CPU>& dArray_in)
      : DynamicArray<U, Architecture::Generic>(dArray_in.size())
    {
      CUDA_CALL( cudaMalloc(&dArrayPtr, dArray_in.size()*sizeof(U)); )
      copyFrom(dArray_in);
    }


    DynamicArray(const DynamicArray<U, Architecture::GPU>& dArray_in)
      : DynamicArray<U, Architecture::Generic>(dArray_in.size())
    {
      CUDA_CALL( cudaMalloc(&dArrayPtr, numberElements*sizeof(U)); )
      copyFrom(dArray_in);
  }

    ~DynamicArray(){
      if(dArrayPtr) {
        CUDA_CALL( cudaFree(dArrayPtr); )
        dArrayPtr = NULL;
      }
    }

    void copyFrom(const DynamicArray<U, Architecture::CPU>& other) {
      CUDA_CALL( cudaMemcpy(dArrayPtr, other.data(), other.size()*sizeof(U),
                            cudaMemcpyHostToDevice); )
    }

    void copyFrom(const DynamicArray<U, Architecture::GPU>& other) {
      CUDA_CALL( cudaMemcpy(dArrayPtr, other.data(), other.size()*sizeof(U),
                            cudaMemcpyDeviceToDevice); )
    }

    void copyTo(DynamicArray<U, Architecture::CPU>& other) const {
      CUDA_CALL( cudaMemcpy(other.data(), dArrayPtr, numberElements*sizeof(U),
                            cudaMemcpyDeviceToHost); )
    }

    void copyTo(DynamicArray<U, Architecture::GPU>& other) const {
      CUDA_CALL( cudaMemcpy(other.data(), dArrayPtr, numberElements*sizeof(U),
                            cudaMemcpyDeviceToDevice); )
    }

    DEVICE HOST
    void clear() {
      numberElements = 0;
      CUDA_CALL( cudaFree(dArrayPtr) );
      CUDA_CALL_ERROR( cudaMalloc(&dArrayPtr, numberElements*sizeof(U)) );
    }
  };

}

#endif // DYNAMICARRAY_CUH