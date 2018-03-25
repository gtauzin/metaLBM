#ifndef DYNAMICARRAY_CUH
#define DYNAMICARRAY_CUH

#include <cstring>
#include <stdlib.h>

#include "DynamicArray.h"
#include "Commons.h"
#include "Options.h"

#ifdef USE_NVSHMEM
  #include<shmem.h>
  #include<shmemx.h>

  #define MALLOC_GPU shmem_malloc
  #define FREE_GPU shmem_free
#else
  #define MALLOC_GPU cudaMalloc
  #define FREE_GPU cudaFree
#endif


namespace lbm {

  template<class U>
  class DynamicArray<U, Architecture::GPU>
    :public DynamicArray<U, Architecture::Generic> {
  private:
    using Base = DynamicArray<U, Architecture::Generic>;

    using Base::numberElements;
    using Base::dArrayPtr;

  public:
    using Base::operator[];
    using Base::data;
    using Base::size;

    DynamicArray(const unsigned int numberElements_in)
      : Base(numberElements_in)
    {
      CUDA_CALL( MALLOC_GPU((void**)&dArrayPtr, numberElements*sizeof(U)); )
    }

    DynamicArray(const DynamicArray& dArray_in)
      : Base(dArray_in.size())
    {
      CUDA_CALL( MALLOC_GPU((void**)&dArrayPtr, numberElements*sizeof(U)); )
      copyFrom(dArray_in);
    }

    ~DynamicArray(){
      if(dArrayPtr) {
        CUDA_CALL( FREE_GPU(dArrayPtr); )
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

  };

 template<class U>
  class DynamicArray<U, Architecture::CPUPinned>
    :public DynamicArray<U, Architecture::CPU> {
  protected:
   using Base = DynamicArray<U, Architecture::CPU>;

    using Base::numberElements;
    using Base::dArrayPtr;

  public:
    using Base::operator[];
    using Base::data;
    using Base::size;

    DynamicArray(const unsigned int numberElements_in)
      : Base()
    {
      numberElements = numberElements_in;
      CUDA_CALL( cudaMallocHost((void**)&dArrayPtr, numberElements*sizeof(U)); )
    }

    DynamicArray(const DynamicArray& dArray_in)
      : Base()
    {
      numberElements = dArray_in.size();
      CUDA_CALL( cudaMallocHost((void**)&dArrayPtr, numberElements*sizeof(U)); )
      copyFrom(dArray_in);
    }

    ~DynamicArray(){
      if(dArrayPtr) {
        CUDA_CALL( cudaFreeHost(dArrayPtr); )
	dArrayPtr = NULL;
      }
    }

    using Base::copyFrom;
    using Base::copyTo;

  }; // end class DynamicArray<U, Architecture::CPUPinned>


  template<class U, unsigned int NumberComponents>
  class MultiDynamicArray<U, Architecture::CPUPinned, NumberComponents>
    : public DynamicArray<U, Architecture::CPUPinned> {
  private:
    using Base = DynamicArray<U, Architecture::CPUPinned>;

  protected:
    using Base::dArrayPtr;
    U * sMultiArrayPtr[NumberComponents];
    unsigned int numberElements;

  public:
    MultiDynamicArray(const unsigned int numberElements_in)
      : Base(NumberComponents*numberElements_in)
      , numberElements(numberElements_in)
    {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        sMultiArrayPtr[iC] = Base::data(iC*numberElements);
      }
    }

    MultiDynamicArray(const MultiDynamicArray& multiArray_in)
      : Base(multiArray_in)
      , numberElements(multiArray_in.numberElements)
    {
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        sMultiArrayPtr[iC] = Base::data(iC*numberElements);
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
      return sMultiArrayPtr[iC];
    }

    DEVICE HOST
    const U * operator[] (int iC) const {
      return sMultiArrayPtr[iC];
    }

    DEVICE HOST
    U ** multiData() {
      return sMultiArrayPtr;
    }

  };



} // end namespace lbm

#endif // DYNAMICARRAY_CUH
