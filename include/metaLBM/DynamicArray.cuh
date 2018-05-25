#pragma once

#include <cstring>
#include <stdlib.h>

#include "DynamicArray.h"
#include "Commons.h"
#include "Options.h"

#ifdef USE_NVSHMEM
  #include<shmem.h>
  #include<shmemx.h>
#endif


namespace lbm {

  template<class U>
  class DynamicArray<U, Architecture::GPU>
    : public DynamicArray<U, Architecture::Generic> {
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
    #ifdef USE_NVSHMEM
      dArrayPtr = (double *) shmem_malloc(numberElements*sizeof(U));
    #else
      LBM_CUDA_CALL( cudaMalloc((void**)&dArrayPtr, numberElements*sizeof(U)); )
    #endif
	}

    DynamicArray(const DynamicArray& dArray_in)
      : Base(dArray_in.size())
    {
    #ifdef USE_NVSHMEM
      dArrayPtr = (double *) shmem_malloc(numberElements*sizeof(U));
    #else
      LBM_CUDA_CALL( cudaMalloc((void**)&dArrayPtr, numberElements*sizeof(U)); )
    #endif

	copyFrom(dArray_in);
    }

    ~DynamicArray(){
      if(dArrayPtr) {
        #ifdef USE_NVSHMEM
          shmem_free(dArrayPtr);
        #else
          LBM_CUDA_CALL( cudaFree(dArrayPtr); )
        #endif

	dArrayPtr = NULL;
      }
    }

    void copyFrom(const DynamicArray<U, Architecture::CPU>& other) {
      LBM_CUDA_CALL( cudaMemcpy(dArrayPtr, other.data(), other.size()*sizeof(U),
				  cudaMemcpyHostToDevice); )
    }

    void copyFrom(const DynamicArray<U, Architecture::GPU>& other) {
      LBM_CUDA_CALL( cudaMemcpy(dArrayPtr, other.data(), other.size()*sizeof(U),
				  cudaMemcpyDeviceToDevice); )
    }

    void copyTo(DynamicArray<U, Architecture::CPU>& other) const {
      LBM_CUDA_CALL( cudaMemcpy(other.data(), dArrayPtr, numberElements*sizeof(U),
				  cudaMemcpyDeviceToHost); )
    }

    void copyTo(DynamicArray<U, Architecture::GPU>& other) const {
      LBM_CUDA_CALL( cudaMemcpy(other.data(), dArrayPtr, numberElements*sizeof(U),
				  cudaMemcpyDeviceToDevice); )
    }

  };

  template<class U>
  class DynamicArray<U, Architecture::CPUPinned>
    : public DynamicArray<U, Architecture::CPU> {
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
      LBM_CUDA_CALL( cudaMallocHost((void**)&dArrayPtr, numberElements*sizeof(U)); )
    }

    DynamicArray(const DynamicArray& dArray_in)
      : Base()
    {
      numberElements = dArray_in.size();
      LBM_CUDA_CALL( cudaMallocHost((void**)&dArrayPtr, numberElements*sizeof(U)); )
      copyFrom(dArray_in);
    }

    ~DynamicArray(){
      if(dArrayPtr) {
        LBM_CUDA_CALL( cudaFreeHost(dArrayPtr); )
        dArrayPtr = NULL;
      }
    }

    using Base::copyFrom;
    using Base::copyTo;

  }; // end class DynamicArray<U, Architecture::CPUPinned>

} // end namespace lbm
