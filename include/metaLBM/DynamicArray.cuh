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
    using DynamicArray<U, Architecture::Generic>::numberElements;
    using DynamicArray<U, Architecture::Generic>::dArrayPtr;
    Computation<Architecture::GPU, 1> computation;

  public:
    using DynamicArray<U, Architecture::Generic>::operator[];
    using DynamicArray<U, Architecture::Generic>::operator();
    using DynamicArray<U, Architecture::Generic>::data;
    using DynamicArray<U, Architecture::Generic>::size;

    DynamicArray()
      : DynamicArray<U, Architecture::Generic>()
      , computation(MathVector<unsigned int, 3>{{0}},
                    MathVector<unsigned int, 3>{{0}})
    {}

    DynamicArray(const unsigned int numberElements_in,
                 const U value_in = (U) 0)
      : DynamicArray<U, Architecture::Generic>(numberElements_in)
      , computation(MathVector<unsigned int, 3>{{0}},
                    MathVector<unsigned int, 3>{{numberElements_in}})
    {
      CUDA_CALL( MALLOC_GPU((void**)&dArrayPtr, numberElements_in*sizeof(U)); )
      computation.Do(*this, value_in);
    }

    DynamicArray(const DynamicArray& dArray_in)
      : DynamicArray<U, Architecture::Generic>(dArray_in.size())
      , computation(MathVector<unsigned int, 3>{{0}},
                    MathVector<unsigned int, 3>{{numberElements}})
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
  private:
   using Base = DynamicArray<U, Architecture::CPU>;

    using Base::numberElements;
    using Base::dArrayPtr;
    Computation<Architecture::GPU, 1> computation;

  public:
    using Base::operator[];
    using Base::operator();
    using Base::data;
    using Base::size;

    DynamicArray()
      : Base()
    {}

    DynamicArray(const unsigned int numberElements_in,
                 const U value_in = (U) 0)
      : Base()
      , computation(MathVector<unsigned int, 3>{{0}},
                    MathVector<unsigned int, 3>{{numberElements_in}})
    {
      numberElements = numberElements_in;
      CUDA_CALL( cudaMallocHost((void**)&dArrayPtr, numberElements_in*sizeof(U)); )
      computation.Do(*this, value_in);
    }

    DynamicArray(const DynamicArray& dArray_in)
      : Base()
      , computation(MathVector<unsigned int, 3>{{0}},
                    MathVector<unsigned int, 3>{{numberElements}})
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

} // end namespace lbm

#endif // DYNAMICARRAY_CUH
