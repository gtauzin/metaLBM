#pragma once

#include "Commons.h"
#include "Options.h"
#include "Stream.h"

namespace lbm {

  template <>
  class Stream<Architecture::GPU> : public Stream<Architecture::Generic> {
  private:
    cudaStream_t stream;

  public:
    Stream(bool isDefault_in = true) {
      if (isDefault_in) {
        LBM_CUDA_CALL(cudaStreamCreateWithFlags(&stream, cudaStreamDefault));
      } else {
        LBM_CUDA_CALL(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
      }
    }

    ~Stream() { LBM_CUDA_CALL(cudaStreamDestroy(stream)); }

    void synchronize() { LBM_CUDA_CALL(cudaStreamSynchronize(stream)); }

    cudaStream_t get() { return stream; }

    cudaStream_t get() const { return stream; }
  };

  template <>
  class Stream<Architecture::GPU_SHMEM> : public Stream<Architecture::GPU> {
  private:
    using Base = Stream<Architecture::GPU>;

  public:
    Stream(bool isDefault_in = true)
      : Base(isDefault_in)
    {}

    ~Stream() {
      Base::~Stream();
    }

    using Base::synchronize;
    using Base::get;
  };


}  // namespace lbm
