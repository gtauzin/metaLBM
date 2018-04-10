#pragma once

#include "Stream.h"
#include "Commons.h"
#include "Options.h"

namespace lbm {

  template<>
  class Stream<Architecture::GPU>
    : public Stream<Architecture::Generic> {
  private:
    cudaStream_t stream;

public:
    Stream() {
      CUDA_CALL( cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) );
    }

    Stream(const cudaStream_t& stream_in)
      : stream(stream_in)
    {}

    ~Stream() {
      CUDA_CALL( cudaStreamDestroy(stream) );
    }

    void synchronize() {
      CUDA_CALL( cudaStreamSynchronize(stream) );
    }

    cudaStream_t get() {
      return stream;
    }

  };

  template<>
  class DefaultStream<Architecture::GPU>
    : public Stream<Architecture::GPU> {
  private:
    using Base = Stream<Architecture::GPU>;

  public:
    DefaultStream()
      : Base((cudaStream_t) 0)
    {}

    using Base::get;

  };
}
