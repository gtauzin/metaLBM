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
    Stream(bool isDefault_in = true) {
      if(isDefault_in) {
	CUDA_CALL( cudaStreamCreateWithFlags(&stream, cudaStreamDefault) );
      }
      else {
	CUDA_CALL( cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) )
      }
    }

    ~Stream() {
      CUDA_CALL( cudaStreamDestroy(stream) );
    }

    void synchronize() {
      CUDA_CALL( cudaStreamSynchronize(stream) );
    }

    cudaStream_t get() {
      return stream;
    }

    cudaStream_t get() const {
      return stream;
    }

  };
}
