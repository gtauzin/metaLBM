#pragma once

#include "Stream.h"
#include "Commons.h"
#include "Options.h"

namespace lbm {

  template<>
  class Stream<Architecture::GPU>
    : public Stream<Architecture::Generic> {
  private:
    bool isDefault;
    cudaStream_t stream;

public:
    Stream(bool isDefault_in = true)
      : isDefault(isDefault_in)
    {
      if(isDefault_in) {
	stream = (cudaStream_t) 0;
      }
      else {
	CUDA_CALL( cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) )
      }
    }

    ~Stream() {
      if(!isDefault) {
	CUDA_CALL( cudaStreamDestroy(stream) );
      }
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
