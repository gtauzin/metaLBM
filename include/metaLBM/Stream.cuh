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
    Stream(bool isDefault_in) {
      if(isDefault_in) {
	CUDA_CALL( cudaStreamCreateWithFlags(&stream, cudaStreamDefault) );
	std::cout << "Stream created default: " << stream << std::endl;
      }
      else {
	CUDA_CALL( cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) )
	std::cout << "Stream created non-default: " << stream << std::endl;
      }
    }

    ~Stream() {
      std::cout << "Stream desctructor: " << stream << std::endl;
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
