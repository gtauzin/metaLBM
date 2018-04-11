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
    Stream(bool isDefault_in = false) {
      std::cout << "Stream beginning of constructor: " << stream << std::endl;
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
      std::cout << "Stream beginning desctructor: " << stream << std::endl;
      CUDA_CALL( cudaStreamDestroy(stream) );
      std::cout << "Stream end desctructor: " << stream << std::endl;
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
      : Base(true)
    {}

    using Base::get;
  };
}
