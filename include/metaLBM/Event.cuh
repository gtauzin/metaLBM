#pragma once

#include "Event.h"
#include "Commons.h"
#include "Options.h"

namespace lbm {

  template<>
  class Event<Architecture::GPU>
    : public Event<Architecture::Generic> {
  private:
    cudaEvent_t event;

public:
    Event() {
      CUDA_CALL( cudaEventCreate(&event) );
    }

    ~Event() {
      CUDA_CALL( cudaEventDestroy(event) );
    }

    void record(Stream<Architecture::GPU> stream) {
      CUDA_CALL( cudaEventRecord(event, stream) );
    }

    void wait(Stream<Architecture::GPU> stream) {
      CUDA_CALL( cudaEventWaitEvent ( stream, event, 0 ) );
    }

  };

}
