#pragma once

#include "Commons.h"
#include "Event.h"
#include "Options.h"

namespace lbm {

  template <>
  class Event<Architecture::GPU> : public Event<Architecture::Generic> {
  protected:
    cudaEvent_t event;

  public:
    Event() { LBM_CUDA_CALL(cudaEventCreate(&event)); }

    ~Event() { LBM_CUDA_CALL(cudaEventDestroy(event)); }

    void record(Stream<Architecture::GPU>& stream) {
      LBM_CUDA_CALL(cudaEventRecord(event, stream.get()));
    }

    void wait(Stream<Architecture::GPU>& stream) {
      LBM_CUDA_CALL(cudaStreamWaitEvent(stream.get(), event, 0));
    }
  };

  template <>
  class Event<Architecture::GPU_SHMEM> : public Event<Architecture::GPU> {
  private:
    using Base = Event<Architecture::GPU>;

  public:
    Event()
      : Base()
    {}

    ~Event() {
      Base::~Event();
    }

    using Base::record;
    using Base::wait;
  };


}  // namespace lbm
