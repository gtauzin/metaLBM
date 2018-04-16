#pragma once

#include "Commons.h"
#include "Options.h"

namespace lbm {

template <Architecture architecture>
class Event {};

template <>
class Event<Architecture::Generic> {
 protected:
 public:
};

template <>
class Event<Architecture::GPU> : public Event<Architecture::Generic> {
 private:
 public:
  Event() {}

  ~Event() {}

  void synchronize() {}

  void record(Stream<Architecture::CPU> stream) {
    CUDA_CALL(cudaEventRecord(event, stream));
  }

  void wait(Stream<Architecture::CPU> stream) {
    CUDA_CALL(cudaEventWaitEvent(stream, event, 0));
  }
};

}  // namespace lbm
