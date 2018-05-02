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
class Event<Architecture::CPU> : public Event<Architecture::Generic> {
 private:
 public:
  Event() {}

  ~Event() {}

  void synchronize() {}

  void record(Stream<Architecture::CPU> stream) {
  }

  void wait(Stream<Architecture::CPU> stream) {
  }
};

}  // namespace lbm
