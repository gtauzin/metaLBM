#pragma once

#include "Commons.h"
#include "Options.h"

namespace lbm {

  template<Architecture architecture>
  class Stream {};

  template<>
  class Stream<Architecture::Generic> {
  protected:

  public:

  };

  template<>
  class Stream<Architecture::CPU>
    : public Stream<Architecture::Generic> {
  private:

  public:
    Stream() {
    }

    ~Stream() {
    }

    void synchronize() {
    }

  };

  template<Architecture architecture>
  class DefaultStream {};

  template<>
  class DefaultStream<Architecture::CPU>
    : public Stream<Architecture::CPU> {
  private:
    using Base = Stream<Architecture::CPU>;

  public:
    DefaultStream()
      : Base()
    {}

  };


}
