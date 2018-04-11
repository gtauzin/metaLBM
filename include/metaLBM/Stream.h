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
    Stream(bool isDefault_in) {
    }

    ~Stream() {
    }

    void synchronize() {
    }

  };

}
