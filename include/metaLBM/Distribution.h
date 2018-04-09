#pragma once

#include <iostream>

#include "Options.h"
#include "Commons.h"
#include "Field.h"
#include "Domain.h"
#include "Lattice.h"
#include "DynamicArray.cuh"

namespace lbm {

  template <class T, Architecture architecture>
  class Distribution
    : public Field<T, L::dimQ, architecture, true> {
  private:
    using Base = Field<T, L::dimQ, architecture, true>;
    DynamicArray<T, architecture> haloArrayPrevious;
    DynamicArray<T, architecture> haloArrayNext;

  protected:
    using Base::localArray;

  public:
    using Base::fieldName;

    Distribution(const unsigned int numberElements_in)
      : Base("distribution", numberElements_in)
      , haloArrayPrevious(hSD::volume()*L::dimQ)
      , haloArrayNext(hSD::volume()*L::dimQ)
    {}

    using Base::getLocalData;

    DEVICE HOST
    T * getHaloDataPrevious() {
      return haloArrayPrevious.data();
    }

    DEVICE HOST
    T * getHaloDataNext() {
      return haloArrayNext.data();
    }

    DynamicArray<T, architecture>& getHaloArrayPrevious() {
      return haloArrayPrevious;
    }
  };

}
