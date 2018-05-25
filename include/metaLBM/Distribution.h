#pragma once

#include <iostream>

#include "Commons.h"
#include "Domain.h"
#include "DynamicArray.cuh"
#include "Field.h"
#include "Lattice.h"
#include "Options.h"

namespace lbm {

template <class T, Architecture architecture>
class Distribution : public Field<T, L::dimQ, architecture, true> {
 private:
  using Base = Field<T, L::dimQ, architecture, true>;
  DynamicArray<T, architecture> haloArrayPrevious;
  DynamicArray<T, architecture> haloArrayNext;

 protected:
  using Base::localArray;

 public:
  using Base::fieldName;

  Distribution()
      : Base("distribution"),
        haloArrayPrevious(hSD::volume() * L::dimQ),
        haloArrayNext(hSD::volume() * L::dimQ) {}

  using Base::getLocalData;

  LBM_DEVICE LBM_HOST T* getHaloDataPrevious() {
    return haloArrayPrevious.data();
  }

  LBM_DEVICE LBM_HOST T* getHaloDataNext() { return haloArrayNext.data(); }

  DynamicArray<T, architecture>& getHaloArrayPrevious() {
    return haloArrayPrevious;
  }
};

}  // namespace lbm
