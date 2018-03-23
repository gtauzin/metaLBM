#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <iostream>

#include "Options.h"
#include "Commons.h"
#include "Field.h"
#include "Domain.h"
#include "Lattice.h"
#include "DynamicArray.cuh"

namespace lbm {

  template <class T, DomainType initDomainType, Architecture architecture>
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

    Distribution(const std::string& fieldName_in,
                 const unsigned int numberElements_in)
      : Base(fieldName_in, numberElements_in)
      , haloArrayPrevious(hSD::volume()*L::dimQ)
      , haloArrayNext(hSD::volume()*L::dimQ)
    {}

    Distribution(const std::string& fieldName_in,
                 const MultiDynamicArray<T, Architecture::CPU, L::dimQ>& initArray_in)
      : Base(fieldName_in, initArray_in)
      , haloArrayPrevious(hSD::volume()*L::dimQ)
      , haloArrayNext(hSD::volume()*L::dimQ)
    {}

    using Base::getLocalData;

    DEVICE HOST
    T * getHaloDataPrevious(const unsigned int offset = 0) {
      return haloArrayPrevious.data(offset);
    }

    DEVICE HOST
    T * getHaloDataNext(const unsigned int offset = 0) {
      return haloArrayNext.data(offset);
    }

    DynamicArray<T, architecture>& getHaloArrayPrevious() {
      return haloArrayPrevious;
    }
  };

}

#endif // DISTRIBUTION_H
