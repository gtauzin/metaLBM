#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <memory>

#include "Options.h"
#include "Field.h"
#include "DynamicArray.h"

namespace lbm {

  template <class T, MemoryLayout memoryLayout>
  class Distribution: public Field<T, L::dimQ, true> {
  protected:
    using Field<T, L::dimQ, true>::field;
    LocalizedField<T, L::dimQ, true, DomainType::Halo, memoryLayout> haloField;

  public:
    Distribution(const std::string& fieldName_in)
      : Field<T, L::dimQ, true>(fieldName_in, lengths_in)
      {}

    Distribution(const std::string& fieldName_in,
                 const DynamicArray<T>& dArray_in)
      : Field<T, L::dimQ, true>(fieldName_in, field_in)
      {}

  };
}

#endif // DISTRIBUTION_H
