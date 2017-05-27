#ifndef FIELD_H
#define FIELD_H

#include <string>
#include <fstream>
#include <ostream>

#include "Options.h"
#include "Domain.h"
#include "DynamicArray.h"
#include "MathVector.h"

namespace lbm {

  /**
   * Field containing data sets to be communicated and dumped.
   *
   * @tparam T datatype.
   * @tparam unsigned int NumberComponents of the field.
   * @tparam bool IsWritten whether the current field should be dumped
   * and therefore allocated.
   */


  template <class T, unsigned int NumberComponents, DomainType domain, bool IsWritten>
  class LocalizedField {};

  template <class T, unsigned int NumberComponents, Domain domain>
  class LocalizedField<T, NumberComponents, true>: public DynamicArray<T> {
  protected:
    using field = DynamicArray<T>::dArray;
    static constexpr unsigned int numberComponents = NumberComponents;
    static constexpr bool isWritten = IsWritten;

  public:
    const std::string fieldName;

    LocalizedField(const std::string& fieldName_in)
      : fieldName(fieldName_in)
      , field(NumberComponents_in*lengths[d::X]*lengths[d::Y]*lengths[d::Z])
      {}

    LocalizedField(const std::string& fieldName_in,
                   const DynamicArray<T>& field_in)
      : fieldName(fieldName_in)
      , field(field_in)
    {}

    T* getRawField(const unsigned int iC) {
      return field[domainVolume * iC];
    }

    T getField(const int idx, const unsigned int iC) {
      return field[iC * domainVolume + idx];
    } // operator overloading????

    MathVector<T, NumberComponents> getVector(const int idx) {
      MathVector<T, NumberComponents> vectorR;
      UnrolledFor<0, NumberComponents>::Do([&] (int iC) {
          vectorR[iC] = getField(idx, iC);
      });

      return vectorR;
    }

    void setField(const int idx, const T value, const unsigned int iC) {
      field[iC * domainVolume + idx] = value;
    }

    void setField(const int idx, const MathVector<T, NumberComponents> value) {
      UnrolledFor<0, NumberComponents>::Do([&] (int iC) {
          setField(idx, value[iC], iC);
      });
    }
  };

  template <class T, unsigned int NumberComponents>
  class LocalizedField<T, NumberComponents, false>: public DynamicArray<T> {
  protected:
    using field = DynamicArray<T>::dArray;
    Index index;

  public:
    const std::string fieldName;

    LocalizedField(const std::string& fieldName_in,
          const MathVector<unsigned int, 3>& lengths_in)
      : fieldName(fieldName_in)
      , field(NULL)
    {}

    LocalizedField(const std::string& fieldName_in,
                   const DynamicArray<T>& field_in)
      : fieldName(fieldName_in)
      , field(NULL)
    {}

    T* getRawField(const unsigned int iC) {
      return NULL;
    }

    T getField(const int idx, const unsigned int iC) {
      return -1;
    }

    void setField(const int idx, const T value, const unsigned int iC) {
    }

  };

template <class T, unsigned int NumberComponents, bool IsWritten>
  class Field {};

  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, true>: public DynamicArray<T> {
  protected:
    using field = DynamicArray<T>::dArray;
    static constexpr unsigned int numberComponents = NumberComponents;
    static constexpr bool isWritten = IsWritten;

  public:
    const std::string fieldName;

    LocalizedField(const std::string& fieldName_in)
      : fieldName(fieldName_in)
      , field(NumberComponents_in*lengths[d::X]*lengths[d::Y]*lengths[d::Z])
      {}

    LocalizedField(const std::string& fieldName_in,
                   const DynamicArray<T>& field_in)
      : fieldName(fieldName_in)
      , field(field_in)
    {}

    T* getRawField(const unsigned int iC) {
      return field[domainVolume * iC];
    }

    T getField(const int idx, const unsigned int iC) {
      return field[iC * domainVolume + idx];
    } // operator overloading????

    MathVector<T, NumberComponents> getVector(const int idx) {
      MathVector<T, NumberComponents> vectorR;
      UnrolledFor<0, NumberComponents>::Do([&] (int iC) {
          vectorR[iC] = getField(idx, iC);
      });

      return vectorR;
    }

    void setField(const int idx, const T value, const unsigned int iC) {
      field[iC * domainVolume + idx] = value;
    }

    void setField(const int idx, const MathVector<T, NumberComponents> value) {
      UnrolledFor<0, NumberComponents>::Do([&] (int iC) {
          setField(idx, value[iC], iC);
      });
    }
  };


}
