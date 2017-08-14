#ifndef FIELD_H
#define FIELD_H

#include <string>
#include <fstream>
#include <ostream>

#include "Options.h"
#include "DynamicArray.h"
#include "MathVector.h"
#include "Domain.h"

namespace lbm {

  /**
   * LocalizedField containing domain-defined data sets.
   *
   * @tparam T datatype.
   * @tparam unsigned int NumberComponents of the field.
   */

  template <class T, unsigned int NumberComponents>
  class LocalizedField {
  protected:
    DynamicArray<T> field;
    unsigned int volume;

  public:
    static constexpr unsigned int numberComponents = NumberComponents;
    const std::string fieldName;

    LocalizedField(const std::string& fieldName_in)
      : field(0)
      , volume(0)
      , fieldName(fieldName_in)
    {}

    LocalizedField(const std::string& fieldName_in,
                   const unsigned int volume_in,
                   const MathVector<T, NumberComponents>& value_in)
      : field(NumberComponents*volume_in)
      , volume(volume_in)
      , fieldName(fieldName_in)
    {
      for(unsigned int i = 0; i < volume_in; ++i) setField(i, value_in);
    }

    LocalizedField(const std::string& fieldName_in,
                   const unsigned int volume_in,
                   const T& value_in = (T) 0)
      : field(NumberComponents*volume_in, value_in)
      , volume(volume_in)
      , fieldName(fieldName_in)
    {}

    T * __restrict__ data(const unsigned int iC = 0) {
      return field.data() + volume * iC;
    }

    unsigned int getVolume() {
      return volume;
    }

    DynamicArray<T>& getArray() {
      return field;
    }

    void swap(LocalizedField<T, NumberComponents>& field_in) {
      field.swap(field_in.getArray());
    }

    T getField(const unsigned int index, const unsigned int iC = 0) const {
      return field[iC * volume + index];
    }

    MathVector<T, NumberComponents> getVector(const unsigned int index) const {
      MathVector<T, NumberComponents> vectorR;
      UnrolledFor<0, NumberComponents>::Do([&] (unsigned int iC) {
          vectorR[iC] = getField(index, iC);
      });

      return vectorR;
    }

    void setField(const unsigned int index, const T value, const unsigned int iC = 0) {
      field[iC * volume + index] = value;
    }

    void setField(LocalizedField<T, NumberComponents> field_in) {
      field.copy(field_in.getArray());
      volume = field_in.getVolume();
    }

    void setField(const unsigned int index, const MathVector<T, NumberComponents>& value) {
      UnrolledFor<0, NumberComponents>::Do([&] (unsigned int iC) {
          setField(index, value[iC], iC);
      });
    }

    T& operator[] (int i) {
      return field[i];
    }

    const T& operator[] (int i) const {
      return field[i];
    }

  };


    /**
   * Field containing data sets to be communicated and dumped.
   *
   * @tparam T datatype.
   * @tparam unsigned int NumberComponents of the field.
   * @tparam bool IsWritten whether the current field should be dumped
   * and therefore allocated.
   */

  template <class T, unsigned int NumberComponents, bool IsWritten>
  class Field {};

  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, true> {
  protected:
    LocalizedField<T, NumberComponents> globalField;
    LocalizedField<T, NumberComponents> localField;

  public:
    static constexpr unsigned int numberComponents = NumberComponents;
    static constexpr bool IsWritten = true;

    const std::string fieldName;

    Field(const std::string& fieldName_in)
      : globalField(fieldName_in, gD::volume())
      , localField(fieldName_in, lD::volume())
      , fieldName(fieldName_in)
    {}


    T * __restrict__ localData(const unsigned int iC = 0) {
      return localField.data(iC);
    }

    T * __restrict__ globalData(const unsigned int iC = 0) {
      return globalField.data(iC);
    }

    LocalizedField<T, NumberComponents> getGlobalField() {
      return globalField;
    }

    T getLocalField(const unsigned int index, const unsigned int iC = 0) {
      return localField.getField(index, iC);
    }

    MathVector<T, NumberComponents> getLocalVector(const unsigned int index) {
      return localField.getVector(index);
    }

    void setGlobalField(LocalizedField<T, NumberComponents> globalField_in) {
      globalField.setField(globalField_in);
    }

    void setLocalField(const unsigned int index, const T value) {
      localField.setField(index, value);
    }


    void setLocalField(const unsigned int index, const MathVector<T, NumberComponents> value) {
      localField.setField(index, value);
    }

    T getGlobalField(const unsigned int index, const unsigned int iC = 0) {
      return globalField.getField(index, iC);
    }

    MathVector<T, NumberComponents> getGlobalVector(const unsigned int index) {
      return globalField.getVector(index);
    }

    void setGlobalField(const unsigned int index, const MathVector<T, NumberComponents> value) {
      globalField.setField(index, value);
    }

  };


  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, false>  {
  public:
    static constexpr unsigned int numberComponents = NumberComponents;
    static constexpr bool IsWritten = false;

    const std::string fieldName;

    Field(const std::string& fieldName_in)
    {}

    T * __restrict__ localData(const unsigned int iC = 0) {
      return NULL;
    }

    T * __restrict__ globalData(const unsigned int iC = 0) {
      return NULL;
    }


    T getLocalField(const unsigned int index, const unsigned int iC = 0) {
      return (T) -1;
    }

    MathVector<T, NumberComponents> getLocalVector(const unsigned int index) {
      return MathVector<T, NumberComponents>{{(T) -1}};
    }

    void setLocalField(const unsigned int index, const MathVector<T, NumberComponents> value) {
    }

    T getGlobalField(const unsigned int index, const unsigned int iC = 0) {
      return (T) -1;
    }

    MathVector<T, NumberComponents> getGlobalVector(const unsigned int index) {
      return MathVector<T, NumberComponents>{{(T) -1}};
    }

    void setGlobalField(const unsigned int index,
                        const MathVector<T, NumberComponents> value) {
    }

  };

  template<class T>
  class Fields {};

}

#endif // FIELD_H
