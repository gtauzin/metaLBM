#ifndef FIELD_H
#define FIELD_H

#include <string>
#include <fstream>
#include <ostream>

#include "Commons.h"
#include "Options.h"
#include "Domain.h"
#include "DynamicArray.h"
#include "MathVector.h"

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
                   const T& value_in = (T) 0)
      : field(NumberComponents*volume_in, value_in)
      , volume(volume_in)
      , fieldName(fieldName_in)
    {}

    LocalizedField(const std::string& fieldName_in,
                   const LocalizedField<T, NumberComponents>& localizedField_in)
      : field(localizedField_in.getArray())
      , volume(localizedField_in.getVolume())
      , fieldName(fieldName_in)
    {}


    T * RESTRICT data(const unsigned int iC = 0) {
      return field.data() + volume * iC;
    }

    unsigned int getVolume() const {
      return volume;
    }

    const DynamicArray<T>& getArray() const {
      return field;
    }

    DynamicArray<T>& getArray() {
      return field;
    }

    void swap(LocalizedField<T, NumberComponents>& field_in) {
      field.swap(field_in.getArray());
    }

    void setField(LocalizedField<T, NumberComponents> field_in) {
      field.copy(field_in.getArray());
      volume = field_in.getVolume();
    }

    T& operator[] (int i) {
      return field[i];
    }

    const T& operator[] (int i) const {
      return field[i];
    }

  };


    /**
   * Field containing data sets to be stored, communicated, and dumped.
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
  private:
    typedef Domain<DomainType::Global, partitionningT,
                   MemoryLayout::Generic, NumberComponents> gNCD;
    typedef Domain<DomainType::Local, PartitionningType::Generic,
                   MemoryLayout::Generic, NumberComponents> lNCD;

  protected:
    LocalizedField<T, NumberComponents> globalField;
    LocalizedField<T, NumberComponents> localField;

  public:
    static constexpr unsigned int numberComponents = NumberComponents;
    static constexpr bool IsWritten = true;

    const std::string fieldName;

    Field(const std::string& fieldName_in,
          const MathVector<T, NumberComponents>& vector_in)
      : globalField(fieldName_in, gD::volume())
      , localField(fieldName_in, lD::volume())
      , fieldName(fieldName_in)
    {
      MathVector<unsigned int, 3> iP;
      for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++) {
        for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++) {
          for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++) {
            iP = {iX, iY, iZ};
            setGlobalVector(iP, vector_in);
          }
        }
      }
    }

    Field(const std::string& fieldName_in,
          const T& value_in = (T) 0)
      : globalField(fieldName_in, gD::volume())
      , localField(fieldName_in, lD::volume())
      , fieldName(fieldName_in)
    {
      for(unsigned int i = 0; i < gD::volume(); ++i) {
        setGlobalValue(i, value_in);
      }
    }

    Field(const std::string& fieldName_in,
          const LocalizedField<T, NumberComponents>& globalField_in)
      : globalField(fieldName_in, globalField_in)
      , localField(fieldName_in, lD::volume())
      , fieldName(fieldName_in)
    {}

    T * RESTRICT localData() {
      return localField.data();
    }

    void setLocalValue(const unsigned int index, const T value, unsigned int iC = 0) {
      localField[lNCD::getIndex(index, iC)] = value;
    }

    void setLocalVector(const unsigned int index,
                       const MathVector<T, NumberComponents> vector) {
      UnrolledFor<0, NumberComponents>::Do([&] (unsigned int iC) {
          setLocalValue(index, vector[iC], iC);
      });
    }

    T getLocalValue(const unsigned int index, const unsigned int iC = 0) {
      return localField[lD::getIndex(index, iC)];
    }

    MathVector<T, NumberComponents> getLocalVector(const unsigned int index) {
      MathVector<T, NumberComponents> vectorR;
      UnrolledFor<0, NumberComponents>::Do([&] (unsigned int iC) {
          vectorR[iC] = getLocalValue(index, iC);
      });

      return vectorR;
    }

    T * RESTRICT globalData() {
      return globalField.data();
    }

    LocalizedField<T, NumberComponents> getGlobalField() {
      return globalField;
    }

    void setGlobalField(LocalizedField<T, NumberComponents> globalField_in) {
      globalField.setField(globalField_in);
    }

    void setGlobalValue(const unsigned int index, const T value) {
      globalField[index] = value;
    }

    void setGlobalValue(const MathVector<unsigned int, 3>& iP,
                        const T value,
                        const unsigned int iC) {
      globalField[gNCD::getIndex(iP, iC)] = value;
    }

    void setGlobalVector(const MathVector<unsigned int, 3>& iP,
                         const MathVector<T, NumberComponents> vector) {
      UnrolledFor<0, NumberComponents>::Do([&] (unsigned int iC) {
          setGlobalValue(iP, vector[iC], iC);
      });
    }

    T getGlobalValue(const MathVector<unsigned int, 3>& iP, const unsigned int iC = 0) const {
      return globalField[gNCD::getIndex(iP, iC)];
    }

    MathVector<T, NumberComponents> getGlobalVector(const MathVector<unsigned int, 3>& iP) const {
      MathVector<T, NumberComponents> vectorR;
      UnrolledFor<0, NumberComponents>::Do([&] (unsigned int iC) {
          vectorR[iC] = getGlobalValue(iP, iC);
      });

      return vectorR;
    }

  };

  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, false>  {
  public:
    static constexpr unsigned int numberComponents = NumberComponents;
    static constexpr bool IsWritten = false;

    const std::string fieldName;

    Field(const std::string& fieldName_in,
          const MathVector<T, NumberComponents>& vector_in)
      : fieldName(fieldName_in)
    {}

    Field(const std::string& fieldName_in,
          const T& value_in = (T) 0)
      : fieldName(fieldName_in)
    {}

    Field(const std::string& fieldName_in,
          const LocalizedField<T, NumberComponents>& globalField_in)
      : fieldName(fieldName_in)
    {}

    T * RESTRICT localData() {
      return NULL;
    }

    void setLocalValue(const unsigned int index, const T value, unsigned int iC = 0) {
    }

    void setLocalVector(const unsigned int index,
                       const MathVector<T, NumberComponents> vector) {
    }

    T getLocalValue(const unsigned int index, const unsigned int iC = 0) {
      return (T) -1;
    }

    MathVector<T, NumberComponents> getLocalVector(const unsigned int index) {
      return MathVector<T, NumberComponents>{{(T) -1}};
    }

    T * RESTRICT globalData() {
      return NULL;
    }

    void setGlobalField(LocalizedField<T, NumberComponents> globalField_in) {
    }

    void setGlobalValue(const unsigned int index, const T value) {
    }

    void setGlobalValue(const MathVector<unsigned int, 3>& iP,
                        const T value,
                        const unsigned int iC) {
    }

    void setGlobalVector(const MathVector<unsigned int, 3>& iP,
                         const MathVector<T, NumberComponents> vector) {
    }

    LocalizedField<T, NumberComponents> getGlobalField() {
      return LocalizedField<T, NumberComponents>("empty");
    }

    T getGlobalValue(const MathVector<unsigned int, 3>& iP, const unsigned int iC = 0) const {
      return (T) -1;
    }

    MathVector<T, NumberComponents> getGlobalVector(const MathVector<unsigned int, 3>& iP) const {
      return MathVector<T, NumberComponents>{{(T) -1}};
    }
  };

}

#endif // FIELD_H
