#ifndef FIELD_H
#define FIELD_H

#include <string>
#include <fstream>
#include <ostream>

#include "Commons.h"
#include "Options.h"
#include "Domain.h"
#include "DynamicArray.cuh"
#include "MathVector.h"

namespace lbm {

   /**
   * Field containing data sets to be stored, communicated, and dumped.
   *
   * @tparam T datatype.
   * @tparam unsigned int NumberComponents of the field.
   * @tparam Architecture on which the code is executed.
   * @tparam bool IsWritten whether the current field should be dumped
   * and therefore allocated.
   */

  template <class T, unsigned int NumberComponents, DomainType initDomainType,
            Architecture architecture, bool IsWritten>
  class Field {};


  template <class T, unsigned int NumberComponents, DomainType initDomainType,
            Architecture architecture>
  class Field<T, NumberComponents, initDomainType, architecture, false>  {
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
          const DynamicArray<T, architecture>& globalArray_in)
      : fieldName(fieldName_in)
    {}

    T * RESTRICT getGlobalData() {
      return NULL;
    }

    void setGlobalValue(const unsigned int index, const T value) {}

    void setGlobalValue(const MathVector<unsigned int, 3>& iP,
                        const T value,
                        const unsigned int iC) {}

    void setGlobalVector(const MathVector<unsigned int, 3>& iP,
                         const MathVector<T, NumberComponents> vector) {}

    DynamicArray<T, Architecture::CPU> getGlobalArray() {
      return DynamicArray<T, Architecture::CPU>();
    }

    T getGlobalValue(const MathVector<unsigned int, 3>& iP,
                     const unsigned int iC = 0) const {
      return (T) -1;
    }

    MathVector<T, NumberComponents> getGlobalVector(const MathVector<unsigned int, 3>& iP) const {
      return MathVector<T, NumberComponents>{{(T) -1}};
    }

    T * RESTRICT getLocalData() {
      return NULL;
    }

    DynamicArray<T, Architecture::CPU> getLocalArray() {
      return DynamicArray<T, Architecture::CPU>();
    }

    DEVICE HOST
    void setLocalValue(const unsigned int index, const T value, const unsigned int iC = 0) {}

    DEVICE HOST
    void setLocalVector(const unsigned int index,
                       const MathVector<T, NumberComponents> vector) {}

    DEVICE HOST
    T getLocalValue(const unsigned int index, const unsigned int iC = 0) {
      return (T) -1;
    }

    DEVICE HOST
    MathVector<T, NumberComponents> getLocalVector(const unsigned int index) {
      return MathVector<T, NumberComponents>{{(T) -1}};
    }

  };

  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, DomainType::Generic, Architecture::Generic, true> {
  public:
    static constexpr unsigned int numberComponents = NumberComponents;
    static constexpr bool IsWritten = true;

    const std::string fieldName;

    Field(const std::string& fieldName_in)
      : fieldName(fieldName_in)
    {}

  };




  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, DomainType::Generic,
              Architecture::CPU, true>
    : public Field<T, NumberComponents, DomainType::Generic,
                 Architecture::Generic, true> {
  private:
    using Base = Field<T, NumberComponents, DomainType::Generic,
                       Architecture::Generic, true>;

  protected:
    DynamicArray<T, Architecture::CPU> localArray;

  public:
    using Base::numberComponents;
    using Base::IsWritten;
    using Base::fieldName;

    Field(const std::string& fieldName_in)
      : Base(fieldName_in)
      , localArray(NumberComponents*lSD::volume())
    {}

    Field(const std::string& fieldName_in,
          const T& value_in)
      : Base(fieldName_in)
      , localArray(NumberComponents*lSD::volume(), value_in)
    {}

    Field(const std::string& fieldName_in,
          const MathVector<T, NumberComponents>& vector_in)
      : Base(fieldName_in)
      , localArray(NumberComponents*lSD::volume(), vector_in[0])
    {
      //computationLocal = Computation
    }

    Field(const std::string& fieldName_in,
          const DynamicArray<T, Architecture::CPU>& localArray_in)
      : Base(fieldName_in)
      , localArray(localArray_in)
    {}

  };


  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, DomainType::Generic,
              Architecture::GPU, true>
    : public Field<T, NumberComponents, DomainType::Generic,
                 Architecture::Generic, true> {
  private:
    using Base = Field<T, NumberComponents, DomainType::Generic,
                       Architecture::Generic, true>;
  protected:
    DynamicArray<T, Architecture::CPUPinned> localArray;

  public:
    using Base::numberComponents;
    using Base::IsWritten;

    using Base::fieldName;

    Field(const std::string& fieldName_in)
      : Base(fieldName_in)
      , localArray(NumberComponents*lSD::volume())
    {}

    Field(const std::string& fieldName_in,
          const T& value_in)
      : Base(fieldName_in)
      , localArray(NumberComponents*lSD::volume(), value_in)
    {}

    Field(const std::string& fieldName_in,
          const MathVector<T, NumberComponents>& vector_in)
      : Base(fieldName_in)
      , localArray(NumberComponents*lSD::volume(), vector_in[0])
    {}

    Field(const std::string& fieldName_in,
          const DynamicArray<T, Architecture::CPU>& localArray_in)
      : Base(fieldName_in)
      , localArray(localArray_in)
    {}

  };


  template <class T, unsigned int NumberComponents, Architecture architecture>
  class Field<T, NumberComponents, DomainType::LocalSpace,
              architecture, true>
    : public Field<T, NumberComponents, DomainType::Generic,
                   architecture, true> {
  private:
    using Base = Field<T, NumberComponents, DomainType::Generic,
                       architecture, true>;
  protected:
    using Base::localArray;

  public:
    using Base::numberComponents;
    using Base::IsWritten;
    using Base::fieldName;

    using Base::Field;

    DynamicArray<T, Architecture::CPU>& getLocalArray() {
      return localArray;
    }

    DEVICE HOST
    T * RESTRICT getLocalData() {
      return localArray.data();
    }

    DEVICE HOST
    void setLocalValue(const unsigned int index,
                       const T value,
                       const unsigned int iC = 0) {
      localArray[lSD::getIndex(index, iC)] = value;
    }

    DEVICE HOST
    void setLocalVector(const unsigned int index,
                        const MathVector<T, NumberComponents> vector) {
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        setLocalValue(index, vector[iC], iC);
      }
    }

    DEVICE HOST
      T getLocalValue(const MathVector<unsigned int, 3>& iP,
                          const unsigned int iC = 0) const {
      return localArray[lSD::getIndex(iP, iC)];
    }

    DEVICE HOST
    MathVector<T, NumberComponents> getLocalVector(const MathVector<unsigned int, 3>& iP) const {
      MathVector<T, NumberComponents> vectorR;
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        vectorR[iC] = getLocalValue(iP, iC);
      }

      return vectorR;
    }

  };


  template <class T, unsigned int NumberComponents, Architecture architecture>
  class Field<T, NumberComponents, DomainType::GlobalSpace,
              architecture, true>
    : public Field<T, NumberComponents, DomainType::LocalSpace,
                   architecture, true> {
  private:
    using Base = Field<T, NumberComponents, DomainType::LocalSpace,
                       architecture, true>;
    typedef Domain<DomainType::GlobalSpace, partitionningT,
                   MemoryLayout::Generic, NumberComponents> gNCD;

  protected:
    DynamicArray<T, Architecture::CPU> globalArray;

  public:
    using Base::numberComponents;
    using Base::IsWritten;
    using Base::fieldName;

    Field(const std::string& fieldName_in)
      : Base(fieldName_in)
      , globalArray(gSD::volume()*NumberComponents)
    {}

    Field(const std::string& fieldName_in,
          const MathVector<T, NumberComponents>& vector_in)
      : Base(fieldName_in)
      , globalArray(gSD::volume()*NumberComponents, vector_in[0])
    {}

    Field(const std::string& fieldName_in,
          const T& value_in)
      : Base(fieldName_in)
      , globalArray(gSD::volume()*NumberComponents, value_in)
    {}

    Field(const std::string& fieldName_in,
          const DynamicArray<T, Architecture::CPU>& globalArray_in)
      : Base(fieldName_in)
      , globalArray(globalArray_in)
    {}

    using Base::getLocalArray;
    using Base::getLocalData;
    using Base::setLocalValue;
    using Base::setLocalVector;
    using Base::getLocalValue;
    using Base::getLocalVector;

    T * RESTRICT getGlobalData() {
      return globalArray.data();
    }

    DynamicArray<T, Architecture::CPU>& getGlobalArray() {
      return globalArray;
    }

    void setGlobalField(DynamicArray<T, Architecture::CPU> globalArray_in) {
      globalArray.copyFrom(globalArray_in);
    }

    void setGlobalValue(const unsigned int index, const T value) {
      globalArray[index] = value;
    }

    void setGlobalValue(const MathVector<unsigned int, 3>& iP,
                        const T value,
                        const unsigned int iC) {
      setGlobalValue(gNCD::getIndex(iP, iC), value);
    }

    void setGlobalVector(const MathVector<unsigned int, 3>& iP,
                         const MathVector<T, NumberComponents> vector) {
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        setGlobalValue(iP, vector[iC], iC);
      }
    }

    T getGlobalValue(const MathVector<unsigned int, 3>& iP,
                     const unsigned int iC = 0) const {
      return globalArray[gNCD::getIndex(iP, iC)];
    }

    MathVector<T, NumberComponents> getGlobalVector(const MathVector<unsigned int, 3>& iP) const {
      MathVector<T, NumberComponents> vectorR;
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        vectorR[iC] = getGlobalValue(iP, iC);
      }

      return vectorR;
    }

  };

} // namespace lbm

#endif // FIELD_H
