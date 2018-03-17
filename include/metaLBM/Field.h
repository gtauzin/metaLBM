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
#include "Helpers.h"

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
    const unsigned int numberElements;
    static constexpr bool IsWritten = false;

    const std::string fieldName;

  Field(const std::string& fieldName_in, const unsigned int numberElements_in,
        const MathVector<T, NumberComponents>& vector_in)
      : fieldName(fieldName_in)
      , numberElements(numberElements_in)
    {}

    Field(const std::string& fieldName_in, const unsigned int numberElements_in,
          const T& value_in = (T) 0)
      : fieldName(fieldName_in)
      , numberElements(numberElements_in)
    {}

    Field(const std::string& fieldName_in, const unsigned int numberElements_in,
          const DynamicArray<T, architecture>& globalArray_in)
      : fieldName(fieldName_in)
      , numberElements(numberElements_in)
    {}

    T * getGlobalData(const unsigned int offset = 0) {
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

    T * getLocalData(const unsigned int offset = 0) {
      return NULL;
    }

    MultiDynamicArray<T, Architecture::CPU, NumberComponents> getLocalArray() {
      return MultiDynamicArray<T, Architecture::CPU, NumberComponents>();
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
    const unsigned int numberElements;
    static constexpr bool IsWritten = true;

    const std::string fieldName;

    Field(const std::string& fieldName_in, const unsigned int numberElements_in)
      : fieldName(fieldName_in)
      , numberElements(numberElements_in)

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
    MultiDynamicArray<T, Architecture::CPU, NumberComponents> localArray;

  public:
    using Base::numberElements;
    using Base::IsWritten;
    using Base::fieldName;

  Field(const std::string& fieldName_in, const unsigned int numberElements_in)
    : Base(fieldName_in, numberElements_in)
      , localArray(numberElements_in)
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
    MultiDynamicArray<T, Architecture::CPUPinned, NumberComponents> localArray;

  public:
    using Base::numberElements;
    using Base::IsWritten;

    using Base::fieldName;

    Field(const std::string& fieldName_in, const unsigned int numberElements_in)
      : Base(fieldName_in, numberElements_in)
      , localArray(numberElements_in)
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
    using Base::numberElements;
    using Base::IsWritten;
    using Base::fieldName;

    Field(const std::string& fieldName_in, const unsigned int numberElements_in)
    : Base(fieldName_in, numberElements_in)
    {}

  Field(const std::string& fieldName_in, const unsigned int numberElements_in,
        const T& value_in)
      : Base(fieldName_in, numberElements_in)
    {
      Computation<Architecture::CPU,
                  L::dimD> computationLocal(lSD::sStart(),
                                            lSD::sEnd());
      computationLocal.Do
        ([&] HOST (const MathVector<unsigned int, 3>& iP) {
          for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
            setLocalValue(iP, value_in, iC);
          }
        });
    }

  Field(const std::string& fieldName_in, const unsigned int numberElements_in,
        const MathVector<T, NumberComponents>& vector_in)
      : Base(fieldName_in, numberElements_in)
    {
      Computation<Architecture::CPU,
                  L::dimD> computationLocal(lSD::sStart(),
                                            lSD::sEnd());
      computationLocal.Do
        ([&] HOST (const MathVector<unsigned int, 3>& iP) {
          setLocalVector(iP, vector_in);
        });
    }

    Field(const std::string& fieldName_in,
          const MultiDynamicArray<T, Architecture::CPU, NumberComponents>& localArray_in)
      : Base(fieldName_in, localArray_in.size()/NumberComponents)
    {
      localArray.copyFrom(localArray_in);
    }

    inline MultiDynamicArray<T, Architecture::CPU, NumberComponents>& getLocalArray() {
      return localArray;
    }

    DEVICE HOST
    T ** getMultiData() {
      return localArray.multiData();
    }

    DEVICE HOST
    T * getLocalData(const unsigned int offset = 0) {
      return localArray.data(offset);
    }


    DEVICE HOST
    inline void setLocalValue(const MathVector<unsigned int, 3> iP,
                              const T value,
                              const unsigned int iC = 0) {
      localArray[iC][lSD::getIndex(iP)] = value;
    }

    DEVICE HOST
    inline void setLocalVector(const MathVector<unsigned int, 3> iP,
                               const MathVector<T, NumberComponents> vector) {
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        setLocalValue(iP, vector[iC], iC);
      }
    }

    DEVICE HOST
    inline T getLocalValue(const MathVector<unsigned int, 3>& iP,
                           const unsigned int iC = 0) const {
      return localArray[iC][lSD::getIndex(iP)];
    }

    DEVICE HOST
    inline MathVector<T, NumberComponents> getLocalVector(const MathVector<unsigned int, 3>& iP) const {
      MathVector<T, NumberComponents> vectorR;
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        vectorR[iC] = getLocalValue(iP, iC);
      }

      return vectorR;
    }

  };


} // namespace lbm

#endif // FIELD_H
