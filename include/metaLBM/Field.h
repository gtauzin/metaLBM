#pragma once

#include <fstream>
#include <ostream>
#include <string>

#include "Commons.h"
#include "Computation.h"
#include "Domain.h"
#include "DynamicArray.cuh"
#include "Helpers.h"
#include "MathVector.h"
#include "Options.h"

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

template <class T,
          unsigned int NumberComponents,
          Architecture architecture,
          bool IsWritten>
class FieldAllocator {};

template <class T, unsigned int NumberComponents>
class FieldAllocator<T, NumberComponents, Architecture::Generic, true> {
 public:
  const unsigned int numberElements;
  static constexpr bool IsWritten = true;

  const std::string fieldName;

  FieldAllocator(const std::string& fieldName_in,
                 const unsigned int numberElements_in)
      : fieldName(fieldName_in),
        numberElements(numberElements_in)

  {}
};

template <class T, unsigned int NumberComponents>
class FieldAllocator<T, NumberComponents, Architecture::CPU, true>
    : public FieldAllocator<T, NumberComponents, Architecture::Generic, true> {
 private:
  using Base = FieldAllocator<T, NumberComponents, Architecture::Generic, true>;

 protected:
  DynamicArray<T, Architecture::CPU> localArray;

 public:
  using Base::fieldName;
  using Base::IsWritten;
  using Base::numberElements;

  FieldAllocator(const std::string& fieldName_in,
                 const unsigned int numberElements_in)
      : Base(fieldName_in, numberElements_in),
        localArray(numberElements_in * NumberComponents) {}
};

template <class T, unsigned int NumberComponents>
class FieldAllocator<T, NumberComponents, Architecture::GPU, true>
    : public FieldAllocator<T, NumberComponents, Architecture::Generic, true> {
 private:
  using Base = FieldAllocator<T, NumberComponents, Architecture::Generic, true>;

 protected:
  DynamicArray<T, Architecture::CPUPinned> localArray;

 public:
  using Base::IsWritten;
  using Base::numberElements;

  using Base::fieldName;

  FieldAllocator(const std::string& fieldName_in,
                 const unsigned int numberElements_in)
      : Base(fieldName_in, numberElements_in),
        localArray(numberElements_in * NumberComponents) {}
};

template <class T,
          unsigned int NumberComponents,
          Architecture architecture,
          bool isWritten>
class Field {};

template <class T, unsigned int NumberComponents, Architecture architecture>
class Field<T, NumberComponents, architecture, true>
    : public FieldAllocator<T, NumberComponents, architecture, true> {
 private:
  using Base = FieldAllocator<T, NumberComponents, architecture, true>;

 protected:
  using Base::localArray;

 public:
  using Base::fieldName;
  using Base::IsWritten;
  using Base::numberElements;

  Field(const std::string& fieldName_in, const unsigned int numberElements_in)
      : Base(fieldName_in, numberElements_in) {}

  Field(const std::string& fieldName_in,
        const unsigned int numberElements_in,
        const T& value_in)
      : Base(fieldName_in, numberElements_in) {
    Computation<architecture, L::dimD> computationLocal(lSD::sStart(),
                                                        lSD::sEnd());
    computationLocal.Do([&] LBM_HOST(const Position& iP) {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        setLocalValue(iP, value_in, iC);
      }
    });
    computationLocal.synchronize();

  }

  Field(const std::string& fieldName_in,
        const unsigned int numberElements_in,
        const MathVector<T, NumberComponents>& vector_in)
      : Base(fieldName_in, numberElements_in) {
    Computation<architecture, L::dimD> computationLocal(lSD::sStart(),
                                                             lSD::sEnd());
    computationLocal.Do(
        [&] LBM_HOST(const Position& iP) { setLocalVector(iP, vector_in); });
    computationLocal.synchronize();

  }

  Field(const std::string& fieldName_in,
        const DynamicArray<T, Architecture::CPU>& localArray_in)
      : Base(fieldName_in, localArray_in.size()) {
    localArray.copyFrom(localArray_in);
  }

  inline DynamicArray<T, Architecture::CPU>& getLocalArray() {
    return localArray;
  }

  LBM_DEVICE LBM_HOST T* getLocalData(const unsigned int iC = 0) {
    return localArray.data(iC * numberElements);
  }

  LBM_DEVICE LBM_HOST inline void setLocalValue(const Position iP,
                                                const T value,
                                                const unsigned int iC = 0) {
    (localArray.data(iC * numberElements))[lSD::getIndex(iP)] = value;
  }

  LBM_DEVICE LBM_HOST inline void setLocalVector(
      const Position iP,
      const MathVector<T, NumberComponents> vector) {
    for (auto iC = 0; iC < NumberComponents; ++iC) {
      setLocalValue(iP, vector[iC], iC);
    }
  }

  LBM_DEVICE LBM_HOST inline T getLocalValue(const Position& iP,
                                             const unsigned int iC = 0) const {
    return (localArray.data(iC * numberElements))[lSD::getIndex(iP)];
  }

  LBM_DEVICE LBM_HOST inline MathVector<T, NumberComponents> getLocalVector(
      const Position& iP) const {
    MathVector<T, NumberComponents> vectorR;
    for (auto iC = 0; iC < NumberComponents; ++iC) {
      vectorR[iC] = getLocalValue(iP, iC);
    }

    return vectorR;
  }
};

template <class T, unsigned int NumberComponents, Architecture architecture>
class Field<T, NumberComponents, architecture, false> {
 public:
  const unsigned int numberElements;
  static constexpr bool IsWritten = false;

  const std::string fieldName;

  Field(const std::string& fieldName_in,
        const unsigned int numberElements_in,
        const MathVector<T, NumberComponents>& vector_in)
      : fieldName(fieldName_in), numberElements(numberElements_in) {}

  Field(const std::string& fieldName_in,
        const unsigned int numberElements_in,
        const T& value_in = (T)0)
      : fieldName(fieldName_in), numberElements(numberElements_in) {}

  Field(const std::string& fieldName_in,
        const unsigned int numberElements_in,
        const DynamicArray<T, architecture>& globalArray_in)
      : fieldName(fieldName_in), numberElements(numberElements_in) {}

  T* getLocalData(const unsigned int iC = 0) { return NULL; }

  DynamicArray<T, Architecture::CPU> getLocalArray() {
    return DynamicArray<T, Architecture::CPU>();
  }

  LBM_DEVICE LBM_HOST void setLocalValue(const unsigned int index,
                                         const T value,
                                         const unsigned int iC = 0) {}

  LBM_DEVICE LBM_HOST void setLocalVector(
      const unsigned int index,
      const MathVector<T, NumberComponents> vector) {}

  LBM_DEVICE LBM_HOST T getLocalValue(const unsigned int index,
                                      const unsigned int iC = 0) {
    return (T)-1;
  }

  LBM_DEVICE LBM_HOST MathVector<T, NumberComponents> getLocalVector(
      const unsigned int index) {
    return MathVector<T, NumberComponents>{{(T)-1}};
  }
};

}  // namespace lbm
