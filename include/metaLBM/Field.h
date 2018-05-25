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

  template <class T, unsigned int NumberComponents, Architecture architecture,
          bool IsWritten>
  class FieldAllocator {};

  template <class T, unsigned int NumberComponents>
  class FieldAllocator<T, NumberComponents, Architecture::Generic, true> {
  public:
    static constexpr bool IsWritten = true;

    const std::string fieldName;

    FieldAllocator(const std::string& fieldName_in)
      : fieldName(fieldName_in)
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

    FieldAllocator(const std::string& fieldName_in)
      : Base(fieldName_in)
      , localArray(FFTWInit::numberElements * NumberComponents)
    {}
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
    using Base::fieldName;

    FieldAllocator(const std::string& fieldName_in)
      : Base(fieldName_in)
      , localArray(FFTWInit::numberElements * NumberComponents)
    {}
  };

  template <class T, unsigned int NumberComponents, Architecture architecture,
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
    Computation<architecture, L::dimD> computationLocal;

    Field(const std::string& fieldName_in)
      : Base(fieldName_in)
      , computationLocal(lSD::sStart(), lSD::sEnd())
    {}

    Field(const std::string& fieldName_in,
          const T& value_in, const Stream<architecture>& stream_in)
      : Base(fieldName_in)
      , computationLocal(lSD::sStart(), lSD::sEnd())
    {
      T * localArrayPtr = localArray.data();
      computationLocal.Do(stream_in, *this, localArrayPtr, value_in,
                          FFTWInit::numberElements);
      computationLocal.synchronize();
    }

    Field(const std::string& fieldName_in,
          const MathVector<T, NumberComponents>& vector_in,
          const Stream<architecture>& stream_in)
      : Base(fieldName_in)
      , computationLocal(lSD::sStart(), lSD::sEnd())
  {
    T * localArrayPtr = localArray.data();
    computationLocal.Do(stream_in, *this, localArrayPtr, vector_in,
                        FFTWInit::numberElements);
    computationLocal.synchronize();
  }

    Field(const std::string& fieldName_in,
          const DynamicArray<T, Architecture::CPU>& localArray_in)
      : Base(fieldName_in) {
      localArray.copyFrom(localArray_in);
    }

    LBM_HOST
    ~Field() {}

    LBM_DEVICE void operator()(const Position& iP, T * localArrayPtr,
                               const T initValue, const unsigned int numberElements) {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        localArrayPtr[iC * numberElements + lSD::getIndex(iP)]
        = initValue;
      }
    }

    LBM_DEVICE void operator()(const Position& iP, T * localArrayPtr,
                               const MathVector<T, L::dimD> initVector,
                               const unsigned int numberElements) {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        localArrayPtr[iC * numberElements + lSD::getIndex(iP)]
          = initVector[iC];
      }
    }


    void initByVector(const MathVector<T, NumberComponents>& vector_in,
                      const Stream<architecture>& stream_in) {
      Computation<architecture, L::dimD> computationLocal(lSD::sStart(),
                                                          lSD::sEnd());

      T * localArrayPtr = localArray.data();
      unsigned int numberElements = FFTWInit::numberElements;
      computationLocal.Do(stream_in, [=] LBM_HOST LBM_DEVICE (const Position& iP) {
          for (auto iC = 0; iC < NumberComponents; ++iC) {
            localArrayPtr[iC * numberElements + lSD::getIndex(iP)]
              = vector_in[iC];
          }
        });
      computationLocal.synchronize();
    }

    inline DynamicArray<T, Architecture::CPU>& getLocalArray() {
      return localArray;
    }

    LBM_DEVICE LBM_HOST T*
      getLocalData(const unsigned int numberElements, const unsigned int iC = 0) {
      return localArray.data(iC * numberElements);
    }

    LBM_DEVICE LBM_HOST inline
      void setLocalValue(const Position iP, const T value,
                         const unsigned int numberElements,
                         const unsigned int iC = 0) {
      (localArray.data(iC * numberElements))[lSD::getIndex(iP)]
        = value;
    }

    LBM_DEVICE LBM_HOST inline
    void setLocalVector(const Position iP,
                        const MathVector<T, NumberComponents> vector,
                        const unsigned int numberElements) {
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        setLocalValue(iP, vector[iC], numberElements, iC);
      }
    }

    LBM_DEVICE LBM_HOST inline
      T getLocalValue(const Position& iP,
                      const unsigned int numberElements,
                      const unsigned int iC = 0) const {
      return (localArray.data(iC * numberElements))[lSD::getIndex(iP)];
    }

    LBM_DEVICE LBM_HOST inline
    MathVector<T, NumberComponents> getLocalVector(const Position& iP,
                                                   const unsigned int numberElements) const {
      MathVector<T, NumberComponents> vectorR;
      for (auto iC = 0; iC < NumberComponents; ++iC) {
        vectorR[iC] = getLocalValue(iP, numberElements, iC);
      }

      return vectorR;
    }
  };

  template <class T, unsigned int NumberComponents, Architecture architecture>
  class Field<T, NumberComponents, architecture, false> {
  public:
    static constexpr bool IsWritten = false;

    const std::string fieldName;

    Field(const std::string& fieldName_in,
        const MathVector<T, NumberComponents>& vector_in)
      : fieldName(fieldName_in)
    {}

    Field(const std::string& fieldName_in,
          const T& value_in = (T)0)
      : fieldName(fieldName_in)
    {}

    Field(const std::string& fieldName_in,
          const DynamicArray<T, architecture>& globalArray_in)
      : fieldName(fieldName_in)
    {}

    T* getLocalData(const unsigned int numberElements,
                    const unsigned int iC = 0) { return NULL; }

    DynamicArray<T, Architecture::CPU> getLocalArray() {
      return DynamicArray<T, Architecture::CPU>();
    }

    LBM_DEVICE LBM_HOST T getLocalValue(const unsigned int index,
                                        const unsigned int numberElements,
                                        const unsigned int iC = 0) {
      return (T)-1;
    }

    LBM_DEVICE LBM_HOST
    MathVector<T, NumberComponents> getLocalVector(const unsigned int index,
                               const unsigned int numberElements) {
      return MathVector<T, NumberComponents>{{(T)-1}};
    }
  };

}  // namespace lbm
