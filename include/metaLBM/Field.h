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

  //TODO: I am pretty sure the

   /**
   * Field containing data sets to be stored, communicated, and dumped.
   *
   * @tparam T datatype.
   * @tparam unsigned int NumberComponents of the field.
   * @tparam Architecture on which the code is executed.
   * @tparam bool IsWritten whether the current field should be dumped
   * and therefore allocated.
   */

  template <class T, unsigned int NumberComponents,
            Architecture architecture, bool IsWritten>
  class Field {};


  template <class T, unsigned int NumberComponents, Architecture architecture>
  class Field<T, NumberComponents, architecture, false>  {
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
          const DynamicArray<T, architecture>& globalArrayHost_in)
      : fieldName(fieldName_in)
    {}

    T * RESTRICT globalData() {
      return NULL;
    }

    void setGlobalField(DynamicArray<T, Architecture::CPU> globalArrayHost_in) {
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

    DynamicArray<T, Architecture::CPU> globalArray() {
      return DynamicArray<T, Architecture::CPU>();
    }

    T getGlobalValue(const MathVector<unsigned int, 3>& iP,
                     const unsigned int iC = 0) const {
      return (T) -1;
    }

    MathVector<T, NumberComponents> getGlobalVector(const MathVector<unsigned int, 3>& iP) const {
      return MathVector<T, NumberComponents>{{(T) -1}};
    }

    DynamicArray<T, Architecture::CPU> localHostArray() {
      return DynamicArray<T, Architecture::CPU>();
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

    DynamicArray<T, Architecture::GPU> localDeviceArray() {
      return DynamicArray<T, Architecture::GPU>();
    }

  };


  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, Architecture::CPU, true> {
  private:
    typedef Domain<DomainType::Global, partitionningT,
                   MemoryLayout::Generic, NumberComponents> gNCD;

  protected:
    DynamicArray<T, Architecture::CPU> globalArrayHost;
    DynamicArray<T, Architecture::CPU> localArrayHost;
    DynamicArray<T, Architecture::GPU> localArrayDevice;

    typedef Domain<DomainType::Local, PartitionningType::Generic,
                   MemoryLayout::Generic, NumberComponents> lNCD;


  public:
    static constexpr unsigned int numberComponents = NumberComponents;
    static constexpr bool IsWritten = true;

    const std::string fieldName;

    Field(const std::string& fieldName_in,
          const MathVector<T, NumberComponents>& vector_in)
      : globalArrayHost(gD::volume()*NumberComponents)
      , localArrayHost(lD::volume()*NumberComponents)
      , localArrayDevice()
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
      : globalArrayHost(gD::volume()*NumberComponents)
      , localArrayHost(lD::volume()*NumberComponents)
      , localArrayDevice()
      , fieldName(fieldName_in)
    {
      for(unsigned int i = 0; i < gD::volume(); ++i) {
        setGlobalValue(i, value_in);
      }
    }

    Field(const std::string& fieldName_in,
          const DynamicArray<T, Architecture::CPU>& globalArrayHost_in)
      : globalArrayHost(globalArrayHost_in)
      , localArrayHost(lD::volume()*NumberComponents)
      , localArrayDevice()
      , fieldName(fieldName_in)
    {}

    T * RESTRICT globalData() {
      return globalArrayHost.data();
    }

    DynamicArray<T, Architecture::CPU>& globalArray() {
      return globalArrayHost;
    }

    void setGlobalField(DynamicArray<T, Architecture::CPU> globalArray_in) {
      globalArrayHost.copyFrom(globalArray_in);
    }

    void setGlobalValue(const unsigned int index, const T value) {
      globalArrayHost[index] = value;
    }

    void setGlobalValue(const MathVector<unsigned int, 3>& iP,
                            const T value,
                            const unsigned int iC) {
      setGlobalValue(gNCD::getIndex(iP, iC), value);
    }

    void setGlobalVector(const MathVector<unsigned int, 3>& iP,
                         const MathVector<T, NumberComponents> vector) {
      UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (unsigned int iC) {
          setGlobalValue(iP, vector[iC], iC);
      });
    }

    T getGlobalValue(const MathVector<unsigned int, 3>& iP,
                         const unsigned int iC = 0) const {
      return globalArrayHost[gNCD::getIndex(iP, iC)];
    }

    MathVector<T, NumberComponents> getGlobalVector(const MathVector<unsigned int, 3>& iP) const {
      MathVector<T, NumberComponents> vectorR;
      UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (unsigned int iC) {
          vectorR[iC] = getGlobalValue(iP, iC);
      });

      return vectorR;
    }

    DynamicArray<T, Architecture::CPU>& localHostArray() {
      return localArrayHost;
    }

    void setLocalValue(const unsigned int index,
                           const T value, unsigned int iC = 0) {
      localArrayHost[lNCD::getIndex(index, iC)] = value;
    }

    void setLocalVector(const unsigned int index,
                            const MathVector<T, NumberComponents> vector) {
      UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (unsigned int iC) {
          setLocalValue(index, vector[iC], iC);
      });
    }

    DynamicArray<T, Architecture::GPU>& localDeviceArray() {
      return localArrayDevice;
    }
  };


  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, Architecture::GPU, true>
    : public Field<T, NumberComponents, Architecture::CPU, true> {
  private:
    typedef Domain<DomainType::Local, PartitionningType::Generic,
                   MemoryLayout::Generic, NumberComponents> lNCD;

  protected:
    DynamicArray<T, Architecture::GPU> localArrayDevice;

  public:
    static constexpr unsigned int numberComponents = NumberComponents;
    static constexpr bool IsWritten = true;

    const std::string fieldName;

    Field(const std::string& fieldName_in,
          const MathVector<T, NumberComponents>& vector_in)
      : Field<T, NumberComponents, Architecture::CPU, true>(fieldName_in, vector_in)
      , localArrayDevice(lD::volume()*NumberComponents)
    {}

    Field(const std::string& fieldName_in,
          const T& value_in = (T) 0)
      : Field<T, NumberComponents, Architecture::CPU, true>(fieldName_in, value_in)
      , localArrayDevice(lD::volume()*NumberComponents)
    {}

    Field(const std::string& fieldName_in,
          const DynamicArray<T, Architecture::CPU>& globalArrayHost_in)
      : Field<T, NumberComponents, Architecture::CPU, true>(fieldName_in,
                                                            globalArrayHost_in)
      , localArrayDevice(lD::volume()*NumberComponents)
    {}

    using Field<T, NumberComponents, Architecture::CPU, true>::globalData;
    using Field<T, NumberComponents, Architecture::CPU, true>::globalArray;
    using Field<T, NumberComponents, Architecture::CPU, true>::setGlobalField;
    using Field<T, NumberComponents, Architecture::CPU, true>::setGlobalValue;
    using Field<T, NumberComponents, Architecture::CPU, true>::setGlobalVector;
    using Field<T, NumberComponents, Architecture::CPU, true>::getGlobalValue;
    using Field<T, NumberComponents, Architecture::CPU, true>::getGlobalVector;

    using Field<T, NumberComponents, Architecture::CPU, true>::localHostArray;

    DynamicArray<T, Architecture::GPU>& localDeviceArray() {
      return localArrayDevice;
    }

    void setLocalValue(const unsigned int index,
                       const T value, unsigned int iC = 0) {
      localArrayDevice[lNCD::getIndex(index, iC)] = value;
    }

    void setLocalVector(const unsigned int index,
                        const MathVector<T, NumberComponents> vector) {
      UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (unsigned int iC) {
          setLocalValue(index, vector[iC], iC);
      });
    }

  };

}

#endif // FIELD_H
