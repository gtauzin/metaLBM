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

    DEVICE HOST
    void setLocalValue(const unsigned int index, const T value, const unsigned int iC = 0) {
    }

    DEVICE HOST
    void setLocalVector(const unsigned int index,
                       const MathVector<T, NumberComponents> vector) {
    }

    DEVICE HOST
    T getLocalValue(const unsigned int index, const unsigned int iC = 0) {
      return (T) -1;
    }

    DEVICE HOST
    MathVector<T, NumberComponents> getLocalVector(const unsigned int index) {
      return MathVector<T, NumberComponents>{{(T) -1}};
    }

    DEVICE HOST
    DynamicArray<T, Architecture::GPU> localDeviceArray() {
      return DynamicArray<T, Architecture::GPU>();
    }

  };

  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, DomainType::Generic, Architecture::Generic, true> {
  protected:
    DynamicArray<T, Architecture::CPU> localArrayHost;
    DynamicArray<T, Architecture::GPU> localArrayDevice;

  public:
    static constexpr unsigned int numberComponents = NumberComponents;
    static constexpr bool IsWritten = true;

    const std::string fieldName;

    Field(const std::string& fieldName_in)
      : fieldName(fieldName_in)
      , localArrayDevice(lSD::volume()*NumberComponents)
      , localArrayHost(lSD::volume()*NumberComponents)
    {}

    DynamicArray<T, Architecture::CPU>& localHostArray() {
      return localArrayHost;
    }

    DynamicArray<T, Architecture::GPU>& localDeviceArray() {
      return localArrayDevice;
    }

    DEVICE HOST
    T * RESTRICT localDeviceData() {
      return localArrayDevice.data();
    }

    DEVICE HOST
    T * RESTRICT localHostData() {
      return localArrayHost.data();
    }

  };


  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, DomainType::GlobalSpace, Architecture::Generic, true>
    : public Field<T, NumberComponents, DomainType::Generic, Architecture::Generic, true> {
  public:
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::numberComponents;
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::IsWritten;

    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::fieldName;

    Field(const std::string& fieldName_in,
          const MathVector<T, NumberComponents>& vector_in)
      : Field<T, NumberComponents, DomainType::Generic,
              Architecture::Generic, true>(fieldName_in)
      , globalArrayHost(gSD::volume()*NumberComponents)
    {
      MathVector<unsigned int, 3> iP;
      for(unsigned int iZ = gSD::start()[d::Z]; iZ < gSD::end()[d::Z]; iZ++) {
        for(unsigned int iY = gSD::start()[d::Y]; iY < gSD::end()[d::Y]; iY++) {
          for(unsigned int iX = gSD::start()[d::X]; iX < gSD::end()[d::X]; iX++) {
            iP =  MathVector<unsigned int, 3>({iX, iY, iZ});
            setGlobalVector(iP, vector_in);
          }
        }
      }
    }

    Field(const std::string& fieldName_in,
          const T& value_in = (T) 0)
      : Field<T, NumberComponents, DomainType::Generic,
              Architecture::Generic, true>(fieldName_in)
      , globalArrayHost(gSD::volume()*NumberComponents)
    {
      for(unsigned int i = 0; i < gSD::volume(); ++i) {
        setGlobalValue(i, value_in);
      }
    }

    Field(const std::string& fieldName_in,
          const DynamicArray<T, Architecture::CPU>& globalArrayHost_in)
      : Field<T, NumberComponents, DomainType::Generic,
              Architecture::Generic, true>(fieldName_in)
      , globalArrayHost(gSD::volume()*NumberComponents)
      {
      globalArrayHost.copyFrom(globalArrayHost_in);
    }

  protected:
    DynamicArray<T, Architecture::CPU> globalArrayHost;

    typedef Domain<DomainType::GlobalSpace, partitionningT,
                   MemoryLayout::Generic, NumberComponents> gNCD;

    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localArrayHost;
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localArrayDevice;

  public:
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localHostArray;
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localDeviceArray;

    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localHostData;
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localDeviceData;

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
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        setGlobalValue(iP, vector[iC], iC);
      }
    }

    T getGlobalValue(const MathVector<unsigned int, 3>& iP,
                     const unsigned int iC = 0) const {
      return globalArrayHost[gNCD::getIndex(iP, iC)];
    }

    MathVector<T, NumberComponents> getGlobalVector(const MathVector<unsigned int, 3>& iP) const {
      MathVector<T, NumberComponents> vectorR;
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        vectorR[iC] = getGlobalValue(iP, iC);
      }

      return vectorR;
    }
  };


  template <class T, unsigned int NumberComponents>
  class Field<T, NumberComponents, DomainType::LocalSpace, Architecture::Generic, true>
    : public Field<T, NumberComponents, DomainType::Generic, Architecture::Generic, true> {
  public:
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::numberComponents;
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::IsWritten;

    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::fieldName;

    Field(const std::string& fieldName_in,
          const T& value_in = (T) 0)
      : Field<T, NumberComponents, DomainType::Generic,
              Architecture::Generic, true>(fieldName_in)
      , globalArrayHost()
    {
      MathVector<unsigned int, 3> iP;
      for(unsigned int iZ = lSD::start()[d::Z]; iZ < lSD::end()[d::Z]; ++iZ) {
        for(unsigned int iY = lSD::start()[d::Y]; iY < lSD::end()[d::Y]; ++iY) {
          for(unsigned int iX = lSD::start()[d::X]; iX < lSD::end()[d::X]; ++iX) {
            iP =  MathVector<unsigned int, 3>({iX, iY, iZ});
            for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
              setLocalHostValue(iP, value_in, iC);
            }
          }
        }
      }
    }

    Field(const std::string& fieldName_in,
          const MathVector<T, NumberComponents>& vector_in)
      : Field<T, NumberComponents, DomainType::Generic,
              Architecture::Generic, true>(fieldName_in)
      , globalArrayHost()
    {
      MathVector<unsigned int, 3> iP;
      for(unsigned int iZ = lSD::start()[d::Z]; iZ < lSD::end()[d::Z]; ++iZ) {
        for(unsigned int iY = lSD::start()[d::Y]; iY < lSD::end()[d::Y]; ++iY) {
          for(unsigned int iX = lSD::start()[d::X]; iX < lSD::end()[d::X]; ++iX) {
            iP =  MathVector<unsigned int, 3>({iX, iY, iZ});
            setLocalHostVector(iP, vector_in);
          }
        }
      }
    }

    Field(const std::string& fieldName_in,
          const DynamicArray<T, Architecture::CPU>& localArrayHost_in)
      : Field<T, NumberComponents, DomainType::Generic,
              Architecture::Generic, true>(fieldName_in)
      , globalArrayHost()
    {
      localArrayHost.copyFrom(localArrayHost_in);
    }

  protected:
    DynamicArray<T, Architecture::CPU> globalArrayHost;

    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localArrayDevice;
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localArrayHost;

  public:
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localHostArray;
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localDeviceArray;

    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localHostData;
    using Field<T, NumberComponents, DomainType::Generic,
                Architecture::Generic, true>::localDeviceData;

    DEVICE HOST
      T getLocalHostValue(const MathVector<unsigned int, 3>& iP,
                          const unsigned int iC = 0) const {
      return localArrayHost[lSD::getIndex(iP, iC)];
    }

    DEVICE HOST
    MathVector<T, NumberComponents> getLocalHostVector(const MathVector<unsigned int, 3>& iP) const {
      MathVector<T, NumberComponents> vectorR;
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        vectorR[iC] = getLocalHostValue(iP, iC);
      }

      return vectorR;
    }

    T * RESTRICT globalData() {
      return globalArrayHost.data();
    }

    DynamicArray<T, Architecture::CPU>& globalArray() {
      return globalArrayHost;
    }

    void setGlobalField(DynamicArray<T, Architecture::CPU> globalArray_in) {
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

    T getGlobalValue(const MathVector<unsigned int, 3>& iP,
                     const unsigned int iC = 0) const {
      return (T) -1;
    }

    MathVector<T, NumberComponents> getGlobalVector(const MathVector<unsigned int, 3>& iP) const {
      return MathVector<T, NumberComponents>{(T) -1};
    }


  private:
    DEVICE HOST
    void setLocalHostValue(const MathVector<unsigned int, 3>& iP,
                           const T value,
                           const unsigned int iC) {
      localArrayHost[lSD::getIndex(iP, iC)] = value;
    }

    DEVICE HOST
    void setLocalHostVector(const MathVector<unsigned int, 3>& iP,
                            const MathVector<T, NumberComponents> vector) {
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        setLocalHostValue(iP, vector[iC], iC);
      }
    }

  };


  template <class T, unsigned int NumberComponents, DomainType initDomainType>
  class Field<T, NumberComponents, initDomainType, Architecture::CPU, true>
  : public Field<T, NumberComponents, initDomainType, Architecture::Generic, true> {
  protected:
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::localArrayDevice;

    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::globalArrayHost;
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::localArrayHost;

  public:
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::numberComponents;
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::IsWritten;

    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::fieldName;

    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::Field;

    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::globalData;
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::globalArray;
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::setGlobalField;
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::setGlobalValue;
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::setGlobalVector;
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::getGlobalValue;
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::getGlobalVector;

    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::localHostArray;
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::localDeviceArray;

    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::localHostData;
    using Field<T, NumberComponents, initDomainType,
                Architecture::Generic, true>::localDeviceData;


    DEVICE HOST
    T * RESTRICT localComputedData() {
      return localArrayHost.data();
    }


    DEVICE HOST
    void setLocalValue(const unsigned int index,
                       const T value,
                       const unsigned int iC = 0) {
      localArrayHost[lSD::getIndex(index, iC)] = value;
    }

    DEVICE HOST
    void setLocalVector(const unsigned int index,
                        const MathVector<T, NumberComponents> vector) {
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        setLocalValue(index, vector[iC], iC);
      }
    }

  };


  template <class T, unsigned int NumberComponents, DomainType initDomainType>
  class Field<T, NumberComponents, initDomainType, Architecture::GPU, true>
    : public Field<T, NumberComponents, initDomainType, Architecture::CPU, true> {
  protected:
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::localArrayDevice;

    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::globalArrayHost;
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::localArrayHost;

  public:
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::numberComponents;
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::IsWritten;

    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::fieldName;

    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::Field;

    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::globalData;
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::globalArray;
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::setGlobalField;
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::setGlobalValue;
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::setGlobalVector;
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::getGlobalValue;
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::getGlobalVector;

    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::localHostArray;
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::localDeviceArray;

    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::localHostData;
    using Field<T, NumberComponents, initDomainType,
                Architecture::CPU, true>::localDeviceData;

    DEVICE HOST
    T * RESTRICT localComputedData() {
      return localArrayDevice.data();
    }

    DEVICE HOST
    void setLocalValue(const unsigned int index,
                       const T value,
                       const unsigned int iC = 0) {
      localArrayDevice[lSD::getIndex(index, iC)] = value;
    }

    DEVICE HOST
    void setLocalVector(const unsigned int index,
                        const MathVector<T, NumberComponents> vector) {
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        setLocalValue(index, vector[iC], iC);
      }
    }

  };

}

#endif // FIELD_H
