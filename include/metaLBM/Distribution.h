#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <iostream>

#include "Options.h"
#include "Commons.h"
#include "Field.h"
#include "Domain.h"
#include "Lattice.h"
#include "DynamicArray.cuh"

namespace lbm {

  template <class T, DomainType initDomainType, Architecture architecture>
    class Distribution {};

  template <class T, DomainType initDomainType>
  class Distribution<T, initDomainType, Architecture::CPU>
    : public Field<T, L::dimQ, initDomainType, Architecture::CPU, true> {
  protected:
    using Field<T, L::dimQ, initDomainType, Architecture::CPU, true>::localArrayHost;
    using Field<T, L::dimQ, initDomainType, Architecture::CPU, true>::localArrayDevice;

    DynamicArray<T, Architecture::CPU> haloArrayHost;
    DynamicArray<T, Architecture::GPU> haloArrayDevice;

  public:
    using Field<T, L::dimQ, initDomainType, Architecture::CPU, true>::fieldName;
    using Field<T, L::dimQ, initDomainType, Architecture::CPU, true>::globalArray;
    using Field<T, L::dimQ, initDomainType, Architecture::CPU, true>::setGlobalField;

    Distribution(const std::string& fieldName_in)
      : Field<T, L::dimQ, initDomainType, Architecture::CPU, true>(fieldName_in)
      , haloArrayHost(hSD::volume()*L::dimQ)
      , haloArrayDevice(hSD::volume()*L::dimQ)
    {}

    Distribution(const std::string& fieldName_in,
                 const DynamicArray<T, Architecture::CPU>& globalArray_in)
      : Field<T, L::dimQ, initDomainType,
              Architecture::CPU, true>(fieldName_in, globalArray_in)
      , haloArrayHost(hSD::volume()*L::dimQ)
      , haloArrayDevice(hSD::volume()*L::dimQ)
    {}

    DynamicArray<T, Architecture::CPU>& haloHostArray() {
      return haloArrayHost;
    }

    DEVICE HOST
    void setHaloField(const unsigned int index, const T value) {
      haloArrayHost[index] = value;
    }

    DEVICE HOST
    T * RESTRICT haloComputedData() {
      return haloArrayHost.data();
    }

    DEVICE HOST
    DynamicArray<T, Architecture::GPU>& haloDeviceArray() {
      return haloArrayDevice;
    }

    /* void packLocal(const MathVector<unsigned int, 3>& iP) { */
    /*   for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) { */
    /*     localArrayHost[hSD::getIndexLocal(iP, iQ)] */
    /*       = haloArrayHost[hSD::getIndex(iP, iQ)]; */
    /*   } */
    /* } */

    void packLocal() {
      MathVector<unsigned int, 3> iP;

      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for(unsigned int iZ = lSD::start()[d::Z]+L::halo()[d::Z];
          iZ < lSD::end()[d::Z]+L::halo()[d::Z]; ++iZ) {
        for(unsigned int iY = lSD::start()[d::Y]+L::halo()[d::Y];
            iY < lSD::end()[d::Y]+L::halo()[d::Y]; ++iY) {
         #pragma omp simd
          for(unsigned int iX = lSD::start()[d::X]+L::halo()[d::X];
              iX < lSD::end()[d::X]+L::halo()[d::X]; ++iX) {
            iP =  MathVector<unsigned int, 3>({iX, iY, iZ});

            for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
                localArrayHost[hSD::getIndexLocal(iP, iQ)]
                    = haloArrayHost[hSD::getIndex(iP, iQ)];
            }
          }
        }
      }
    }


    void unpackLocal() {
      MathVector<unsigned int, 3> iP;

      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for(unsigned int iZ = lSD::start()[d::Z]+L::halo()[d::Z];
          iZ < lSD::end()[d::Z]+L::halo()[d::Z]; ++iZ) {
        for(unsigned int iY = lSD::start()[d::Y]+L::halo()[d::Y];
            iY < lSD::end()[d::Y]+L::halo()[d::Y]; ++iY) {
         #pragma omp simd
          for(unsigned int iX = lSD::start()[d::X]+L::halo()[d::X];
              iX < lSD::end()[d::X]+L::halo()[d::X]; ++iX) {
            iP =  MathVector<unsigned int, 3>({iX, iY, iZ});

            for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
              haloArrayHost[hSD::getIndex(iP, iQ)]
                  = localArrayHost[lSD::getIndex(iP-L::halo(), iQ)];
            }

          }
        }
      }
   }

  };


  template <class T, DomainType initDomainType>
  class Distribution<T, initDomainType, Architecture::GPU>
    : public Field<T, L::dimQ, initDomainType, Architecture::GPU, true> {
  private:
    using Base = Field<T, L::dimQ, initDomainType, Architecture::GPU, true>;
    using Field<T, L::dimQ, initDomainType, Architecture::GPU, true>::localArrayHost;
    using Field<T, L::dimQ, initDomainType, Architecture::GPU, true>::localArrayDevice;

    DynamicArray<T, Architecture::CPU> haloArrayHost;
    DynamicArray<T, Architecture::GPU> haloArrayDevice;

  public:
    using Field<T, L::dimQ, initDomainType, Architecture::GPU, true>::globalArray;
    using Field<T, L::dimQ, initDomainType, Architecture::GPU, true>::setGlobalField;

    Distribution(const std::string& fieldName_in)
      : Field<T, L::dimQ, initDomainType, Architecture::GPU, true>(fieldName_in)
      , haloArrayHost(hSD::volume()*L::dimQ)
      , haloArrayDevice(hSD::volume()*L::dimQ)
      {}

    Distribution(const std::string& fieldName_in,
                 const DynamicArray<T, Architecture::CPU>& globalArray_in)
      : Field<T, L::dimQ, initDomainType,
              Architecture::GPU, true>(fieldName_in, globalArray_in)
      , haloArrayHost(hSD::volume()*L::dimQ)
      , haloArrayDevice(hSD::volume()*L::dimQ)
      {}

    DynamicArray<T, Architecture::CPU>& haloHostArray() {
      return haloArrayHost;
    }

    DEVICE HOST
    void setHaloField(const unsigned int index, const T value) {
      haloArrayDevice[index] = value;
    }

    DEVICE HOST
    T * RESTRICT haloComputedData() {
      return haloArrayDevice.data();
    }

    DEVICE HOST
    DynamicArray<T, Architecture::GPU>& haloDeviceArray() {
      return haloArrayDevice;
    }


    void packLocal() {
      MathVector<unsigned int, 3> iP;

      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for(unsigned int iZ = lSD::start()[d::Z]+L::halo()[d::Z];
          iZ < lSD::end()[d::Z]+L::halo()[d::Z]; ++iZ) {
        for(unsigned int iY = lSD::start()[d::Y]+L::halo()[d::Y];
            iY < lSD::end()[d::Y]+L::halo()[d::Y]; ++iY) {
         #pragma omp simd
          for(unsigned int iX = lSD::start()[d::X]+L::halo()[d::X];
              iX < lSD::end()[d::X]+L::halo()[d::X]; ++iX) {
            iP = MathVector<unsigned int, 3>({iX, iY, iZ});

            for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
                localArrayHost[hSD::getIndexLocal(iP, iQ)]
                    = haloArrayHost[hSD::getIndex(iP, iQ)];
            }
          }
        }
      }
    }

    // TODO use Computation.cuh

    void unpackLocal() {
      MathVector<unsigned int, 3> iP;

      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for(unsigned int iZ = lSD::start()[d::Z]+L::halo()[d::Z];
          iZ < lSD::end()[d::Z]+L::halo()[d::Z]; ++iZ) {
        for(unsigned int iY = lSD::start()[d::Y]+L::halo()[d::Y];
            iY < lSD::end()[d::Y]+L::halo()[d::Y]; ++iY) {
         #pragma omp simd
          for(unsigned int iX = lSD::start()[d::X]+L::halo()[d::X];
              iX < lSD::end()[d::X]+L::halo()[d::X]; ++iX) {
            iP = MathVector<unsigned int, 3>({iX, iY, iZ});

            for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
              haloArrayHost[hSD::getIndex(iP, iQ)]
                  = localArrayHost[lSD::getIndex(iP-L::halo(), iQ)];
            }

          }
        }
      }
   }

  };

}

#endif // DISTRIBUTION_H
