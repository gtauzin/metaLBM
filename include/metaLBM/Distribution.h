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

  template <class T, Architecture architecture>
    class Distribution {};

  template <class T>
  class Distribution<T, Architecture::CPU>
    : public Field<T, L::dimQ, Architecture::CPU, true> {
  protected:
    using Field<T, L::dimQ, Architecture::CPU, true>::localArrayHost;
    using Field<T, L::dimQ, Architecture::CPU, true>::localArrayDevice;

    DynamicArray<T, Architecture::CPU> haloArrayHost;
    DynamicArray<T, Architecture::GPU> haloArrayDevice;

  public:
    using Field<T, L::dimQ, Architecture::CPU, true>::globalArray;
    using Field<T, L::dimQ, Architecture::CPU, true>::setGlobalField;

    Distribution(const std::string& fieldName_in)
      : Field<T, L::dimQ, Architecture::CPU, true>(fieldName_in)
      , haloArrayHost(hD::volume()*L::dimQ)
      , haloArrayDevice(hD::volume()*L::dimQ)
    {}

    Distribution(const std::string& fieldName_in,
                 const DynamicArray<T, Architecture::CPU>& globalArray_in)
      : Field<T, L::dimQ, Architecture::CPU, true>(fieldName_in, globalArray_in)
      , haloArrayHost(hD::volume()*L::dimQ)
      , haloArrayDevice(hD::volume()*L::dimQ)
    {}

    DEVICE HOST
    void swapHalo(Distribution<T, Architecture::CPU>& distribution_in) {
      haloArrayHost.swap(distribution_in.haloHostArray());
    }

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

    void packLocal() {
      MathVector<unsigned int, 3> iP;

      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for(unsigned int iZ = lD::start()[d::Z]+L::halo()[d::Z];
          iZ < lD::end()[d::Z]+L::halo()[d::Z]; ++iZ) {
        for(unsigned int iY = lD::start()[d::Y]+L::halo()[d::Y];
            iY < lD::end()[d::Y]+L::halo()[d::Y]; ++iY) {
         #pragma omp simd
          for(unsigned int iX = lD::start()[d::X]+L::halo()[d::X];
              iX < lD::end()[d::X]+L::halo()[d::X]; ++iX) {
            iP =  MathVector<unsigned int, 3>({iX, iY, iZ});

            for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
                localArrayHost[hD::getIndexLocal(iP, iQ)]
                    = haloArrayHost[hD::getIndex(iP, iQ)];
            }
          }
        }
      }
    }

    void unpackLocal() {
      MathVector<unsigned int, 3> iP;

      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for(unsigned int iZ = lD::start()[d::Z]+L::halo()[d::Z];
          iZ < lD::end()[d::Z]+L::halo()[d::Z]; ++iZ) {
        for(unsigned int iY = lD::start()[d::Y]+L::halo()[d::Y];
            iY < lD::end()[d::Y]+L::halo()[d::Y]; ++iY) {
         #pragma omp simd
          for(unsigned int iX = lD::start()[d::X]+L::halo()[d::X];
              iX < lD::end()[d::X]+L::halo()[d::X]; ++iX) {
            iP =  MathVector<unsigned int, 3>({iX, iY, iZ});

            for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
              haloArrayHost[hD::getIndex(iP, iQ)]
                  = localArrayHost[lD::getIndex(iP-L::halo(), iQ)];
            }

          }
        }
      }
   }

  };


  template <class T>
  class Distribution<T, Architecture::GPU>
    : public Field<T, L::dimQ, Architecture::GPU, true> {
  private:
    using Field<T, L::dimQ, Architecture::GPU, true>::localArrayHost;
    using Field<T, L::dimQ, Architecture::GPU, true>::localArrayDevice;

    DynamicArray<T, Architecture::CPU> haloArrayHost;
    DynamicArray<T, Architecture::GPU> haloArrayDevice;

  public:
    using Field<T, L::dimQ, Architecture::GPU, true>::globalArray;
    using Field<T, L::dimQ, Architecture::GPU, true>::setGlobalField;

    Distribution(const std::string& fieldName_in)
      : Field<T, L::dimQ, Architecture::GPU, true>(fieldName_in)
      , haloArrayHost(hD::volume()*L::dimQ)
      , haloArrayDevice(hD::volume()*L::dimQ)
      {}

    Distribution(const std::string& fieldName_in,
                 const DynamicArray<T, Architecture::CPU>& globalArray_in)
      : Field<T, L::dimQ, Architecture::GPU, true>(fieldName_in, globalArray_in)
      , haloArrayHost(hD::volume()*L::dimQ)
      , haloArrayDevice(hD::volume()*L::dimQ)
      {}

    DEVICE HOST
    void swapHalo(Distribution<T, Architecture::GPU>& distribution_in) {
      haloArrayDevice.swap(distribution_in.haloDeviceArray());
    }

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
      for(unsigned int iZ = lD::start()[d::Z]+L::halo()[d::Z];
          iZ < lD::end()[d::Z]+L::halo()[d::Z]; ++iZ) {
        for(unsigned int iY = lD::start()[d::Y]+L::halo()[d::Y];
            iY < lD::end()[d::Y]+L::halo()[d::Y]; ++iY) {
         #pragma omp simd
          for(unsigned int iX = lD::start()[d::X]+L::halo()[d::X];
              iX < lD::end()[d::X]+L::halo()[d::X]; ++iX) {
            iP = MathVector<unsigned int, 3>({iX, iY, iZ});

            for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
                localArrayHost[hD::getIndexLocal(iP, iQ)]
                    = haloArrayHost[hD::getIndex(iP, iQ)];
            }
          }
        }
      }
    }

    void unpackLocal() {
      MathVector<unsigned int, 3> iP;

      #pragma omp parallel for schedule(static) num_threads(NTHREADS)
      for(unsigned int iZ = lD::start()[d::Z]+L::halo()[d::Z];
          iZ < lD::end()[d::Z]+L::halo()[d::Z]; ++iZ) {
        for(unsigned int iY = lD::start()[d::Y]+L::halo()[d::Y];
            iY < lD::end()[d::Y]+L::halo()[d::Y]; ++iY) {
         #pragma omp simd
          for(unsigned int iX = lD::start()[d::X]+L::halo()[d::X];
              iX < lD::end()[d::X]+L::halo()[d::X]; ++iX) {
            iP = MathVector<unsigned int, 3>({iX, iY, iZ});

            for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
              haloArrayHost[hD::getIndex(iP, iQ)]
                  = localArrayHost[lD::getIndex(iP-L::halo(), iQ)];
            }

          }
        }
      }
   }

  };

}

#endif // DISTRIBUTION_H
