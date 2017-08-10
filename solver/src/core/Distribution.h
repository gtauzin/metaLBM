#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <iostream>

#include "Options.h"
#include "Field.h"
#include "Domain.h"
#include "Lattice.h"

namespace lbm {

  template <class T>
  class Distribution
    : public Field<T, L::dimQ, true> {
  protected:
    using Field<T, L::dimQ, true>::globalField;
    using Field<T, L::dimQ, true>::localField;
    LocalizedField<T, L::dimQ> haloField;

  public:
    using Field<T, L::dimQ, true>::setGlobalField;

    Distribution(const std::string& fieldName_in)
      : Field<T, L::dimQ, true>(fieldName_in)
      , haloField(fieldName_in, hD::volume())
      {}

    void swapHalo(Distribution<T>& distribution_in) {
      /* MathVector<unsigned int, 3> iP; */
      /* for (unsigned int iY = L::halo()[d::Y]+lD::start()[d::Y]; iY < L::halo()[d::Y]+lD::end()[d::Y]; ++iY) { */
      /*   for (unsigned int iX = L::halo()[d::X]+lD::start()[d::X]; iX < L::halo()[d::X]+lD::end()[d::X]; ++iX) { */
      /*     iP = {iX, iY, 0}; */
      /*     std::cout << "iP: " << iP << std::endl; */
      /*     for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) { */
      /*       std::cout << "iQ: " << iQ */
      /*                 << ", " << distribution_in.fieldName << ": " */
      /*                 << distribution_in.haloData()[hD::getIndex(iP, iQ)] */
      /*                 << ",  " <<  haloField.fieldName << ": " */
      /*                 << haloField[hD::getIndex(iP, iQ)] << std::endl; */
      /*     } */
      /*   } */
      /* } */

      haloField.swap(distribution_in.getHaloField());
    }

    LocalizedField<T, L::dimQ>& getHaloField() {
      return haloField;
    }

    void setHaloField(const unsigned int index, T value) {
      haloField.setField(index, value);
    }

    T * __restrict__ haloData(const unsigned int iC = 0) {
      return haloField.data(iC);
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
            iP = {iX, iY, iZ};

            UnrolledFor<0, L::dimQ>::Do([&] (unsigned int iQ) {
                localField[hD::getIndexLocal(iP, iQ)]
                    = haloField[hD::getIndex(iP, iQ)];
            });
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
            iP = {iX, iY, iZ};

            UnrolledFor<0, L::dimQ>::Do([&] (unsigned int iQ) {
                haloField[hD::getIndex(iP, iQ)]
                  = localField[hD::getIndexLocal(iP, iQ)];
            });

          }
        }
      }
   }


  };
}

#endif // DISTRIBUTION_H
