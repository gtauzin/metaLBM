#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <iostream>

#include "Options.h"
#include "Commons.h"
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

    Distribution(const std::string& fieldName_in,
                 const LocalizedField<T, L::dimQ>& globalField_in)
      : Field<T, L::dimQ, true>(fieldName_in, globalField_in)
      , haloField(fieldName_in, hD::volume())
      {}

    T * RESTRICT haloData(const unsigned int iC = 0) {
      return haloField.data(iC);
    }

    void swapHalo(Distribution<T>& distribution_in) {
      haloField.swap(distribution_in.getHaloField());
    }

    LocalizedField<T, L::dimQ>& getHaloField() {
      return haloField;
    }

    void setHaloField(const unsigned int index, const T value) {
      haloField[index] = value;
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
                  = localField[lD::getIndex(iP-L::halo(), iQ)];
            });

            /* MathVector<int, 3> rankMPI{1, 0, 0}; */
            /* MPI_Comm_rank(MPI_COMM_WORLD, &rankMPI[d::X]); */
            /* if(rankMPI[d::X] == 1 && iP == MathVector<unsigned int, 3>({3,3,0})) { */
            /*   std::cout << "In unpackLocal" << std::endl; */

            /*   for(int iQ = 0; iQ < L::dimQ; ++iQ) { */
            /*     std::cout << "localField(" << iQ << "): " */
            /*               << localField[lD::getIndex(iP-L::halo(), iQ)] */
            /*               << ", globalField(" << iQ << "): " */
            /*               << globalField[gQD::getIndex(gD::offset({1,0,0})+iP-L::halo(), iQ)] */
            /*               << std::endl; */
            /*   } */
            /* } */


          }
        }
      }
   }


  };
}

#endif // DISTRIBUTION_H
