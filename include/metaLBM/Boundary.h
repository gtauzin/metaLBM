#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "StaticArray.h"
#include "Domain.h"

namespace lbm {

  template<class T, BoundaryType boundaryType, AlgorithmType algorithmType>
  class Boundary{};

  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull> {
  private:
    MathVector<unsigned int, 3> iP_Origin;
    MathVector<unsigned int, 3> iP_Destination;

  protected:
    #pragma omp simd
    DEVICE HOST
    inline void applyX(T * RESTRICT f,
                       const MathVector<unsigned int, 3>& iP) {
      SCOREP_INSTRUMENT_ON("Boundary<T, boundaryType, algorithmType>::applyX")

      iP_Origin =  MathVector<unsigned int, 3>({L::halo()[d::X], iP[d::Y], iP[d::Z]});
      iP_Destination =  MathVector<unsigned int, 3>({L::halo()[d::X] + lD::length()[d::X], iP[d::Y], iP[d::Z]});

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
      }

      iP_Origin =  MathVector<unsigned int, 3>({L::halo()[d::X]+ lD::length()[d::X] -1, iP[d::Y], iP[d::Z]});
      iP_Destination =  MathVector<unsigned int, 3>({0, iP[d::Y], iP[d::Z]});

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
      }

    }

    #pragma omp simd
    DEVICE HOST
    inline void applyY(T * RESTRICT f,
                       const MathVector<unsigned int, 3>& iP) {
      SCOREP_INSTRUMENT_ON("Boundary<T, boundaryType, algorithmType>::applyY")

      iP_Origin =  MathVector<unsigned int, 3>({iP[d::X], L::halo()[d::Y], iP[d::Z]});
      iP_Destination =  MathVector<unsigned int, 3>({iP[d::X], L::halo()[d::Y] + lD::length()[d::Y], iP[d::Z]});

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
      }

      iP_Origin =  MathVector<unsigned int, 3>({iP[d::X],
            L::halo()[d::Y]+ lD::length()[d::Y] -1, iP[d::Z]});
      iP_Destination =  MathVector<unsigned int, 3>({iP[d::X], 0, iP[d::Z]});

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
      }
    }

    #pragma omp simd
    DEVICE HOST
    inline void applyZ(T * RESTRICT f,
                       const MathVector<unsigned int, 3>& iP) {
      SCOREP_INSTRUMENT_ON("Boundary<T, boundaryType, algorithmType>::applyZ")

      iP_Origin =  MathVector<unsigned int, 3>({iP[d::X], iP[d::Y], L::halo()[d::Z]});
      iP_Destination =  MathVector<unsigned int, 3>({iP[d::X], iP[d::Y],
            L::halo()[d::Z] + lD::length()[d::Z]});

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
      }

      iP_Origin =  MathVector<unsigned int, 3>({iP[d::X], iP[d::Y],
            L::halo()[d::Z] + lD::length()[d::Z] - 1});
      iP_Destination =  MathVector<unsigned int, 3>({iP[d::X], iP[d::Y], 0});

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
          f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
      }

    }

  };

  template<class T, AlgorithmType algorithmType>
  class Boundary<T, BoundaryType::Generic, algorithmType> {

  public:
    void apply() {}
  };

  template<class T, AlgorithmType algorithmType>
  class Boundary<T, BoundaryType::BounceBack_Halfway, algorithmType>
    : public Boundary<T, BoundaryType::Generic, algorithmType> {
  protected:

  public:
    void apply() {}

  };

  template<class T, AlgorithmType algorithmType>
  class Boundary<T, BoundaryType::Entropic, algorithmType>
    : public Boundary<T, BoundaryType::Generic, algorithmType> {
  protected:

  public:
    void apply() {}

  };

  typedef Boundary<dataT, boundaryT, algorithmT> Boundary_;

}

#endif // BOUNDARY_H