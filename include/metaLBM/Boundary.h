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
    inline void applyX(T * __restrict__ f,
                       const MathVector<unsigned int, 3>& iP) {
      SCOREP_INSTRUMENT("Boundary<T, boundaryType, algorithmType>::applyX")

      iP_Origin = {L::halo()[d::X], iP[d::Y], iP[d::Z]};
      iP_Destination = {L::halo()[d::X] + lD::length()[d::X], iP[d::Y], iP[d::Z]};

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
      });

      iP_Origin = {L::halo()[d::X]+ lD::length()[d::X] -1, iP[d::Y], iP[d::Z]};
      iP_Destination = {0, iP[d::Y], iP[d::Z]};

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
      });

    }

    inline void applyY(T * __restrict__ f,
                       const MathVector<unsigned int, 3>& iP) {
      SCOREP_INSTRUMENT("Boundary<T, boundaryType, algorithmType>::applyY")

      iP_Origin = {iP[d::X], L::halo()[d::Y], iP[d::Z]};
      iP_Destination = {iP[d::X], L::halo()[d::Y] + lD::length()[d::Y], iP[d::Z]};

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
        });

      iP_Origin = {iP[d::X], L::halo()[d::Y]+ lD::length()[d::Y] -1, iP[d::Z]};
      iP_Destination = {iP[d::X], 0, iP[d::Z]};

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
        });
    }

    inline void applyZ(T * __restrict__ f,
                       const MathVector<unsigned int, 3>& iP) {
      SCOREP_INSTRUMENT("Boundary<T, boundaryType, algorithmType>::applyZ")

      iP_Origin = {iP[d::X], iP[d::Y], L::halo()[d::Z]};
      iP_Destination = {iP[d::X], iP[d::Y], L::halo()[d::Z] + lD::length()[d::Z]};

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
      });

      iP_Origin = {iP[d::X], iP[d::Y], L::halo()[d::Z] + lD::length()[d::Z] - 1};
      iP_Destination = {iP[d::X], iP[d::Y], 0};

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
          f[hD::getIndex(iP_Destination, iQ)] = f[hD::getIndex(iP_Origin, iQ)];
      });

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
