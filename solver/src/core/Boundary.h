#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include "Options.h"
#include "Domain.h"
#include "MathVector.h"
#include "StaticArray.h"

namespace lbm {

  template<class T, BoundaryType boundaryType>
  class Boundary{};

  template<class T, BoundaryType::Periodic>
  class Boundary {
  private:
    Domain<DomainType::Local> D;
    int origin;
    int destination;

  protected:

    template<unsigned int direction>
    inline void apply(T * __restrict__ f,
                      const MathVector<unsigned int, 3>& iP)  {
      // ABORT
    }

    template<>
    inline void apply<d::X>(T * __restrict__ f,
                            const MathVector<unsigned int, 3>& iP)  {

      origin = calculateIndex({D::start()[d::X], iP[d::Y], iP[d::Z]});
      destination = calculateIndex({D::end()[d::X], iP[d::Y], iP[d::Z]});

      srcBtoT  = idxL(L::hX, iY, iZ);
      destBtoT = idxL(L::hX + L::lX_l, iY, iZ);

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
                f[idxPop(destBtoT, iQ)] = f[idxPop(srcBtoT, iQ)];
        });

      srcBtoT  = idxL(L::hX + L::lX_l-1, iY, iZ);
      destBtoT = idxL(0, iY, iZ);

      UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
              f[idxPop(destBtoT, iQ)] = f[idxPop(srcBtoT, iQ)];
        });

    }

  public:

  };

  template<class T, BoundaryType boundaryType::Generic>
  class Boundary{

  public:

  };

  template<class T, BoundaryType::BounceBounceBack_Halfway>
  class Boundary : public Boundary<T, BoundaryType::Generic> {
  protected:

  public:

  };

  template<class T, BoundaryType::Entropic>
  class Boundary : public Boundary<T, BoundaryType::Generic> {
  protected:

  public:

  };


  template<class T>
  class Boundaries {
  private:
    StaticArray<Boundary<T, BoundaryType::Generic>, numberBoundaries> boundariesArray;

  public:
    Boundaries(const StaticArray<Boundary<T, BoundaryType::Generic>, numberBoundaries>&
               boundariesArray_in)
      : boundariesArray(boundariesArray_in)
    {}

    inline void apply(T * __restrict__ f){
      for(auto boundary : boundariesArray){
        boundary.apply(f);
      }
    }
  };

}

#endif // BOUNDARY_H
