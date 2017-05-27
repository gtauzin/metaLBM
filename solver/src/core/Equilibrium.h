#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include <omp.h>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/sources/severity_feature.hpp>
namespace logging = boost::log;

#include "Options.h"
#include "Helpers.h"

namespace lbm {

  template<class T, LatticeType latticeType, EquilibriumType equilibriumType>
  class Equilibrium {};

  template <class T, LatticeType latticeType>
  class Equilibrium<T, latticeType, EquilibriumType::Basic> {
  protected:
    T density;
    MathVector<T, L::dimD> velocity;
    T velocity2;

  public:
    // TODO objects are aware of their types
    static constexpr EquilibriumType Type = equilibriumType;

    Equilibrium()
      : density( (T) 0)
      , velocity(MathVector<T, L::dimD>{{ (T) 0}})
      , velocity2(0)
    {}

    inline void setDensity(T density_in) {
      density = density_in;
    }

    inline void setVelocity(MathVector<T, L::dimD> velocity_in) {
      velocity = velocity_in;
      velocity2 = velocity_in.norm2();
    }

    inline void setVariables(const T density_in,
                             const MathVector<T, L::dimD> velocity_in) {
      setDensity(density_in);
      setVelocity(velocity_in);
  }

    #pragma omp declare simd
    template <unsigned int iQ>
    inline T computeEquilibrium() {
      T cu = L::celerity()[iQ].dot(velocity);

      T fEq_iQ = 1.0 +
        cu*L::inv_cs2 - 0.5 * velocity2*L::inv_cs2 + 0.5*Power<T, 2>(L::inv_cs2)*cu*cu
        - 0.5*Power<T, 2>(L::inv_cs2)*cu*velocity2 + Power<T, 3>(cu)*Power<T, 3>(L::inv_cs2)/6.0
        + 0.125*velocity2*velocity2*Power<T, 2>(L::inv_cs2) - 0.25*cu*cu*velocity2*Power<T, 3>(L::inv_cs2)
        + Power<T, 4>(cu)*Power<T, 4>(L::inv_cs2)/24.0;

      return density * L::weight()[iQ]*fEq_iQ;
    }
  };

  template <class T>
  class Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>
    : public Equilibrium<T, LatticeType::Generic, EquilibriumType::Basic> {
  private:
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Basic>::density;
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Basic>::velocity;
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Basic>::velocity2;

  public:
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Basic>::Equilibrium;
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Basic>::setDensity;
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Basic>::setVelocity;
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Basic>::setVariables;

    #pragma omp declare simd
    template <unsigned int iQ>
    inline T computeEquilibrium() {

      T fEq_iQ = 1.0;

      UnrolledFor<0, L::dimD>::Do([&] (int iD) {
          fEq_iQ *= (2.0 - sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
            * Power<T, L::celerity()[iQ][iD]>((2* velocity[iD]
                                               + sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
                                              /(1.0 - velocity[iD]));
        });

      return density * L::weight()[iQ]*fEq_iQ;
    }
  };

  template <class T>
  class Equilibrium<T, LatticeType::D2Q9, EquilibriumType::Basic>
    : public Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic> {
  private:
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::density;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::velocity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::velocity2;

  public:
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::Equilibrium;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::setDensity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::setVelocity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::setVariables;

    #pragma omp declare simd
    template <unsigned int iQ>
    inline T computeEquilibrium() {

      T fEq_iQ = 1.0;

      UnrolledFor<0, L::dimD>::Do([&] (int iD) {
          fEq_iQ *= (2.0 - sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
            * Power<T, L::celerity()[iQ][iD]>((2* velocity[iD]
                                                  + sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
                                                 /(1.0 - velocity[iD]));
        });

      return density * L::weight()[iQ]*fEq_iQ;
    }
  };

  template <class T>
  class Equilibrium<T, LatticeType::D3Q27, EquilibriumType::Basic>
      : public Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic> {
  private:
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::density;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::velocity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::velocity2;

  public:
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::Equilibrium;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::setDensity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::setVelocity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Basic>::setVariables;

    #pragma omp declare simd
    template <unsigned int iQ>
    inline T computeEquilibrium() {

      T fEq_iQ = 1.0;

      UnrolledFor<0, L::dimD>::Do([&] (int iD) {
          fEq_iQ *= (2.0 - sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
            * Power<T, L::celerity()[iQ][iD]>((2* velocity[iD]
                                                  + sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
                                                 /(1.0 - velocity[iD]));
        });

      return density * L::weight()[iQ]*fEq_iQ;
    }

  };

}

#endif // EQUILIBRIUM_H
