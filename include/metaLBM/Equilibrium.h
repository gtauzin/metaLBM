#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include <omp.h>

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "Lattice.h"
#include "Helpers.h"

namespace lbm {

  template<class T, LatticeType latticeType, EquilibriumType equilibriumType>
  class Equilibrium {};

  template <class T, LatticeType latticeType>
  class Equilibrium<T, latticeType, EquilibriumType::Incompressible> {
  protected:
    T density;
    MathVector<T, L::dimD> velocity;
    T velocity2;

  public:
    static constexpr EquilibriumType Type = EquilibriumType::Incompressible;

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
      SCOREP_INSTRUMENT_OFF("Equilibrium<T, L, EquilibriumType::Incompressible>::setVariables")

      setDensity(density_in);
      setVelocity(velocity_in);
  }

    #pragma omp declare simd
    inline T calculate(const unsigned int iQ) const {
      SCOREP_INSTRUMENT_OFF("Equilibrium<T, latticeType, EquilibriumType::Incompressible>::calculate")

      T cu = L::celerity()[iQ].dot(velocity);

      T fEq_iQ = 1.0
        + cu*L::inv_cs2 - 0.5 * velocity2*L::inv_cs2
        + 0.5*Power<T, 2>::Do(L::inv_cs2)*cu*cu
        - 0.5*Power<T, 2>::Do(L::inv_cs2)*cu*velocity2
        + Power<T, 3>::Do(cu)*Power<T, 3>::Do(L::inv_cs2)/6.0
        + 0.125*velocity2*velocity2*Power<T, 2>::Do(L::inv_cs2)
        - 0.25*cu*cu*velocity2*Power<T, 3>::Do(L::inv_cs2)
        + Power<T, 4>::Do(cu)*Power<T, 4>::Do(L::inv_cs2)/24.0;

      return density * L::weight()[iQ]*fEq_iQ;
    }
  };

  template <class T>
  class Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>
    : public Equilibrium<T, LatticeType::Generic, EquilibriumType::Incompressible> {
  private:
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Incompressible>::density;
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Incompressible>::velocity;
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Incompressible>::velocity2;

  public:
    static constexpr EquilibriumType Type = EquilibriumType::Incompressible;

    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Incompressible>::Equilibrium;
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Incompressible>::setDensity;
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Incompressible>::setVelocity;
    using Equilibrium<T, LatticeType::Generic, EquilibriumType::Incompressible>::setVariables;

    #pragma omp declare simd
    inline T calculate(const unsigned int iQ) const {
      SCOREP_INSTRUMENT_OFF("Equilibrium<T, D1Q3, EquilibriumType::Incompressible>::calculate")

      T fEq_iQ = 1.0;

      UnrolledFor<0, L::dimD>::Do([&] (int iD) {
          fEq_iQ *= (2.0 - sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
            * PowerBase((2* velocity[iD]
                        + sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
                       /(1.0 - velocity[iD]), L::celerity()[iQ][iD]);
        });

      return density * L::weight()[iQ]*fEq_iQ;
    }
  };

  template <class T>
  class Equilibrium<T, LatticeType::D2Q9, EquilibriumType::Incompressible>
    : public Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible> {
  private:
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::density;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::velocity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::velocity2;

  public:
    static constexpr EquilibriumType Type = EquilibriumType::Incompressible;

    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::Equilibrium;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::setDensity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::setVelocity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::setVariables;

    #pragma omp declare simd
    inline T calculate(const unsigned int iQ) const {
      SCOREP_INSTRUMENT_OFF("Equilibrium<T, D2Q9, EquilibriumType::Incompressible>::calculate")

      T fEq_iQ = 1.0;

      UnrolledFor<0, L::dimD>::Do([&] (int iD) {
          fEq_iQ *= (2.0 - sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
            * PowerBase((2* velocity[iD]
                         + sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
                        /(1.0 - velocity[iD]), uiL::celerity()[iQ][iD]);
        });


      return density * L::weight()[iQ]*fEq_iQ;
    }
  };

  template <class T>
  class Equilibrium<T, LatticeType::D3Q27, EquilibriumType::Incompressible>
      : public Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible> {
  private:
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::density;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::velocity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::velocity2;

  public:
    static constexpr EquilibriumType Type = EquilibriumType::Incompressible;

    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::Equilibrium;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::setDensity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::setVelocity;
    using Equilibrium<T, LatticeType::D1Q3, EquilibriumType::Incompressible>::setVariables;

    #pragma omp declare simd
    inline T calculate(const unsigned int iQ) const {
      SCOREP_INSTRUMENT_OFF("Equilibrium<T, D3Q27, EquilibriumType::Incompressible>::calculate")

      T fEq_iQ = (T) 1;

      UnrolledFor<0, L::dimD>::Do([&] (int iD) {
          fEq_iQ *= (2.0 - sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
            * PowerBase((2* velocity[iD]
                         + sqrt(1.0 + 3.0*velocity[iD]*velocity[iD]))
                        /(1.0 - velocity[iD]), L::celerity()[iQ][iD]);
        });

      return density * L::weight()[iQ]*fEq_iQ;
    }

  };

  typedef Equilibrium<dataT, latticeT, equilibriumT> Equilibrium_;

}

#endif // EQUILIBRIUM_H
