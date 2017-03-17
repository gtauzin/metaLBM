#ifndef LATTICE_H
#define LATTICE_H

#include <array>
#include <numeric>
#include <memory>

#include <boost/align/aligned_allocator.hpp>
#include <vector>
template<class T, std::size_t Alignment = 1>
#ifdef ALIGN_MEMORY
  using vector = std::vector<T,
  boost::alignment::aligned_allocator<T, Alignment> >;
#else
  using vector = std::vector<T>;
#endif

#include "input.h"
#include "commons.h"

namespace lbm {

  template <class T, LatticeType L>
  class Lattice {
  public:
    vector<T, CACHE_LINE> f_distribution;

    Lattice(const vector<T, CACHE_LINE>& f_distribution_in)
      : f_distribution(f_distribution_in)
    {}

  };

  template <class T, LatticeType L>
  class Field {
  public:
    vector<T, CACHE_LINE> nextDensity;
    vector<MathVector<T, P::dimD>, CACHE_LINE> nextVelocity;
    vector<T, CACHE_LINE> nextAlpha;
    vector<T, CACHE_LINE> nextDistribution;
    vector<T, CACHE_LINE> previousDensity;
    vector<MathVector<T, P::dimD>, CACHE_LINE> previousVelocity;

    Field(const unsigned int size_in)
      : nextDensity(size_in)
      , nextVelocity{size_in, MathVector<T, P::dimD>()}
      , nextAlpha(size_in)
      , nextDistribution(P::dimQ*size_in)
      , previousDensity(size_in)
      , previousVelocity{size_in, MathVector<T, P::dimD>()}
      {}

    Field(const vector<T, CACHE_LINE>& nextDensity_in,
          const vector<MathVector<T, P::dimD>, CACHE_LINE> nextVelocity_in,
          const vector<T, CACHE_LINE>& nextAlpha_in,
          const vector<T, CACHE_LINE>& nextDistribution_in,
          const vector<T, CACHE_LINE>& previousDensity_in,
          const vector<MathVector<T, P::dimD>, CACHE_LINE> previousVelocity_in)
      : nextDensity(nextDensity_in)
      , nextVelocity(nextVelocity_in)
      , nextAlpha(nextAlpha_in)
      , nextDistribution(nextDistribution_in)
      , previousDensity(previousDensity_in)
      , previousVelocity(previousVelocity_in)
    {}

    virtual int idxF(const int iX, const int iY, const int iZ) = 0;

  };

  template <class T, LatticeType L>
    class LocalField : public Field<T, L> {
  public:
  LocalField()
    : Field<T, L>(s_l())
      {}

  LocalField(const vector<T, CACHE_LINE>& nextDensity_in,
             const vector<MathVector<T, P::dimD>, CACHE_LINE> nextVelocity_in,
             const vector<T, CACHE_LINE>& nextAlpha_in,
             const vector<T, CACHE_LINE>& nextDistribution_in,
             const vector<T, CACHE_LINE>& previousDensity_in,
             const vector<MathVector<T, P::dimD>, CACHE_LINE> previousVelocity_in)
    : Field<T, L>(nextDensity_in, nextVelocity_in,
                  nextAlpha_in, nextDistribution_in,
                  previousDensity_in, previousVelocity_in)
      {}

    inline int idxF(const int iX, const int iY, const int iZ) {
      return idx_lF(iX, iY, iZ);
    }

    inline T getMass() {
      return std::accumulate(this->nextDistribution.cbegin(),
                             this->nextDistribution.cend(), 0.);
    }

    inline void packDistribution(const Lattice<T, L>& l) {
      for(int iZ = P::hZ; iZ < P::hZ + P::lZ_l; ++iZ) {
        for(int iY = P::hY; iY < P::hY + P::lY_l; ++iY) {
          for(int iX = P::hX; iX < P::hX + P::lX_l; ++iX) {

            int idx_L = idxL(iX, iY, iZ);
            int idx_F = idx_inF(iX, iY, iZ);
            UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
                this->nextDistribution[idxPop_lF(idx_F, iQ)] = l.f_distribution[idxPop(idx_L, iQ)];
              });
              }
        }
      }
    }

    inline vector<T, CACHE_LINE> unpackDistribution() {
      vector<T, CACHE_LINE> f_distributionR(P::dimQ
                                            *(P::lX_l+2*P::hX)
                                            *(P::lY_l+2*P::hY)
                                            *(P::lZ_l+2*P::hZ), (T)(-10));
      for(int iZ = P::hZ; iZ < P::hZ + P::lZ_l; ++iZ) {
        for(int iY = P::hY; iY < P::hY + P::lY_l; ++iY) {
          for(int iX = P::hX; iX < P::hX + P::lX_l; ++iX) {

            int idx_F = idx_inF(iX, iY, iZ);
            int idx_L = idxL(iX, iY, iZ);
            UnrolledFor<0, P::dimQ>::Do([&] (int iQ) {
                f_distributionR[idxPop(idx_L, iQ)] = this->nextDistribution[idxPop_lF(idx_F, iQ)];
              });

          }
        }
      }
      return f_distributionR;
    }

  };

  template <class T, LatticeType L>
    class GlobalField : public Field<T, L> {
  public:
    GlobalField()
      : Field<T, L>(s_g())
      {}

    GlobalField(const vector<T, CACHE_LINE>& nextDensity_in,
                const vector<MathVector<T, P::dimD>, CACHE_LINE> nextVelocity_in,
                const vector<T, CACHE_LINE>& nextAlpha_in,
                const vector<T, CACHE_LINE>& nextDistribution_in,
                const vector<T, CACHE_LINE>& previousDensity_in,
                const vector<MathVector<T, P::dimD>, CACHE_LINE> previousVelocity_in)
      : Field<T, L>(nextDensity_in, nextVelocity_in,
                    nextAlpha_in, nextDistribution_in,
                    previousDensity_in, previousVelocity_in)
      {}

    inline int idxF(const int iX, const int iY, int iZ) {
      return idx_gF(iX, iY, iZ);
    }
  };

}
#endif // LATTICE_H
