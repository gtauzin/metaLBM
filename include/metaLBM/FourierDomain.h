#pragma once

#include "Domain.h"

namespace lbm {

/**
 * Domain defining space where DynamicArray lives and providing them
 * with a multi-dimensional index.
 *
 * @tparam latticeType type of lattice used.
 * @tparam domainType type of domain locality used.
 * @tparam partitionningType type of domain partitionning used.
 * @tparam memoryLayout type of memory layout used.
 */

template <PartitionningType partitionningType, unsigned int NumberComponents>
struct Domain<DomainType::LocalFourier,
              partitionningType,
              MemoryLayout::Generic,
              NumberComponents> {
 public:
  LBM_HOST LBM_DEVICE static inline constexpr Position start() {
    return Position({0, 0, 0});
  }

  LBM_HOST LBM_DEVICE static inline constexpr Position end() {
    return ProjectPadComplexAndLeave1<unsigned int, L::dimD>::Do(
        lSD::sLength());
  }

  LBM_HOST LBM_DEVICE static inline constexpr Position length() {
    return end();
  }

  LBM_HOST LBM_DEVICE static inline unsigned int volume() {
    return length()[d::X] * length()[d::Y] * length()[d::Z];
  }

  LBM_HOST LBM_DEVICE static inline unsigned int getIndex(const Position& iP) {
    return length()[d::Z] * (length()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
  }
};

template <PartitionningType partitionningType>
struct Domain<DomainType::GlobalFourier,
              partitionningType,
              MemoryLayout::Generic,
              1> : public Domain<DomainType::LocalFourier,
                                 partitionningType,
                                 MemoryLayout::Generic,
                                 1> {
 private:
  using Base = Domain<DomainType::LocalFourier,
                      partitionningType,
                      MemoryLayout::Generic,
                      1>;

 public:
  LBM_HOST LBM_DEVICE static inline constexpr Position start() {
    return Position({0, 0, 0});
  }

  LBM_HOST LBM_DEVICE static inline constexpr Position end() {
    return ProjectPadComplexAndLeave1<unsigned int, L::dimD>::Do(
        gSD::sLength());
  }

  LBM_HOST LBM_DEVICE static inline constexpr Position length() {
    return end();
  }

  LBM_HOST LBM_DEVICE static inline constexpr unsigned int volume() {
    return length()[d::X] * length()[d::Y] * length()[d::Z];
  }

  LBM_HOST LBM_DEVICE static inline constexpr unsigned int maxWaveNumber() {
    return arrayMax(globalLengthInt) / 2;
  }

  LBM_HOST LBM_DEVICE static inline Position offset(const MathVector<int, 3>& rank) {
    Position offsetR{{0}};
    for (auto iD = 0; iD < L::dimD; ++iD) {
      offsetR[iD] = (unsigned int)Base::length()[iD] * rank[iD];
    }

    return offsetR;
  }

  LBM_HOST LBM_DEVICE static inline unsigned int getIndex(const Position& iP) {
    const unsigned int localLengthX = Base::length()[d::X];
    const unsigned int localVolume = Base::volume();

    const unsigned int indexLocal = Base::getIndex(
        Position({iP[d::X] - iP[d::X] / localLengthX * localLengthX, iP[d::Y],
                  iP[d::Z]}));
    return iP[d::X] / localLengthX * localVolume + indexLocal;
  }
};

typedef Domain<DomainType::GlobalFourier,
               partitionningT,
               MemoryLayout::Generic,
               1>
    gFD;

typedef Domain<DomainType::LocalFourier,
               PartitionningType::Generic,
               MemoryLayout::Generic,
               1>
    lFD;

}  // namespace lbm
