#pragma once

#include "Helpers.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Options.h"

namespace lbm {

constexpr int globalLengthInt[3] = {globalLengthX,
                                    L::dimD > 1 ? globalLengthY : 1,
                                    L::dimD > 2 ? globalLengthZ : 1};

constexpr unsigned int globalLengthUInt[3] = {globalLengthX,
                                              L::dimD > 1 ? globalLengthY : 1,
                                              L::dimD > 2 ? globalLengthZ : 1};

constexpr ptrdiff_t globalLengthPtrdiff_t[3] = {
    globalLengthX, L::dimD > 1 ? ::lbm::globalLengthY : 1,
    L::dimD > 2 ? ::lbm::globalLengthZ : 1};

constexpr Position localLength = {globalLengthX / numProcs,
                                  L::dimD > 1 ? globalLengthY : 1,
                                  L::dimD > 2 ? globalLengthZ : 1};

/**
 * Domain defining space where DynamicArray lives and providing them
 * with a multi-dimensional index.
 *
 * @tparam latticeType type of lattice used.
 * @tparam domainType type of domain locality used.
 * @tparam partitionningType type of domain partitionning used.
 * @tparam memoryLayout type of memory layout used.
 */

template <DomainType domainType,
          PartitionningType partitionningType,
          MemoryLayout memoryLayout,
          unsigned int NumberComponents>
struct Domain {};

template <unsigned int NumberComponents>
struct Domain<DomainType::LocalSpace,
              PartitionningType::Generic,
              MemoryLayout::Generic,
              NumberComponents> {
 public:
  LBM_HOST LBM_DEVICE static inline
  constexpr Position pStart() {
    return Position({0, 0, 0});
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  constexpr Position pEnd() {
    return ProjectPadRealAndLeave1<unsigned int, L::dimD>::Do(
        {{globalLengthX / numProcs, globalLengthY, globalLengthZ}});
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  constexpr Position pLength() {
    return pEnd();
  }

  LBM_HOST LBM_DEVICE static inline unsigned int pVolume() {
    return pLength()[d::X] * pLength()[d::Y] * pLength()[d::Z];
  }

  LBM_HOST LBM_DEVICE static inline constexpr Position sStart() {
    return Position({0, 0, 0});
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  constexpr Position sEnd() {
    return localLength;
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  constexpr Position sLength() {
    return sEnd();
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  unsigned int sVolume() {
    return sLength()[d::X] * sLength()[d::Y] * sLength()[d::Z];
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  unsigned int getIndex(
      const Position& iP) {
    return pLength()[d::Z] * (pLength()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
  }
};

template <PartitionningType partitionningType, unsigned int NumberComponents>
struct Domain<DomainType::GlobalSpace, partitionningType,
              MemoryLayout::Generic, NumberComponents>
  : public Domain<DomainType::LocalSpace, PartitionningType::Generic,
                  MemoryLayout::Generic, NumberComponents> {
 private:
  using Base = Domain<DomainType::LocalSpace, PartitionningType::Generic,
                      MemoryLayout::Generic, NumberComponents>;

 public:
  LBM_HOST LBM_DEVICE static LBM_INLINE
  constexpr Position pStart() {
    return Position({0, 0, 0});
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  constexpr Position pEnd() {
    return ProjectAndLeave1<unsigned int, L::dimD>::Do(
        {{numProcs * Base::pLength()[d::X], Base::pLength()[d::Y],
          Base::pLength()[d::Z]}});
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  constexpr Position pLength() {
    return pEnd();
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  unsigned int pVolume() {
    return pLength()[d::X] * pLength()[d::Y] * pLength()[d::Z];
  }

  LBM_HOST LBM_DEVICE static inline
  Position pOffset(const MathVector<int, 3>& rank) {
    Position offsetR{{0}};
    for (auto iD = 0; iD < L::dimD; ++iD) {
      offsetR[iD] = (unsigned int)Base::pLength()[iD] * rank[iD];
    }
    return offsetR;
  }

  LBM_HOST LBM_DEVICE static inline constexpr Position sStart() {
    return Position({0, 0, 0});
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  constexpr Position sEnd() {
    return ProjectAndLeave1<unsigned int, L::dimD>::Do(
        {{globalLengthX, globalLengthY, globalLengthZ}});
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  constexpr Position sLength() {
    return sEnd();
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  unsigned int sVolume() {
    return sLength()[d::X] * sLength()[d::Y] * sLength()[d::Z];
  }

  LBM_HOST LBM_DEVICE static inline
  Position sOffset(const MathVector<int, 3>& rank) {
    Position offsetR{{0}};
    for (auto iD = 0; iD < L::dimD; ++iD) {
      offsetR[iD] = (unsigned int)Base::sLength()[iD] * rank[iD];
    }
    return offsetR;
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  unsigned int getIndex(const Position& iP) {
    const unsigned int indexLocal = Base::getIndex(
        {iP[d::X] - iP[d::X] / Base::length()[d::X] * Base::length()[d::X],
         iP[d::Y], iP[d::Z]});
    return iP[d::X] / Base::length()[d::X] * Base::volume() + indexLocal;
  }
};

template <>
struct Domain<DomainType::HaloSpace,
              PartitionningType::Generic,
              MemoryLayout::Generic,
              L::dimQ> : public Domain<DomainType::LocalSpace, PartitionningType::Generic,
                                       MemoryLayout::Generic, L::dimQ> {
 private:
  using Base = Domain<DomainType::LocalSpace, PartitionningType::Generic,
                      MemoryLayout::Generic,  L::dimQ>;

 public:
  LBM_HOST LBM_DEVICE static inline
  constexpr Position start() {
    return Position({0, 0, 0});
  }

  LBM_HOST LBM_DEVICE static inline
  constexpr Position end() {
    return Base::sLength() + L::halo() + L::halo();
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  constexpr Position length() {
    return Base::sLength() + L::halo() + L::halo();
  }

  LBM_HOST LBM_DEVICE static inline
  unsigned int volume() {
    return length()[d::X] * length()[d::Y] * length()[d::Z];
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  unsigned int getIndex(const Position& iP) {
    return length()[d::Z] * (length()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE
  unsigned int getIndexLocal(const Position& iP) {
    return Base::getIndex(iP - L::halo());
  }
 };

template <>
struct Domain<DomainType::HaloSpace,
              PartitionningType::Generic,
              MemoryLayout::AoS,
              L::dimQ> : public Domain<DomainType::HaloSpace,
                                       PartitionningType::Generic,
                                       MemoryLayout::Generic,
                                       L::dimQ> {
 private:
  using Base = Domain<DomainType::HaloSpace,
                      PartitionningType::Generic,
                      MemoryLayout::Generic,
                      L::dimQ>;

 public:
  using Base::end;
  using Base::getIndex;
  using Base::getIndexLocal;
  using Base::length;
  using Base::start;
  using Base::volume;

  LBM_HOST LBM_DEVICE static LBM_INLINE unsigned int getIndex(
      const Position& iP,
      const unsigned int iC) {
    return getIndex(iP) * L::dimQ + iC;
  }

  LBM_HOST LBM_DEVICE static LBM_INLINE unsigned int getIndex(
      const unsigned int index,
      const unsigned int iC) {
    return index * L::dimQ + iC;
  }
};

template <>
struct Domain<DomainType::HaloSpace,
              PartitionningType::Generic,
              MemoryLayout::SoA,
              L::dimQ> : public Domain<DomainType::HaloSpace,
                                       PartitionningType::Generic,
                                       MemoryLayout::Generic,
                                       L::dimQ> {
 private:
  using Base = Domain<DomainType::HaloSpace,
                      PartitionningType::Generic,
                      MemoryLayout::Generic,
                      L::dimQ>;

 public:
  using Base::end;
  using Base::getIndex;
  using Base::getIndexLocal;
  using Base::length;
  using Base::start;
  using Base::volume;

  LBM_HOST LBM_DEVICE LBM_INLINE static unsigned int getIndex(
      const Position& iP,
      const unsigned int iC) {
    return iC * volume() + getIndex(iP);
  }

  LBM_HOST LBM_DEVICE LBM_INLINE static unsigned int getIndex(
      const unsigned int index,
      const unsigned int iC) {
    return iC * volume() + index;
  }
};

template <unsigned int NumberComponents>
struct Domain<DomainType::BufferXSpace,
              PartitionningType::Generic,
              MemoryLayout::Generic,
              NumberComponents> : public Domain<DomainType::HaloSpace,
                                                PartitionningType::Generic,
                                                MemoryLayout::Generic,
                                                NumberComponents> {
 private:
  using Base = Domain<DomainType::HaloSpace,
                      PartitionningType::Generic,
                      MemoryLayout::Generic,
                      NumberComponents>;

 public:
  LBM_HOST LBM_DEVICE static inline constexpr Position start() {
    return Position({0, 0, 0});
  }

  LBM_HOST LBM_DEVICE static inline constexpr Position end() {
    return Position(
        {L::halo()[d::X],
         ProjectAndLeave1<unsigned int, L::dimD>::Do(Base::length())[d::Y],
         ProjectAndLeave1<unsigned int, L::dimD>::Do(Base::length())[d::Z]});
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

  LBM_HOST LBM_DEVICE static inline unsigned int getIndex(
      const Position& iP,
      const unsigned int iC) {
    return iC * volume() + getIndex(iP);
  }
};

using BaseDomain_ = Domain<DomainType::Generic,
                           PartitionningType::Generic,
                           MemoryLayout::Generic,
                           1>;

using gSD = Domain<DomainType::GlobalSpace,
                   PartitionningType::Generic,
                   MemoryLayout::Generic,
                   1>;

using lSD = Domain<DomainType::LocalSpace,
                   PartitionningType::Generic,
                   MemoryLayout::Generic,
                   1>;
using hSD =
    Domain<DomainType::HaloSpace, PartitionningType::Generic, memoryL, L::dimQ>;
using bXSD = Domain<DomainType::BufferXSpace,
                    PartitionningType::Generic,
                    MemoryLayout::Generic,
                    L::dimQ>;

}  // namespace lbm
