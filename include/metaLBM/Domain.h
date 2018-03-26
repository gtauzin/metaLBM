#ifndef DOMAIN_H
#define DOMAIN_H

#include "Options.h"
#include "Helpers.h"
#include "Lattice.h"
#include "MathVector.h"

namespace alloc {
  unsigned int numberElements;
}

namespace lbm {

  constexpr int globalLength[] = {globalLengthX, L::dimD>1 ? globalLengthY: 1,
                                  L::dimD>2 ? globalLengthZ: 1};

  /**
   * Domain defining space where DynamicArray lives and providing them
   * with a multi-dimensional index.
   *
   * @tparam latticeType type of lattice used.
   * @tparam domainType type of domain locality used.
   * @tparam partitionningType type of domain partitionning used.
   * @tparam memoryLayout type of memory layout used.
   */

  template <DomainType domainType, PartitionningType partitionningType,
            MemoryLayout memoryLayout, unsigned int NumberComponents>
    struct Domain {};

  template <unsigned int NumberComponents>
  struct Domain<DomainType::LocalSpace, PartitionningType::Generic,
                MemoryLayout::Generic, NumberComponents> {
  public:
    HOST DEVICE
    static inline constexpr Position pStart() {
      return Position({0, 0, 0});
    }

    HOST DEVICE
    static inline constexpr Position pEnd() {
      return ProjectPadRealAndLeave1<unsigned int, L::dimD>::Do({{globalLengthX/NPROCS,
            globalLengthY, globalLengthZ}});
    }

    HOST DEVICE
    static inline constexpr Position pLength() {
      return pEnd();
    }

    HOST DEVICE
    static inline unsigned int pVolume() {
      return pLength()[d::X]*pLength()[d::Y]*pLength()[d::Z];
    }


    HOST DEVICE
    static inline constexpr Position sStart() {
      return Position({0, 0, 0});
    }

    HOST DEVICE
    static inline constexpr Position sEnd() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do({{globalLengthX/NPROCS,
            globalLengthY, globalLengthZ}});
    }

    HOST DEVICE
    static inline constexpr Position sLength() {
      return sEnd();
    }

    HOST DEVICE
    static inline unsigned int sVolume() {
      return sLength()[d::X]*sLength()[d::Y]*sLength()[d::Z];
    }

    HOST DEVICE
    static inline unsigned int getIndex(const Position& iP) {
      return pLength()[d::Z] * (pLength()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    HOST DEVICE
    static inline unsigned int getIndex(const Position& iP,
                                        const unsigned int iC) {
      return iC * alloc::numberElements + getIndex(iP);
    }

    HOST DEVICE
    static inline unsigned int getIndex(const unsigned int index,
                                        const unsigned int iC) {
      return iC * alloc::numberElements + index;
    }

    HOST DEVICE
    static inline unsigned int getIndex(const unsigned int iC) {
      return iC * alloc::numberElements;
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
    HOST DEVICE
    static inline constexpr Position pStart() {
      return Position({0, 0, 0});
    }

    HOST DEVICE
    static inline constexpr Position pEnd() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do({{NPROCS * Base::pLength()[d::X], Base::pLength()[d::Y],Base::pLength()[d::Z]}});
    }

    HOST DEVICE
    static inline constexpr Position pLength() {
      return pEnd();
    }

    HOST DEVICE
    static inline unsigned int pVolume() {
      return pLength()[d::X]*pLength()[d::Y]*pLength()[d::Z];
    }

    HOST DEVICE
    static inline Position pOffset(const MathVector<int, 3>& rankMPI) {
      Position offsetR{{0}};
      for(auto iD = 0; iD < L::dimD; ++iD) {
        offsetR[iD] = (unsigned int) Base::pLength()[iD] * rankMPI[iD];
      }
      return offsetR;
    }

    HOST DEVICE
    static inline constexpr Position sStart() {
      return Position({0, 0, 0});
    }

    HOST DEVICE
    static inline constexpr Position sEnd() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do({{globalLengthX, globalLengthY,
              globalLengthZ}});
    }

    HOST DEVICE
    static inline constexpr Position sLength() {
      return sEnd();
    }

    HOST DEVICE
    static inline unsigned int sVolume() {
      return sLength()[d::X]*sLength()[d::Y]*sLength()[d::Z];
    }

    HOST DEVICE
    static inline Position sOffset(const MathVector<int, 3>& rankMPI) {
      Position offsetR{{0}};
      for(auto iD = 0; iD < L::dimD; ++iD) {
        offsetR[iD] = (unsigned int) Base::sLength()[iD] * rankMPI[iD];
      }
      return offsetR;
    }


    HOST DEVICE
    static inline unsigned int getIndex(const Position& iP) {
      const unsigned int indexLocal = Base::getIndex({iP[d::X] - iP[d::X]/Base::length()[d::X] * Base::length()[d::X], iP[d::Y], iP[d::Z]});
      return iP[d::X]/Base::length()[d::X] * Base::volume() + indexLocal;
    }

  };


  template <>
  struct Domain<DomainType::HaloSpace, PartitionningType::Generic,
                MemoryLayout::Generic, L::dimQ>
    : public Domain<DomainType::LocalSpace, PartitionningType::Generic,
                    MemoryLayout::Generic, L::dimQ> {
  private:
    using Base = Domain<DomainType::LocalSpace, PartitionningType::Generic,
                        MemoryLayout::Generic, L::dimQ>;

  public:
    HOST DEVICE
    static inline constexpr Position start() {
      return Position({0, 0, 0});
    }

    HOST DEVICE
    static inline constexpr Position end() {
      return Base::sLength() + L::halo() + L::halo();
    }

    HOST DEVICE
    static inline constexpr Position length() {
      return Base::sLength() + L::halo() + L::halo();
    }

    HOST DEVICE
    static inline unsigned int volume() {
      return length()[d::X]*length()[d::Y]*length()[d::Z];
    }

    HOST DEVICE
    static inline unsigned int getIndex(const Position& iP) {
      return length()[d::Z] * (length()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    HOST DEVICE
    static inline unsigned int getIndexLocal(const Position& iP) {
      return Base::getIndex(iP - L::halo());
    }

    HOST DEVICE
    static inline unsigned int getIndexLocal(const Position& iP,
                                             const unsigned int iC) {
      return iC*alloc::numberElements + getIndexLocal(iP);
    }

  };


  template <>
  struct Domain<DomainType::HaloSpace, PartitionningType::Generic,
                MemoryLayout::AoS, L::dimQ>
    : public Domain<DomainType::HaloSpace, PartitionningType::Generic,
                    MemoryLayout::Generic, L::dimQ> {
  private:
    using Base = Domain<DomainType::HaloSpace, PartitionningType::Generic,
                        MemoryLayout::Generic, L::dimQ>;

  public:
    using Base::start;
    using Base::end;
    using Base::length;
    using Base::volume;
    using Base::getIndex;
    using Base::getIndexLocal;

    HOST DEVICE
    static inline unsigned int getIndex(const Position& iP,
                                        const unsigned int iC) {
      return getIndex(iP) * L::dimQ + iC;
    }

    HOST DEVICE
    static inline unsigned int getIndex(const unsigned int index,
                                        const unsigned int iC) {
      return index * L::dimQ + iC;
    }

  };

  template <>
  struct Domain<DomainType::HaloSpace, PartitionningType::Generic,
                MemoryLayout::SoA, L::dimQ>
    : public Domain<DomainType::HaloSpace, PartitionningType::Generic,
                    MemoryLayout::Generic, L::dimQ> {
  private:
    using Base = Domain<DomainType::HaloSpace, PartitionningType::Generic,
                        MemoryLayout::Generic, L::dimQ>;

  public:
    using Base::start;
    using Base::end;
    using Base::length;
    using Base::volume;
    using Base::getIndex;
    using Base::getIndexLocal;

    HOST DEVICE
    static inline unsigned int getIndex(const Position& iP,
                                        const unsigned int iC) {
      return iC * volume() + getIndex(iP);
    }

    HOST DEVICE
    static inline unsigned int getIndex(const unsigned int index,
                                        const unsigned int iC) {
      return iC * volume() + index;
    }

  };


  template <unsigned int NumberComponents>
  struct Domain<DomainType::BufferXSpace, PartitionningType::Generic,
                MemoryLayout::Generic, NumberComponents>
    : public Domain<DomainType::HaloSpace, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents> {
  private:
    using Base = Domain<DomainType::HaloSpace, PartitionningType::Generic,
                        MemoryLayout::Generic, NumberComponents>;

  public:
    HOST DEVICE
    static inline constexpr Position start() {
      return Position({0, 0, 0});
    }

    HOST DEVICE
    static inline constexpr Position end() {
      return Position({L::halo()[d::X],
            ProjectAndLeave1<unsigned int, L::dimD>::Do(Base::length())[d::Y],
            ProjectAndLeave1<unsigned int, L::dimD>::Do(Base::length())[d::Z]});
    }

    HOST DEVICE
    static inline constexpr Position length() {
      return end();
    }

    HOST DEVICE
      static inline unsigned int volume() {
      return length()[d::X]*length()[d::Y]*length()[d::Z];
    }

    HOST DEVICE
    static inline unsigned int getIndex(const Position& iP) {
      return length()[d::Z] * (length()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    HOST DEVICE
    static inline unsigned int getIndex(const Position& iP,
                                        const unsigned int iC) {
      return iC * volume() + getIndex(iP);
    }
  };

  typedef Domain<DomainType::Generic, PartitionningType::Generic,
                 MemoryLayout::Generic, 1> BaseDomain_;

  typedef Domain<DomainType::GlobalSpace, PartitionningType::Generic,
                 MemoryLayout::Generic, 1> gSD;

  typedef Domain<DomainType::LocalSpace, PartitionningType::Generic,
                 MemoryLayout::Generic, 1> lSD;
  typedef Domain<DomainType::HaloSpace, PartitionningType::Generic,
                 memoryL, L::dimQ> hSD;
  typedef Domain<DomainType::BufferXSpace, PartitionningType::Generic,
                 MemoryLayout::Generic, L::dimQ> bXSD;

} // namespace lbm

#endif // DOMAIN_H
