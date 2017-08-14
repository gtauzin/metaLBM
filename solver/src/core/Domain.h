#ifndef DOMAIN_H
#define DOMAIN_H

#include "Input.h"
#include "Options.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Helpers.h"

namespace lbm {


  /**
   * Domain defining space where LocalizedField lives and providing them
   * with a multi-dimensional index.
   *
   * @tparam latticeType type of lattice used.
   * @tparam domainType type of domain locality used.
   * @tparam partitionningType type of domain partitionning used.
   * @tparam memoryLayout type of memory layout used.
   */

  template <LatticeType latticeType, DomainType domainType,
            PartitionningType partitionningType,
            MemoryLayout memoryLayout>
  struct Domain {};

  template <LatticeType latticeType, PartitionningType partitionningType,
            MemoryLayout memoryLayout>
  struct Domain<latticeType, DomainType::Local, partitionningType, memoryLayout>  {
  private:
      typedef Lattice<int, latticeType> lattice;

  public:
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3>{0, 0, 0};
    }

    static inline constexpr MathVector<unsigned int, 3> end() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l);
    }

    static inline constexpr MathVector<unsigned int, 3> length() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l);
    }

    static inline constexpr unsigned int volume() {
      return length()[d::X]*length()[d::Y]*length()[d::Z];
    }

    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP) {
      return length()[d::Z] * (length()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      return iC * volume() + getIndex(iP);
    }

  };



  template <LatticeType latticeType, PartitionningType partitionningType,
            MemoryLayout memoryLayout>
  struct Domain<latticeType, DomainType::Global, partitionningType, memoryLayout>
    : public Domain<latticeType, DomainType::Local, partitionningType, memoryLayout> {
  private:
    typedef Lattice<int, latticeType> lattice;

  public:
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3>{0, 0, 0};
    }

    static inline constexpr MathVector<unsigned int, 3> end() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_g);
    }

    static inline constexpr MathVector<unsigned int, 3> length() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_g);
    }

    static inline constexpr unsigned int volume() {
      return length()[d::X]*length()[d::Y]*length()[d::Z];
    }

    static inline MathVector<unsigned int, 3> offset
      (const MathVector<int, 3>& rankMPI = MathVector<int, 3>{0, 0, 0}) {
      MathVector<unsigned int, 3> offsetR{{0}};
      UnrolledFor<0, lattice::dimD>::Do([&] (unsigned int iD) {
          offsetR[iD] = (unsigned int) length_g[iD]*rankMPI[iD];
        });

      return offsetR;
  }

    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP) {
      const int indexLocal = Domain<latticeType, DomainType::Local,
                                    partitionningType, memoryLayout>::getIndex({iP[d::X] - iP[d::X]/length()[d::X] * length()[d::X], iP[d::Y], iP[d::Z]});
      return iP[d::X]/length()[d::X]
        * Domain<latticeType, DomainType::Local,
                 partitionningType, memoryLayout>::volume()
        + indexLocal;
    }

    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      return iC * volume() + getIndex(iP);
    }

    static inline unsigned int getIndex(const unsigned int index,
                                        const unsigned int iC) {
      return iC * volume() + index;
    }

  };


  template <LatticeType latticeType, PartitionningType partitionningType>
  struct Domain<latticeType, DomainType::Halo,
                partitionningType, MemoryLayout::Generic>
    : public Domain<latticeType, DomainType::Local,
                    partitionningType, MemoryLayout::Generic> {
  private:
    typedef Lattice<int, latticeType> lattice;

  public:
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3>{0, 0, 0};
    }

    static inline constexpr MathVector<unsigned int, 3> end() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l) + 2 * lattice::halo();
    }

    static inline constexpr MathVector<unsigned int, 3> length() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l) + 2 * lattice::halo();
    }

    static inline constexpr unsigned int volume() {
      return length()[d::X]*length()[d::Y]*length()[d::Z];
    }

    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP) {
      return length()[d::Z] * (length()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    static inline unsigned int getIndexLocal(const MathVector<unsigned int, 3>& iP,
                                             const unsigned int iC) {
      return Domain<latticeType, DomainType::Local,
                    partitionningType,
                    MemoryLayout::Generic>::getIndex(iP - lattice::halo(), iC);
    }

    static inline unsigned int getIndexLocal(const MathVector<unsigned int, 3>& iP) {
      return Domain<latticeType, DomainType::Local,
                    partitionningType,
                    MemoryLayout::Generic>::getIndex(iP - lattice::halo());
    }

  };


  template <LatticeType latticeType, PartitionningType partitionningType>
  struct Domain<latticeType, DomainType::Halo,
                partitionningType, MemoryLayout::AoS>
    : public Domain<latticeType, DomainType::Halo,
                    partitionningType, MemoryLayout::Generic> {
  private:
    typedef Lattice<int, latticeType> lattice;

  public:
    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::start;
    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::end;

    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::length;
    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::volume;

    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::getIndex;
    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::getIndexLocal;


    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      return getIndex(iP) * lattice::dimQ + iC;
    }

    static inline unsigned int getIndex(const unsigned int index,
                                        const unsigned int iC) {
      return index * lattice::dimQ + iC;
    }

  };

  template <LatticeType latticeType, PartitionningType partitionningType>
  struct Domain<latticeType, DomainType::Halo, partitionningType, MemoryLayout::SoA>
    : public Domain<latticeType, DomainType::Halo, partitionningType, MemoryLayout::Generic> {
  private:
    typedef Lattice<int, latticeType> lattice;

  public:
    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::start;
    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::end;

    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::length;
    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::volume;

    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::getIndex;
    using Domain<latticeType, DomainType::Halo,
                 partitionningType, MemoryLayout::Generic>::getIndexLocal;


    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      return iC * volume() + getIndex(iP);
    }

    static inline unsigned int getIndex(const unsigned int index,
                                        const unsigned int iC) {
      return iC * volume() + index;
    }

  };

  template <LatticeType latticeType, PartitionningType partitionningType,
            MemoryLayout memoryLayout>
  struct Domain<latticeType, DomainType::BufferX,
                partitionningType, memoryLayout>
    : public Domain<latticeType, DomainType::Halo,
                    partitionningType, MemoryLayout::Generic> {
  private:
    typedef Lattice<int, latticeType> lattice;

  public:
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3>{0, 0, 0};
    }

    static inline constexpr MathVector<unsigned int, 3> end() {
      return {L::halo()[d::X],
              ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l)[d::Y]+2*L::halo()[d::Y],
              ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l)[d::Z]+2*L::halo()[d::Z]};
    }

    static inline constexpr MathVector<unsigned int, 3> length() {
      return {L::halo()[d::X],
          ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l)[d::Y]+2*L::halo()[d::Y],
              ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l)[d::Z]+2*L::halo()[d::Z]};
    }

    static inline constexpr unsigned int volume() {
      return length()[d::X]*length()[d::Y]*length()[d::Z];
    }

    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP) {
      return length()[d::Z] * (length()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      return iC * volume() + getIndex(iP);
    }

  };


  typedef Domain<latticeT, DomainType::Global, partitionningT, memoryL> gD;
  typedef Domain<latticeT, DomainType::Local, partitionningT, memoryL> lD;
  typedef Domain<latticeT, DomainType::Halo, partitionningT, memoryL> hD;
  typedef Domain<latticeT, DomainType::BufferX, partitionningT, memoryL> bXD;

}

#endif // DOMAIN_H
