#ifndef DOMAIN_H
#define DOMAIN_H

#include "Input.h"
#include "Options.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Helpers.h"

namespace lbm {


  /**
   * Index required to provide Distributions with multi-dimensional index.
   *
   * @tparam memoryLayout type of memory layout used.
   */


  template <Lattice lattice, DomainType domainType,
            PartionningType partitionningType,
            MemoryLayout memoryLayout>
  struct Domain {);

  template <Lattice lattice, PartitionningType partitionningType,
            MemoryLayout memoryLayout = MemoryLayout::Generic>
  struct Domain<lattice, DomainType::Generic, partionningType, memoryLayout> {
    static inline constexpr MathVector<unsigned int, 3> start() {
      \\ ABORT
      return MathVector<unsigned int, 3> endR{{0}};
    }

    static inline constexpr MathVector<unsigned int, 3> end() {
      \\ ABORT
      return MathVector<unsigned int, 3> endR{{0}};
    }

    static inline constexpr MathVector<unsigned int, 3> length() {
      \\ ABORT
        return end() - start();
    }

    static constexpr unsigned int volume = length()[d::X]*length()[d::Y]*length()[d::Z];

    static unsigned int index;

    static inline unsigned int updateIndex(const MathVector<unsigned int, 3>& iP) {
      index = length()[d::Z] * (length()[d::Z] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    inline unsigned int getIndex(const unsigned int iC) {
      return iC * volume + index;
    }

  };

  template <Lattice lattice, PartitionningType partitionningType,
            MemoryLayout memoryLayout = MemoryLayout::Generic>
  struct Domain<lattice, DomainType::Local, partionningType, memoryLayout>
    : public Domain<lattice, DomainType::Generic, partionningType, memoryLayout> {
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3> endR{{0}};
    }

    static inline constexpr MathVector<unsigned int, 3> end() {
      constexpr MathVector<unsigned int, 3> endR{{1}};
      UnrolledFor<0, lattice::dimD>::Do([&] (int iD) {
          endR[iD] = length_g()[iD]/process()[iD];
        });

      return endR;
    }

    using Domain<lattice, DomainType::Generic, partionningType, memoryLayout>::length;
    using Domain<lattice, DomainType::Generic, partionningType, memoryLayout>::volume;

    using Domain<lattice, DomainType::Generic, partionningType, memoryLayout>::index;

    static inline unsigned int updateIndex(const MathVector<unsigned int, 3>& iP) {
      const int iPLocal = {iP[d::X] - iP[d::X]/length()[d::X] * length()[d::X],
                           iP[d::Y], iP[d::Z]};
      return iP[d::X]/length()[d::X]
        * Domain<lattice, DomainType::Local, partitionningType>::volume
        + Domain<lattice, DomainType::Local, partitionningType>::updateIndex(iPLocal);

      index = length()[d::Z] * (length()[d::Z] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

  };

  template <Lattice lattice,
            PartitionningType partitionningType = ParitionningType::Generic,
            MemoryLayout memoryLayout = MemoryLayout::Generic>
  struct Domain<lattice, DomainType::Global, partionningType, memoryLayout>
    : public Domain<lattice, DomainType::Local, partionningType, memoryLayout> {
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3> endR{{0}};
    }

    static inline constexpr MathVector<unsigned int, 3> end() {
      constexpr MathVector<unsigned int, 3> endR{{1}};
      UnrolledFor<0, lattice::dimD>::Do([&] (int iD) {
          endR[iD] = length_g()[iD];
        });

      return endR;
    }

    static inline MathVector<unsigned int, 3> offset
    (const MathVector<unsigned int, 3> rankMPI = MathVector<unsigned int, 3>{{0}}) {
      constexpr MathVector<unsigned int, 3> offsetR{{0}};
      UnrolledFor<0, lattice::dimD>::Do([&] (int iD) {
          offsetR[iD] = length_g()[iD]*rankMPI[iD];
        });

      return offsetR;
  }


    using Domain<lattice, DomainType::Local, partionningType, memoryLayout>::length;
    using Domain<lattice, DomainType::Local, partionningType, memoryLayout>::volume;

    using Domain<lattice, DomainType::Local, partionningType, memoryLayout>::index;
    using Domain<lattice, DomainType::Local, partionningType, memoryLayout>::updateIndex;

  };

  template <Lattice lattice, PartitionningType partitionningType>
  struct Domain<lattice, DomainType::Halo, partionningType, MemoryLayout::Generic>
    : public Domain<lattice, DomainType::Generic, partionningType, memoryLayout> {
    static inline constexpr MathVector<unsigned int, 3> start() {
      return L::halo();
    }

    static inline constexpr MathVector<unsigned int, 3> end() {
      constexpr MathVector<unsigned int, 3> endR{{1}};
      UnrolledFor<0, lattice::dimD>::Do([&] (int iD) {
          endR[iD] = length_g()[iD]/process()[iD] + 2* L::halo();
        });

      return endR;
    }

    using Domain<lattice, DomainType::Generic, partionningType, memoryLayout>::length;
    using Domain<lattice, DomainType::Generic, partionningType, memoryLayout>::volume;

    using Domain<lattice, DomainType::Generic, partionningType, memoryLayout>::index;
    using Domain<lattice, DomainType::Generic, partionningType, memoryLayout>::updateIndex;

  };


  template <Lattice lattice, PartitionningType partitionningType>
  struct Domain<lattice, DomainType::Halo, partionningType, MemoryLayout::AoS>
    : public Domain<lattice, DomainType::Halo, partionningType, MemoryLayout::Generic> {

    using Domain<lattice, DomainType::Halo, partionningType, MemoryLayout::Generic>::start;
    using Domain<lattice, DomainType::Halo, partionningType, MemoryLayout::Generic>::end;

    using Domain<lattice, DomainType::Halo, partionningType, MemoryLayout::Generic>::length;
    using Domain<lattice, DomainType::Halo, partionningType, MemoryLayout::Generic>::volume;

    using Domain<lattice, DomainType::Halo, partionningType, memoryLayout>::index;

    static inline unsigned int updateIndex(const MathVector<unsigned int, 3>& iP) {
      index = length()[d::Z] * (length()[d::Z] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    inline unsigned int getIndex(const unsigned int iC) {
      return index * lattice::dimQ + iC;
    }


    };

  template <Lattice lattice, PartitionningType partitionningType>
  struct Domain<lattice, DomainType::Halo, partionningType, MemoryLayout::SoA>
    : public Domain<lattice, DomainType::Halo, partionningType, MemoryLayout::Generic> {

    using Domain<lattice, DomainType::Halo,
                 partionningType, MemoryLayout::Generic>::start;
    using Domain<lattice, DomainType::Halo,
                 partionningType, MemoryLayout::Generic>::end;

    using Domain<lattice, DomainType::Halo,
                 partionningType, MemoryLayout::Generic>::length;
    using Domain<lattice, DomainType::Halo,
                 partionningType, MemoryLayout::Generic>::volume;

    using Domain<lattice, DomainType::Halo,
                 partionningType, MemoryLayout::Generic>::index;

    static inline unsigned int updateIndex(const MathVector<unsigned int, 3>& iP) {
      index = length()[d::Z] * (length()[d::Z] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    inline unsigned int getIndex(const unsigned int iC) {
      return iC * domainVolume + index;
    }

    };

}

#endif // DOMAIN_H
