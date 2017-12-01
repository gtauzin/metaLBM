#ifndef DOMAIN_H
#define DOMAIN_H

#include "Options.h"
#include "Helpers.h"
#include "MathVector.h"
#include "Lattice.h"

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

  template <DomainType domainType, PartitionningType partitionningType,
            MemoryLayout memoryLayout, unsigned int NumberComponents>
  struct Domain {};

  template <unsigned int NumberComponents>
  struct Domain<DomainType::Generic, PartitionningType::Generic,
                MemoryLayout::Generic, NumberComponents> {
  public:
    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3>({0, 0, 0});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> end() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l());
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> length() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l());
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr unsigned int volume() {
      return length()[d::X]*length()[d::Y]*length()[d::Z];
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP) {
      return length()[d::Z] * (length()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const unsigned int index) {
      return index;
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      return iC * volume() + getIndex(iP);
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const unsigned int index,
                                        const unsigned int iC) {
      return iC * volume() + index;
    }

  };

  template <unsigned int NumberComponents>
  struct Domain<DomainType::Local, PartitionningType::Generic,
                MemoryLayout::Generic, NumberComponents> {
  public:
    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3>({0, 0, 0});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> end() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l());
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> length() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l());
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr unsigned int volume() {
      return length()[d::X]*length()[d::Y]*length()[d::Z];
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP) {
      return length()[d::Z] * (length()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const unsigned int index) {
      return index;
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      return iC * volume() + getIndex(iP);
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const unsigned int index,
                                        const unsigned int iC) {
      return iC * volume() + index;
    }

  };


  template <PartitionningType partitionningType>
  struct Domain<DomainType::Global, partitionningType,
                MemoryLayout::Generic, 1>
    : public Domain<DomainType::Local, partitionningType,
                    MemoryLayout::Generic, 1>
  {
  public:
    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3>({0, 0, 0});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> end() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_g());
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> length() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_g());
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr unsigned int volume() {
      return length()[d::X]*length()[d::Y]*length()[d::Z];
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline MathVector<unsigned int, 3> offset(const MathVector<int, 3>& rankMPI) {
      MathVector<unsigned int, 3> offsetR{{0}};
      for(unsigned int iD = 0; iD < L::dimD; ++iD) {
          offsetR[iD] = (unsigned int) length_g()[iD]*rankMPI[iD];
      }

      return offsetR;
  }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP) {
      const unsigned int localLengthX = Domain<DomainType::Local, PartitionningType::Generic,
                                               MemoryLayout::Generic, 1>::length()[d::X];
      const unsigned int localVolume = Domain<DomainType::Local, PartitionningType::Generic,
                                              MemoryLayout::Generic, 1>::volume();

      const unsigned int indexLocal
        = Domain<DomainType::Local, PartitionningType::Generic,
                 MemoryLayout::Generic, 1>::getIndex(MathVector<unsigned int, 3>({iP[d::X] - iP[d::X]/localLengthX * localLengthX, iP[d::Y], iP[d::Z]}));
      return iP[d::X]/localLengthX * localVolume + indexLocal;
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      const unsigned int localLengthX = Domain<DomainType::Local, PartitionningType::Generic,
                                               MemoryLayout::Generic, 1>::length()[d::X];
      const unsigned int localVolume = Domain<DomainType::Local, PartitionningType::Generic,
                                              MemoryLayout::Generic, 1>::volume();

      const unsigned int indexLocal
        = Domain<DomainType::Local, PartitionningType::Generic,
                 MemoryLayout::Generic, 1>::getIndex(MathVector<unsigned int, 3>({iP[d::X] - iP[d::X]/localLengthX * localLengthX, iP[d::Y], iP[d::Z]}), iC);
      return iP[d::X]/localLengthX * localVolume + indexLocal;
    }

  };

  template <PartitionningType partitionningType, unsigned int NumberComponents>
  struct Domain<DomainType::Global, partitionningType,
                MemoryLayout::Generic, NumberComponents>
    : public Domain<DomainType::Global, partitionningType,
                    MemoryLayout::Generic, 1>
  {
  public:
    using Domain<DomainType::Global, partitionningType,
                 MemoryLayout::Generic, 1>::start;
    using Domain<DomainType::Global, partitionningType,
                 MemoryLayout::Generic, 1>::end;
    using Domain<DomainType::Global, partitionningType,
                 MemoryLayout::Generic, 1>::length;
    using Domain<DomainType::Global, partitionningType,
                 MemoryLayout::Generic, 1>::volume;
    using Domain<DomainType::Global, partitionningType,
                 MemoryLayout::Generic, 1>::offset;

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP) {
      const unsigned int localLengthX = Domain<DomainType::Local, PartitionningType::Generic,
                                               MemoryLayout::Generic, NumberComponents>::length()[d::X];
      const unsigned int localVolume = Domain<DomainType::Local, PartitionningType::Generic,
                                              MemoryLayout::Generic, NumberComponents>::volume();

      const unsigned int indexLocal = Domain<DomainType::Local, PartitionningType::Generic,
                                             MemoryLayout::Generic, NumberComponents>::getIndex({iP[d::X] - iP[d::X]/localLengthX * localLengthX, iP[d::Y], iP[d::Z]});
      return iP[d::X]/localLengthX * NumberComponents * localVolume + indexLocal;
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      const unsigned int localLengthX = Domain<DomainType::Local, PartitionningType::Generic,
                                               MemoryLayout::Generic, NumberComponents>::length()[d::X];
      const unsigned int localVolume = Domain<DomainType::Local, PartitionningType::Generic,
                                              MemoryLayout::Generic, NumberComponents>::volume();

      const unsigned int indexLocal = Domain<DomainType::Local, PartitionningType::Generic,
                                             MemoryLayout::Generic, NumberComponents>::getIndex({iP[d::X] - iP[d::X]/localLengthX * localLengthX, iP[d::Y], iP[d::Z]}, iC);
      return iP[d::X]/localLengthX * NumberComponents * localVolume + indexLocal;
    }
  };



  template <unsigned int NumberComponents>
  struct Domain<DomainType::Halo, PartitionningType::Generic,
                MemoryLayout::Generic, NumberComponents>
    : public Domain<DomainType::Local, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents> {
  public:
    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3>({0, 0, 0});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> end() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l()) + 2 * L::halo();
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> length() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l()) + 2 * L::halo();
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr unsigned int volume() {
      return length()[d::X]*length()[d::Y]*length()[d::Z];
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP) {
      return length()[d::Z] * (length()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndexLocal(const MathVector<unsigned int, 3>& iP,
                                             const unsigned int iC) {
      return Domain<DomainType::Local, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents>::getIndex(iP - L::halo(), iC);
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndexLocal(const MathVector<unsigned int, 3>& iP) {
      return Domain<DomainType::Local, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents>::getIndex(iP - L::halo());
    }

  };


  template <unsigned int NumberComponents>
  struct Domain<DomainType::Halo, PartitionningType::Generic,
                MemoryLayout::AoS, NumberComponents>
    : public Domain<DomainType::Halo, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents> {
  public:
    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::start;
    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::end;

    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::length;
    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::volume;

    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::getIndex;
    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::getIndexLocal;


    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      return getIndex(iP) * L::dimQ + iC;
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const unsigned int index,
                                        const unsigned int iC) {
      return index * L::dimQ + iC;
    }

  };

  template <unsigned int NumberComponents>
  struct Domain<DomainType::Halo, PartitionningType::Generic,
                MemoryLayout::SoA, NumberComponents>
    : public Domain<DomainType::Halo, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents> {
  public:
    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::start;
    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::end;

    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::length;
    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::volume;

    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::getIndex;
    using Domain<DomainType::Halo, PartitionningType::Generic,
                 MemoryLayout::Generic, NumberComponents>::getIndexLocal;


    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      return iC * volume() + getIndex(iP);
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const unsigned int index,
                                        const unsigned int iC) {
      return iC * volume() + index;
    }

  };

  template <unsigned int NumberComponents>
  struct Domain<DomainType::BufferX, PartitionningType::Generic,
                MemoryLayout::Generic, NumberComponents>
    : public Domain<DomainType::Halo, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents> {
  public:
    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3>({0, 0, 0});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> end() {
      return MathVector<unsigned int, 3>({L::halo()[d::X],
              ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l())[d::Y]+2*L::halo()[d::Y],
            ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l())[d::Z]+2*L::halo()[d::Z]});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> length() {
      return MathVector<unsigned int, 3>({L::halo()[d::X],
          ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l())[d::Y]+2*L::halo()[d::Y],
            ProjectAndLeave1<unsigned int, L::dimD>::Do(length_l())[d::Z]+2*L::halo()[d::Z]});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr unsigned int volume() {
      return length()[d::X]*length()[d::Y]*length()[d::Z];
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP) {
      return length()[d::Z] * (length()[d::Y] * iP[d::X] + iP[d::Y]) + iP[d::Z];
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      return iC * volume() + getIndex(iP);
    }

  };

  typedef Domain<DomainType::Global, partitionningT,
                 MemoryLayout::Generic, 1> gD;
  typedef Domain<DomainType::Global, partitionningT,
                 MemoryLayout::Generic, L::dimD> gDD;
  typedef Domain<DomainType::Global, partitionningT,
                 MemoryLayout::Generic, L::dimQ> gQD;
  typedef Domain<DomainType::Local, PartitionningType::Generic,
                 MemoryLayout::Generic, 1> lD;
  typedef Domain<DomainType::Halo, PartitionningType::Generic,
                 memoryL, L::dimQ> hD;
  typedef Domain<DomainType::BufferX, PartitionningType::Generic,
                 MemoryLayout::Generic, L::dimQ> bXD;

}

#endif // DOMAIN_H
