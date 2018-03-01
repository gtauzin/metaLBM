#ifndef DOMAIN_H
#define DOMAIN_H

#include "Options.h"
#include "Helpers.h"
#include "Lattice.h"
#include "MathVector.h"

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
  struct Domain<DomainType::LocalSpace, PartitionningType::Generic,
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
      return ProjectAndLeave1<unsigned int, L::dimD>::Do({lengthX_g/NPROCS,
            lengthY_g, lengthZ_g});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> length() {
      return end();
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


  template <PartitionningType partitionningType, unsigned int NumberComponents>
  struct Domain<DomainType::GlobalSpace, partitionningType,
                MemoryLayout::Generic, NumberComponents>
    : public Domain<DomainType::LocalSpace, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents> {
  private:
    using Base = Domain<DomainType::LocalSpace, PartitionningType::Generic,
                        MemoryLayout::Generic, NumberComponents>;

  public:
  public:
    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3>({0, 0, 0});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> end() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do({lengthX_g, lengthY_g, lengthZ_g});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> length() {
      return end();
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
        offsetR[iD] = (unsigned int) Base::length()[iD] * rankMPI[iD];
      }
      return offsetR;
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC = 0) {
      const unsigned int localLengthX = Base::length()[d::X];
      //if(iC != 0) std::cout << "iC " << iC << " of " << NumberComponents << std::endl;
      const unsigned int indexLocal = Base::getIndex({iP[d::X] - iP[d::X]/localLengthX * localLengthX, iP[d::Y], iP[d::Z]}, iC);
      return iP[d::X]/localLengthX * NumberComponents * Base::volume() + indexLocal;
    }
  };


  template <unsigned int NumberComponents>
  struct Domain<DomainType::HaloSpace, PartitionningType::Generic,
                MemoryLayout::Generic, NumberComponents>
    : public Domain<DomainType::LocalSpace, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents> {
  private:
    using Base = Domain<DomainType::LocalSpace, PartitionningType::Generic,
                        MemoryLayout::Generic, NumberComponents>;

  public:
    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> start() {
      return MathVector<unsigned int, 3>({0, 0, 0});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> end() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(Base::length()) + 2 * L::halo();
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> length() {
      return ProjectAndLeave1<unsigned int, L::dimD>::Do(Base::length()) + 2 * L::halo();
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
      return Base::getIndex(iP - L::halo(), iC);
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndexLocal(const MathVector<unsigned int, 3>& iP) {
      return Base::getIndex(iP - L::halo());
    }

  };


  template <unsigned int NumberComponents>
  struct Domain<DomainType::HaloSpace, PartitionningType::Generic,
                MemoryLayout::AoS, NumberComponents>
    : public Domain<DomainType::HaloSpace, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents> {
  private:
    using Base = Domain<DomainType::HaloSpace, PartitionningType::Generic,
                        MemoryLayout::Generic, NumberComponents>;

  public:
    using Base::start;
    using Base::end;
    using Base::length;
    using Base::volume;
    using Base::getIndex;
    using Base::getIndexLocal;

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
  struct Domain<DomainType::HaloSpace, PartitionningType::Generic,
                MemoryLayout::SoA, NumberComponents>
    : public Domain<DomainType::HaloSpace, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents> {
  private:
    using Base = Domain<DomainType::HaloSpace, PartitionningType::Generic,
                        MemoryLayout::Generic, NumberComponents>;

  public:
    using Base::start;
    using Base::end;
    using Base::length;
    using Base::volume;
    using Base::getIndex;
    using Base::getIndexLocal;

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
  struct Domain<DomainType::BufferXSpace, PartitionningType::Generic,
                MemoryLayout::Generic, NumberComponents>
    : public Domain<DomainType::HaloSpace, PartitionningType::Generic,
                    MemoryLayout::Generic, NumberComponents> {
  private:
    using Base = Domain<DomainType::HaloSpace, PartitionningType::Generic,
                        MemoryLayout::Generic, NumberComponents>;

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
            ProjectAndLeave1<unsigned int, L::dimD>::Do(Base::length())[d::Y],
            ProjectAndLeave1<unsigned int, L::dimD>::Do(Base::length())[d::Z]});
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> length() {
      return end();
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

  typedef Domain<DomainType::GlobalSpace, PartitionningType::Generic,
                 MemoryLayout::Generic, 1> gSD;
  typedef Domain<DomainType::GlobalSpace, PartitionningType::Generic,
                 MemoryLayout::Generic, L::dimD> gSDD;
  typedef Domain<DomainType::GlobalSpace, PartitionningType::Generic,
                 MemoryLayout::Generic, L::dimQ> gSQD;

  typedef Domain<DomainType::LocalSpace, PartitionningType::Generic,
                 MemoryLayout::Generic, 1> lSD;
  typedef Domain<DomainType::LocalSpace, PartitionningType::Generic,
                 MemoryLayout::Generic, L::dimD> lSDD;
  typedef Domain<DomainType::HaloSpace, PartitionningType::Generic,
                 memoryL, L::dimQ> hSD;
  typedef Domain<DomainType::BufferXSpace, PartitionningType::Generic,
                 MemoryLayout::Generic, L::dimQ> bXSD;

}

#endif // DOMAIN_H
