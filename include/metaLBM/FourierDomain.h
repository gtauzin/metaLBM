#ifndef FOURIERDOMAIN_H
#define FOURIERDOMAIN_H

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


  template <unsigned int NumberComponents>
  struct Domain<DomainType::LocalFourier, PartitionningType::Generic,
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
      return ProjectAndLeaveHalfAnd1<unsigned int, L::dimD>::Do(length_l());
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> length() {
      return ProjectAndLeaveHalfAnd1<unsigned int, L::dimD>::Do(length_l());
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
  struct Domain<DomainType::GlobalFourier, partitionningType,
                MemoryLayout::Generic, 1>
    : public Domain<DomainType::LocalFourier, partitionningType,
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
      return ProjectAndLeaveHalfAnd1<unsigned int, L::dimD>::Do(length_g());
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline constexpr MathVector<unsigned int, 3> length() {
      return ProjectAndLeaveHalfAnd1<unsigned int, L::dimD>::Do(length_g());
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
      const unsigned int localLengthX = Domain<DomainType::LocalFourier,
                                               PartitionningType::Generic,
                                               MemoryLayout::Generic, 1>::length()[d::X];
      const unsigned int localVolume = Domain<DomainType::LocalFourier,
                                              PartitionningType::Generic,
                                              MemoryLayout::Generic, 1>::volume();

      const unsigned int indexLocal
        = Domain<DomainType::LocalFourier, PartitionningType::Generic,
                 MemoryLayout::Generic, 1>::getIndex(MathVector<unsigned int, 3>({iP[d::X] - iP[d::X]/localLengthX * localLengthX, iP[d::Y], iP[d::Z]}));
      return iP[d::X]/localLengthX * localVolume + indexLocal;
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      const unsigned int localLengthX = Domain<DomainType::LocalFourier,
        PartitionningType::Generic,
                                               MemoryLayout::Generic, 1>::length()[d::X];
      const unsigned int localVolume = Domain<DomainType::LocalFourier,
                                              PartitionningType::Generic,
                                              MemoryLayout::Generic, 1>::volume();

      const unsigned int indexLocal
        = Domain<DomainType::LocalFourier, PartitionningType::Generic,
                 MemoryLayout::Generic, 1>::getIndex(MathVector<unsigned int, 3>({iP[d::X] - iP[d::X]/localLengthX * localLengthX, iP[d::Y], iP[d::Z]}), iC);
      return iP[d::X]/localLengthX * localVolume + indexLocal;
    }

  };

  template <PartitionningType partitionningType, unsigned int NumberComponents>
  struct Domain<DomainType::GlobalFourier, partitionningType,
                MemoryLayout::Generic, NumberComponents>
    : public Domain<DomainType::GlobalFourier, partitionningType,
                    MemoryLayout::Generic, 1>
  {
  public:
    using Domain<DomainType::GlobalFourier, partitionningType,
                 MemoryLayout::Generic, 1>::start;
    using Domain<DomainType::GlobalFourier, partitionningType,
                 MemoryLayout::Generic, 1>::end;
    using Domain<DomainType::GlobalFourier, partitionningType,
                 MemoryLayout::Generic, 1>::length;
    using Domain<DomainType::GlobalFourier, partitionningType,
                 MemoryLayout::Generic, 1>::volume;
    using Domain<DomainType::GlobalFourier, partitionningType,
                 MemoryLayout::Generic, 1>::offset;

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP) {
      const unsigned int localLengthX = Domain<DomainType::LocalFourier,
                                               PartitionningType::Generic,
                                               MemoryLayout::Generic, NumberComponents>::length()[d::X];
      const unsigned int localVolume = Domain<DomainType::LocalFourier,
                                              PartitionningType::Generic,
                                              MemoryLayout::Generic, NumberComponents>::volume();

      const unsigned int indexLocal = Domain<DomainType::LocalFourier,
                                             PartitionningType::Generic,
                                             MemoryLayout::Generic, NumberComponents>::getIndex({iP[d::X] - iP[d::X]/localLengthX * localLengthX, iP[d::Y], iP[d::Z]});
      return iP[d::X]/localLengthX * NumberComponents * localVolume + indexLocal;
    }

    #pragma omp declare simd
    HOST DEVICE
    static inline unsigned int getIndex(const MathVector<unsigned int, 3>& iP,
                                        const unsigned int iC) {
      const unsigned int localLengthX = Domain<DomainType::LocalFourier,
                                               PartitionningType::Generic,
                                               MemoryLayout::Generic, NumberComponents>::length()[d::X];
      const unsigned int localVolume = Domain<DomainType::LocalFourier,
                                              PartitionningType::Generic,
                                              MemoryLayout::Generic, NumberComponents>::volume();

      const unsigned int indexLocal = Domain<DomainType::LocalFourier,
                                             PartitionningType::Generic,
                                             MemoryLayout::Generic, NumberComponents>::getIndex({iP[d::X] - iP[d::X]/localLengthX * localLengthX, iP[d::Y], iP[d::Z]}, iC);
      return iP[d::X]/localLengthX * NumberComponents * localVolume + indexLocal;
    }
  };


  typedef Domain<DomainType::GlobalFourier, partitionningT,
                 MemoryLayout::Generic, 1> gFD;
  typedef Domain<DomainType::GlobalFourier, partitionningT,
                 MemoryLayout::Generic, L::dimD> gFDD;
  typedef Domain<DomainType::GlobalFourier, partitionningT,
                 MemoryLayout::Generic, L::dimQ> gFQD;

  typedef Domain<DomainType::LocalFourier, PartitionningType::Generic,
                 MemoryLayout::Generic, 1> lFD;
  typedef Domain<DomainType::LocalFourier, PartitionningType::Generic,
    MemoryLayout::Generic, L::dimD> lFDD;

}

#endif // DOMAIN_H
