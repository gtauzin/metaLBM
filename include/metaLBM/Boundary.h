#pragma once

#include "Commons.h"
#include "Domain.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Options.h"

namespace lbm {

  template <class T>
  struct Packer {
  public:
    LBM_DEVICE LBM_HOST inline void operator()(const Position& iP, T*
                                               const local, T* halo,
                                               const unsigned int numberElements) {
    auto indexLocal = hSD::getIndexLocal(iP);
    for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
      (local + iQ * numberElements)[indexLocal]
        = halo[hSD::getIndex(iP, iQ)];
    }
    }
  };

  template <class T>
  struct Unpacker {
  public:
    LBM_DEVICE LBM_HOST inline void operator()(const Position& iP,
                                               T* halo, T* const local,
                                               const unsigned int numberElements) {
      auto indexLocal = hSD::getIndexLocal(iP);
      for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
        halo[hSD::getIndex(iP, iQ)]
          = (local + iQ * numberElements)[indexLocal];
      }
    }
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType,
            PartitionningType partitionningType, Implementation implementation,
            unsigned int Dimension>
  class Boundary {};

  template <class T, unsigned int Dimension>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::Generic, Implementation::Generic,
                 Dimension> {
 public:
  LBM_DEVICE LBM_HOST LBM_INLINE static void applyYBottom(const Position& iP,
                                                          T* haloDistributionPtr) {
    LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyYBottom", 5)

    Position iP_Destination = {iP[d::X], iP[d::Y] + lSD::sLength()[d::Y],
                               iP[d::Z]};

    #pragma unroll
    for (auto iQ = 0; iQ < L::faceQ; ++iQ) {
      haloDistributionPtr[hSD::getIndex(iP_Destination, L::iQ_Bottom()[iQ])] =
          haloDistributionPtr[hSD::getIndex(iP, L::iQ_Bottom()[iQ])];
    }
  }

  LBM_DEVICE LBM_HOST LBM_INLINE static void applyYTop(const Position& iP,
                                                       T* haloDistributionPtr) {
    LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyYTop", 5)

    Position iP_Destination{
        iP[d::X], iP[d::Y] - (L::halo()[d::Y] + lSD::sLength()[d::Y] - 1),
        iP[d::Z]};

    #pragma unroll
    for (auto iQ = 0; iQ < L::faceQ; ++iQ) {
      haloDistributionPtr[hSD::getIndex(iP_Destination, L::iQ_Top()[iQ])] =
          haloDistributionPtr[hSD::getIndex(iP, L::iQ_Top()[iQ])];
    }
  }

  LBM_DEVICE LBM_HOST LBM_INLINE static void applyZFront(const Position& iP,
                                                         T* haloDistributionPtr) {
    LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyZFront", 5)

    Position iP_Destination = {iP[d::X], iP[d::Y],
                               iP[d::Z] + lSD::sLength()[d::Z]};

    #pragma unroll
    for (auto iQ = 0; iQ < L::faceQ; ++iQ) {
      haloDistributionPtr[hSD::getIndex(iP_Destination, L::iQ_Front()[iQ])] =
          haloDistributionPtr[hSD::getIndex(iP, L::iQ_Front()[iQ])];
    }
  }

  LBM_DEVICE LBM_HOST LBM_INLINE static void applyZBack(const Position& iP,
                                                        T* haloDistributionPtr) {
    LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyZBack", 5)

    Position iP_Destination{
        iP[d::X], iP[d::Y],
        iP[d::Z] - (L::halo()[d::Z] + lSD::sLength()[d::Z] - 1)};

    #pragma unroll
    for (auto iQ = 0; iQ < L::faceQ; ++iQ) {
      haloDistributionPtr[hSD::getIndex(iP_Destination, L::iQ_Back()[iQ])] =
          haloDistributionPtr[hSD::getIndex(iP, L::iQ_Back()[iQ])];
    }
  }
};

  template <class T, Implementation implementation>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
               PartitionningType::OneD, implementation, 1>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 1> {
  public:
    LBM_DEVICE LBM_HOST LBM_INLINE static void applyYBottom(const Position& iP,
                                                            T* haloDistributionPtr) {}

    LBM_DEVICE LBM_HOST LBM_INLINE static void applyYTop(const Position& iP,
                                                         T* haloDistributionPtr) {}

    LBM_DEVICE LBM_HOST LBM_INLINE static void applyZFront(const Position& iP,
                                                           T* haloDistributionPtr) {}

    LBM_DEVICE LBM_HOST LBM_INLINE static void applyZBack(const Position& iP,
                                                          T* haloDistributionPtr) {}
  };

  template <class T, Implementation implementation>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
               PartitionningType::OneD, implementation, 2>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 2> {
 private:
  using Base = Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                        PartitionningType::Generic, Implementation::Generic, 2>;

 public:
  using Base::applyYBottom;
  using Base::applyYTop;

    LBM_DEVICE LBM_HOST LBM_INLINE static void applyZFront(const Position& iP,
                                                           T* haloDistributionPtr) {}

    LBM_DEVICE LBM_HOST LBM_INLINE static void applyZBack(const Position& iP,
                                                          T* haloDistributionPtr) {}
  };

  template <class T, Implementation implementation>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::TwoD, implementation, 2>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 2> {
  public:
    LBM_DEVICE LBM_HOST LBM_INLINE static void applyYBottom(const Position& iP,
                                                            T* haloDistributionPtr) {}

    LBM_DEVICE LBM_HOST LBM_INLINE static void applyYTop(const Position& iP,
                                                         T* haloDistributionPtr) {
    }

    LBM_DEVICE LBM_HOST LBM_INLINE static void applyZFront(const Position& iP,
                                                           T* haloDistributionPtr) {}

    LBM_DEVICE LBM_HOST LBM_INLINE static void applyZBack(const Position& iP,
                                                          T* haloDistributionPtr) {}
  };

  template <class T, Implementation implementation>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::OneD, implementation, 3>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 3> {
 private:
    using Base = Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                          PartitionningType::Generic, Implementation::Generic, 3>;

 public:
  using Base::applyYBottom;
  using Base::applyYTop;
  using Base::applyZBack;
  using Base::applyZFront;
};

 template <class T, Implementation implementation>
 class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                PartitionningType::TwoD, implementation, 3>
   : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                     PartitionningType::Generic, Implementation::Generic, 3> {
 private:
   using Base = Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                         PartitionningType::Generic, Implementation::Generic, 3>;

 public:
  LBM_DEVICE LBM_HOST LBM_INLINE static void applyYBottom(const Position& iP,
                                                          T* haloDistributionPtr) {}

  LBM_DEVICE LBM_HOST LBM_INLINE static void applyYTop(const Position& iP,
                                                       T* haloDistributionPtr) {
  }

  using Base::applyZBack;
  using Base::applyZFront;
};

 template <class T, Implementation implementation>
 class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                PartitionningType::ThreeD, implementation, 3>
   : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                     PartitionningType::Generic, Implementation::Generic, 3> {
 public:
   LBM_DEVICE LBM_HOST LBM_INLINE static void applyYBottom(const Position& iP,
                                                           T* haloDistributionPtr) {}

   LBM_DEVICE LBM_HOST LBM_INLINE static void applyYTop(const Position& iP,
                                                        T* haloDistributionPtr) {}

   LBM_DEVICE LBM_HOST LBM_INLINE static void applyZFront(const Position& iP,
                                                          T* haloDistributionPtr) {}

   LBM_DEVICE LBM_HOST LBM_INLINE static void applyZBack(const Position& iP,
                                                         T* haloDistributionPtr) {}
 };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType,
            PartitionningType partitionningType, Implementation implementation,
            unsigned int Dimension>
  class BottomBoundary : public Boundary<T, boundaryType, algorithmType, partitionningType,
                                         implementation, Dimension> {
  private:
    using Base = Boundary<T, boundaryType, algorithmType, partitionningType,
                          implementation, Dimension>;

  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP,
                                        T* haloDistributionPtr) {
      Base::applyYBottom(iP, haloDistributionPtr);
    }
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType,
            PartitionningType partitionningType, Implementation implementation,
            unsigned int Dimension>
class TopBoundary : public Boundary<T, boundaryType, algorithmType, partitionningType,
                                    implementation, Dimension> {
 private:
  using Base = Boundary<T, boundaryType, algorithmType, partitionningType,
                        implementation, Dimension>;

  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP,
                                        T* haloDistributionPtr) {
      Base::applyYTop(iP, haloDistributionPtr);
    }
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType,
            PartitionningType partitionningType, Implementation implementation,
            unsigned int Dimension>
  class FrontBoundary : public Boundary<T, boundaryType, algorithmType, partitionningType,
                                        implementation, Dimension> {
 private:
  using Base = Boundary<T, boundaryType, algorithmType, partitionningType,
                        implementation, Dimension>;

  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP,
                                        T* haloDistributionPtr) {
      Base::applyZFront(iP, haloDistributionPtr);
    }
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType,
            PartitionningType partitionningType, Implementation implementation,
            unsigned int Dimension>
class BackBoundary : public Boundary<T, boundaryType, algorithmType, partitionningType,
                                     implementation, Dimension> {
 private:
    using Base = Boundary<T, boundaryType, algorithmType, partitionningType,
                        implementation, Dimension>;

  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP,
                                        T* haloDistributionPtr) {
      Base::applyZBack(iP, haloDistributionPtr);
    }
  };

  template <class T, AlgorithmType algorithmType, PartitionningType partitionningType,
            Implementation implementation, unsigned int Dimension>
  class Boundary<T, BoundaryType::BounceBack_Halfway, algorithmType,
                 partitionningType, implementation, Dimension>
    : public Boundary<T, BoundaryType::Generic, algorithmType, PartitionningType::Generic,
                      Implementation::Generic, Dimension> {
  public:
    void apply() {}
  };

  template <class T, AlgorithmType algorithmType, PartitionningType partitionningType,
            Implementation implementation, unsigned int Dimension>
  class Boundary<T, BoundaryType::Entropic, algorithmType, partitionningType,
                 implementation, Dimension>
    : public Boundary<T, BoundaryType::Generic, algorithmType, PartitionningType::Generic,
                      Implementation::Generic, Dimension> {
  public:
    void apply() {}
  };

}  // namespace lbm
