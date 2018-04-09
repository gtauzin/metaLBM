#pragma once

#include "Commons.h"
#include "Options.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Domain.h"

namespace lbm {

  template<class T>
  struct Packer {
  public:
    DEVICE HOST
    inline void operator()(const Position& iP,
                           T * const local, T * halo,
                           const unsigned int numberElements) {
      auto indexLocal = hSD::getIndexLocal(iP);
      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        (local+iQ*numberElements)[indexLocal] = halo[hSD::getIndex(iP, iQ)];
      }
    }
  };

  template<class T>
  struct Unpacker {
  public:
    DEVICE HOST
    inline void operator()(const Position& iP,
                           T * halo, T * const local,
                           const unsigned int numberElements) {
      auto indexLocal = hSD::getIndexLocal(iP);
      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        halo[hSD::getIndex(iP, iQ)] = (local+iQ*numberElements)[indexLocal];
      }
    }
  };



  template<class T, BoundaryType boundaryType, AlgorithmType algorithmType,
           PartitionningType partitionningType, Implementation implementation,
           unsigned int Dimension>
  class Boundary {};

  template<class T, unsigned int Dimension>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::Generic, Implementation::Generic, Dimension> {
  public:
    DEVICE HOST INLINE
    static void applyYBottom(const Position& iP, T * haloDistributionPtr) {
      { INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyYBottom",5) }

      Position iP_Destination = {iP[d::X], iP[d::Y] + lSD::sLength()[d::Y], iP[d::Z]};

      #pragma unroll
      for(auto iQ = 0; iQ < L::faceQ; ++iQ) {
        haloDistributionPtr[hSD::getIndex(iP_Destination, L::iQ_Bottom()[iQ])]
          = haloDistributionPtr[hSD::getIndex(iP, L::iQ_Bottom()[iQ])];
      }
    }

    DEVICE HOST INLINE
    static void applyYTop(const Position& iP, T * haloDistributionPtr) {
      { INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyYTop",5) }

      Position iP_Destination{iP[d::X],
          iP[d::Y] - (L::halo()[d::Y]+lSD::sLength()[d::Y] - 1), iP[d::Z]};

      #pragma unroll
      for(auto iQ = 0; iQ < L::faceQ; ++iQ) {
        haloDistributionPtr[hSD::getIndex(iP_Destination, L::iQ_Top()[iQ])]
          = haloDistributionPtr[hSD::getIndex(iP, L::iQ_Top()[iQ])];
      }
    }

    DEVICE HOST INLINE
    static void applyZFront(const Position& iP, T * haloDistributionPtr) {
      { INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyZFront",5) }

      Position iP_Destination = {iP[d::X], iP[d::Y], iP[d::Z] + lSD::sLength()[d::Z]};

      #pragma unroll
      for(auto iQ = 0; iQ < L::faceQ; ++iQ) {
        haloDistributionPtr[hSD::getIndex(iP_Destination, L::iQ_Front()[iQ])]
          = haloDistributionPtr[hSD::getIndex(iP, L::iQ_Front()[iQ])];
      }
    }

    DEVICE HOST INLINE
    static void applyZBack(const Position& iP, T * haloDistributionPtr) {
      { INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyZBack",5) }

      Position iP_Destination{iP[d::X], iP[d::Y],
          iP[d::Z] - (L::halo()[d::Z]+lSD::sLength()[d::Z] - 1)};

      #pragma unroll
      for(auto iQ = 0; iQ < L::faceQ; ++iQ) {
        haloDistributionPtr[hSD::getIndex(iP_Destination, L::iQ_Back()[iQ])]
          = haloDistributionPtr[hSD::getIndex(iP, L::iQ_Back()[iQ])];
      }
    }
  };

  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::OneD, Implementation::MPI, 1>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 1> {
  public:
    DEVICE HOST INLINE
    static void applyYBottom(const Position& iP, T * haloDistributionPtr) {}

    DEVICE HOST INLINE
    static void applyYTop(const Position& iP, T * haloDistributionPtr) {}

    DEVICE HOST INLINE
    static void applyZFront(const Position& iP, T * haloDistributionPtr) {}

    DEVICE HOST INLINE
    static void applyZBack(const Position& iP, T * haloDistributionPtr) {}
  };

  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::OneD, Implementation::MPI, 2>
 : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 2> {
  private:
    using Base = Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                          PartitionningType::Generic, Implementation::Generic, 2>;

  public:
    using Base::applyYTop;
    using Base::applyYBottom;

    DEVICE HOST INLINE
    static void applyZFront(const Position& iP, T * haloDistributionPtr) {}

    DEVICE HOST INLINE
    static void applyZBack(const Position& iP, T * haloDistributionPtr) {}
  };

  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::TwoD, Implementation::MPI, 2>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 2> {
  public:
    DEVICE HOST INLINE
    static void applyYBottom(const Position& iP, T * haloDistributionPtr) {}

    DEVICE HOST INLINE
    static void applyYTop(const Position& iP, T * haloDistributionPtr) {}

    DEVICE HOST INLINE
    static void applyZFront(const Position& iP, T * haloDistributionPtr) {}

    DEVICE HOST INLINE
    static void applyZBack(const Position& iP, T * haloDistributionPtr) {}
  };


  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::OneD, Implementation::MPI, 3>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 3> {
  private:
    using Base = Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                   PartitionningType::Generic, Implementation::Generic, 3>;

  public:
    using Base::applyYTop;
    using Base::applyYBottom;
    using Base::applyZFront;
    using Base::applyZBack;

  };

  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::TwoD, Implementation::MPI, 3>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 3> {
  private:
    using Base = Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                          PartitionningType::Generic, Implementation::Generic, 3>;

  public:
    DEVICE HOST INLINE
    static void applyYBottom(const Position& iP, T * haloDistributionPtr) {}

    DEVICE HOST INLINE
    static void applyYTop(const Position& iP, T * haloDistributionPtr) {}

    using Base::applyZFront;
    using Base::applyZBack;
  };


  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::ThreeD, Implementation::MPI, 3>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 3> {
  public:
    DEVICE HOST INLINE
    static void applyYBottom(const Position& iP, T * haloDistributionPtr) {}

    DEVICE HOST INLINE
    static void applyYTop(const Position& iP, T * haloDistributionPtr) {}

    DEVICE HOST INLINE
    static void applyZFront(const Position& iP, T * haloDistributionPtr) {}

    DEVICE HOST INLINE
    static void applyZBack(const Position& iP, T * haloDistributionPtr) {}
  };


  template<class T, AlgorithmType algorithmType, PartitionningType partitionningType,
           Implementation implementation, unsigned int Dimension>
  class Boundary<T, BoundaryType::BounceBack_Halfway, algorithmType,
                 partitionningType, implementation, Dimension>
    : public Boundary<T, BoundaryType::Generic, algorithmType,
                      PartitionningType::Generic, Implementation::Generic, Dimension> {
  protected:

  public:
    void apply() {}

  };


  template<class T, AlgorithmType algorithmType, PartitionningType partitionningType,
           Implementation implementation, unsigned int Dimension>
  class Boundary<T, BoundaryType::Entropic, algorithmType,
                 partitionningType, implementation, Dimension>
    : public Boundary<T, BoundaryType::Generic, algorithmType,
                      PartitionningType::Generic, Implementation::Generic, Dimension> {
  protected:

  public:
    void apply() {}

  };

}
