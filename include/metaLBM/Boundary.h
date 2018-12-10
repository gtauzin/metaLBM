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
    LBM_DEVICE LBM_HOST inline
    void operator()(const Position& iP, T* const local, T* halo,
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
    LBM_DEVICE LBM_HOST inline
    void operator()(const Position& iP, T* halo, T* const local,
                    const unsigned int numberElements) {
      auto indexLocal = hSD::getIndexLocal(iP);
      for (auto iQ = 0; iQ < L::dimQ; ++iQ) {
        halo[hSD::getIndex(iP, iQ)]
          = (local + iQ * numberElements)[indexLocal];
      }
    }
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType,
            unsigned int Dimension>
  class Boundary {};

  template <class T, unsigned int Dimension>
  class Boundary<T, BoundaryType::Periodic_OUT, AlgorithmType::Pull, Dimension> {
  public:
    LBM_DEVICE LBM_HOST LBM_INLINE
    static void applyYBottom(const Position& iP, T* haloDistributionDestinationPtr,
                           T* haloDistributionSourcePtr) {
      LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyYTop", 5)

      Position iP_Source{iP[d::X], iP[d::Y] + lSD::sLength()[d::Y], iP[d::Z]};

      #pragma unroll
      for(auto iQ = 0; iQ < L::faceQ; ++iQ) {
          haloDistributionDestinationPtr[hSD::getIndex(iP, L::iQ_Top()[iQ])] =
            haloDistributionDestinationPtr[hSD::getIndex(iP_Source, L::iQ_Top()[iQ])];
      }
    }

    LBM_DEVICE LBM_HOST LBM_INLINE
    static void applyYTop(const Position& iP, T* haloDistributionDestinationPtr,
                          T* haloDistributionSourcePtr) {
      LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyYBottom", 5)

        Position iP_Source{iP[d::X], iP[d::Y] - lSD::sLength()[d::Y], iP[d::Z]};

      #pragma unroll
      for(auto iQ = 0; iQ < L::faceQ; ++iQ) {
        haloDistributionDestinationPtr[hSD::getIndex(iP, L::iQ_Bottom()[iQ])] =
          haloDistributionDestinationPtr[hSD::getIndex(iP_Source, L::iQ_Bottom()[iQ])];
      }
    }

    LBM_DEVICE LBM_HOST LBM_INLINE
    static void applyZFront(const Position& iP, T* haloDistributionDestinationPtr,
                            T* haloDistributionSourcePtr) {
      LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyZBack", 5)

      Position iP_Source{iP[d::X], iP[d::Y], iP[d::Z] + lSD::sLength()[d::Z]};

      #pragma unroll
      for(auto iQ = 0; iQ < L::faceQ; ++iQ) {
        haloDistributionDestinationPtr[hSD::getIndex(iP, L::iQ_Back()[iQ])] =
          haloDistributionDestinationPtr[hSD::getIndex(iP_Source, L::iQ_Back()[iQ])];
      }
    }

    LBM_DEVICE LBM_HOST LBM_INLINE
    static void applyZBack(const Position& iP, T* haloDistributionDestinationPtr,
                           T* haloDistributionSourcePtr) {
      LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyZFront", 5)

      Position iP_Source{iP[d::X], iP[d::Y], iP[d::Z] - lSD::sLength()[d::Z]};

      #pragma unroll
      for(auto iQ = 0; iQ < L::faceQ; ++iQ) {
        haloDistributionDestinationPtr[hSD::getIndex(iP, L::iQ_Front()[iQ])] =
          haloDistributionDestinationPtr[hSD::getIndex(iP_Source, L::iQ_Front()[iQ])];
      }
    }

  };

  template <class T, unsigned int Dimension>
  class Boundary<T, BoundaryType::Periodic_IN, AlgorithmType::Pull, Dimension> {
  public:
    LBM_DEVICE LBM_HOST LBM_INLINE
    static void applyXLeft(const Position& iP, T* haloDistributionDestinationPtr,
                           T* haloDistributionSourcePtr) {
      LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyXLeft", 5)

      Position iP_Destination, iP_Source;
      unsigned int cumulativeLevel = 0;

      for(auto iH = 0; iH < 2 * L::halo()[d::X] - iP[d::X]; ++iH) {
        cumulativeLevel += L::level()[L::dimH-1-iH];

        #pragma unroll
        for(auto iQ = 1 + 2*L::faceQ - cumulativeLevel;
            iQ < 1 + 2*L::faceQ - (cumulativeLevel - L::level()[L::dimH-1-iH]); ++iQ) {
          iP_Destination = iP - uiL::celerity()[iQ];

          iP_Source =
            Position{{iP_Destination[d::X] < L::halo()[d::X]
                      ? iP_Destination[d::X] + lSD::sLength()[d::X] : iP_Destination[d::X],
                      iP_Destination[d::Y] < L::halo()[d::Y]
                      ? iP_Destination[d::Y] + lSD::sLength()[d::Y]
                      : (iP_Destination[d::Y] > L::halo()[d::Y] + lSD::sLength()[d::Y] - 1
                         ? iP_Destination[d::Y] - lSD::sLength()[d::Y] : iP_Destination[d::Y]),
                      iP_Destination[d::Z] < L::halo()[d::Z]
                      ? iP_Destination[d::Z] + lSD::sLength()[d::Z]
                      : (iP_Destination[d::Z] > L::halo()[d::Z] + lSD::sLength()[d::Z] - 1
                         ? iP_Destination[d::Z] - lSD::sLength()[d::Z] : iP_Destination[d::Z])}};

          haloDistributionDestinationPtr[hSD::getIndex(iP_Destination, iQ)] =
            haloDistributionSourcePtr[hSD::getIndex(iP_Source, iQ)];
        }
      }
    }

    LBM_DEVICE LBM_HOST LBM_INLINE
    static void applyXRight(const Position& iP, T* haloDistributionDestinationPtr,
                            T* haloDistributionSourcePtr) {
      LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyXRight", 5)

      Position iP_Destination, iP_Source;
      unsigned int cumulativeLevel = 0;

      for(auto iH = 0; iH < iP[d::X] - lSD::sLength()[d::X] + 1; ++iH) {
        cumulativeLevel += L::level()[L::dimH-1-iH];

        #pragma unroll
        for(auto iQ = 1 + L::faceQ - cumulativeLevel;
            iQ < 1 + L::faceQ - (cumulativeLevel - L::level()[L::dimH-1-iH]); ++iQ) {
          iP_Destination = iP - uiL::celerity()[iQ];

          iP_Source =
            Position{{iP_Destination[d::X] > L::halo()[d::X] + lSD::sLength()[d::X] - 1
                      ? iP_Destination[d::X] - lSD::sLength()[d::X] : iP_Destination[d::X],
                      iP_Destination[d::Y] < L::halo()[d::Y]
                      ? iP_Destination[d::Y] + lSD::sLength()[d::Y]
                      : (iP_Destination[d::Y] > L::halo()[d::Y] + lSD::sLength()[d::Y] - 1
                         ? iP_Destination[d::Y] - lSD::sLength()[d::Y] : iP_Destination[d::Y]),
                      iP_Destination[d::Z] < L::halo()[d::Z]
                      ? iP_Destination[d::Z] + lSD::sLength()[d::Z]
                      : (iP_Destination[d::Z] > L::halo()[d::Z] + lSD::sLength()[d::Z] - 1
                         ? iP_Destination[d::Z] - lSD::sLength()[d::Z] : iP_Destination[d::Z])}};

          haloDistributionDestinationPtr[hSD::getIndex(iP_Destination, iQ)] =
            haloDistributionSourcePtr[hSD::getIndex(iP_Source, iQ)];

        }
      }
    }


    LBM_DEVICE LBM_HOST LBM_INLINE
    static void applyYBottom(const Position& iP, T* haloDistributionDestinationPtr,
                             T* haloDistributionSourcePtr) {
      LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyYBottom", 5)

      Position iP_Destination, iP_Source;
      unsigned int cumulativeLevel = 0;

      for(auto iH = 0; iH < 2 * L::halo()[d::Y] - iP[d::Y]; ++iH) {
        cumulativeLevel += L::level()[L::dimH-1-iH];

        #pragma unroll
        for(auto iQ = L::faceQ - cumulativeLevel;
            iQ < L::faceQ - (cumulativeLevel - L::level()[L::dimH-1-iH]); ++iQ) {
          iP_Destination = iP - uiL::celerity()[L::iQ_Top()[iQ]];

          iP_Source =
            Position{{iP_Destination[d::X],
                      iP_Destination[d::Y] < L::halo()[d::Y]
                      ? iP_Destination[d::Y] + lSD::sLength()[d::Y]
                      : iP_Destination[d::Y],
                      iP_Destination[d::Z] < L::halo()[d::Z]
                      ? iP_Destination[d::Z] + lSD::sLength()[d::Z]
                      : (iP_Destination[d::Z] > L::halo()[d::Z] + lSD::sLength()[d::Z] - 1
                         ? iP_Destination[d::Z] - lSD::sLength()[d::Z]
                         : iP_Destination[d::Z])}};

          haloDistributionDestinationPtr[hSD::getIndex(iP_Destination, L::iQ_Top()[iQ])] =
            haloDistributionSourcePtr[hSD::getIndex(iP_Source, L::iQ_Top()[iQ])];
        }
      }
    }

    LBM_DEVICE LBM_HOST LBM_INLINE
    static void applyYTop(const Position& iP, T* haloDistributionDestinationPtr,
                          T* haloDistributionSourcePtr) {
      LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyYTop", 5)

      Position iP_Destination, iP_Source;
      unsigned int cumulativeLevel = 0;

      //printf("Top iP: %d %d\n", iP[d::X], iP[d::Y]);

      for(auto iH = 0; iH < iP[d::Y] - lSD::sLength()[d::Y] + 1; ++iH) {
        cumulativeLevel += L::level()[L::dimH-1-iH];

        #pragma unroll
        for(auto iQ = L::faceQ - cumulativeLevel;
            iQ < L::faceQ - (cumulativeLevel - L::level()[L::dimH-1-iH]); ++iQ) {
          iP_Destination = iP - uiL::celerity()[L::iQ_Bottom()[iQ]];

          iP_Source =
            Position{{iP_Destination[d::X],
                      iP_Destination[d::Y] > L::halo()[d::Y] + lSD::sLength()[d::Y] - 1
                      ? iP_Destination[d::Y] - lSD::sLength()[d::Y]
                      : iP_Destination[d::Y],
                      iP_Destination[d::Z] < L::halo()[d::Z]
                      ? iP_Destination[d::Z] + lSD::sLength()[d::Z]
                      : (iP_Destination[d::Z] > L::halo()[d::Z] + lSD::sLength()[d::Z] - 1
                         ? iP_Destination[d::Z] - lSD::sLength()[d::Z]
                         : iP_Destination[d::Z])}};

          haloDistributionDestinationPtr[hSD::getIndex(iP_Destination, L::iQ_Bottom()[iQ])] =
            haloDistributionSourcePtr[hSD::getIndex(iP_Source, L::iQ_Bottom()[iQ])];
        }
      }
    }

    LBM_DEVICE LBM_HOST LBM_INLINE
    static void applyZFront(const Position& iP, T* haloDistributionDestinationPtr,
                            T* haloDistributionSourcePtr) {
      LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyZBack", 5)

      Position iP_Destination, iP_Source;
      unsigned int cumulativeLevel = 0;

      for(auto iH = 0; iH < 2 * L::halo()[d::Z] - iP[d::Z]; ++iH) {
        cumulativeLevel += L::level()[L::dimH-1-iH];

        #pragma unroll
        for(auto iQ = L::faceQ - cumulativeLevel;
            iQ < L::faceQ - (cumulativeLevel - L::level()[L::dimH-1-iH]); ++iQ) {
          iP_Destination = iP - uiL::celerity()[L::iQ_Back()[iQ]];

          iP_Source =
            Position{{iP_Destination[d::X], iP_Destination[d::Y],
                      iP_Destination[d::Z] < L::halo()[d::Z]
                      ? iP_Destination[d::Z] + lSD::sLength()[d::Z]
                      : (iP_Destination[d::Z] > L::halo()[d::Z] + lSD::sLength()[d::Z] - 1
                         ? iP_Destination[d::Z] - lSD::sLength()[d::Z]
                         : iP_Destination[d::Z])}};

          haloDistributionDestinationPtr[hSD::getIndex(iP_Destination, L::iQ_Back()[iQ])] =
            haloDistributionSourcePtr[hSD::getIndex(iP_Source, L::iQ_Back()[iQ])];
        }
      }
    }

    LBM_DEVICE LBM_HOST LBM_INLINE
    static void applyZBack(const Position& iP, T* haloDistributionDestinationPtr,
                           T* haloDistributionSourcePtr) {
      LBM_INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyZFront", 5)

      Position iP_Destination, iP_Source;
      unsigned int cumulativeLevel = 0;

      for(auto iH = 0; iH < iP[d::Z] - lSD::sLength()[d::Z] + 1; ++iH) {
        cumulativeLevel += L::level()[L::dimH-1-iH];

        #pragma unroll
        for(auto iQ = L::faceQ - cumulativeLevel;
            iQ < L::faceQ - (cumulativeLevel - L::level()[L::dimH-1-iH]); ++iQ) {
          iP_Destination = iP - uiL::celerity()[L::iQ_Front()[iQ]];

          iP_Source =
            Position{{iP_Destination[d::X], iP_Destination[d::Y],
                      iP_Destination[d::Z] < L::halo()[d::Z]
                      ? iP_Destination[d::Z] + lSD::sLength()[d::Z]
                      : (iP_Destination[d::Z] > L::halo()[d::Z] + lSD::sLength()[d::Z] - 1
                         ? iP_Destination[d::Z] - lSD::sLength()[d::Z]
                         : iP_Destination[d::Z])}};

          haloDistributionDestinationPtr[hSD::getIndex(iP_Destination, L::iQ_Front()[iQ])] =
            haloDistributionSourcePtr[hSD::getIndex(iP_Source, L::iQ_Front()[iQ])];
        }
      }
    }
  };


  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType, unsigned int Dimension>
  class LeftBoundary : public Boundary<T, boundaryType, algorithmType,Dimension> {
  private:
    using Base = Boundary<T, boundaryType, algorithmType, Dimension>;

  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {
      Base::applyXLeft(iP, haloDistributionDestinationPtr, haloDistributionSourcePtr);
    }
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType, unsigned int Dimension>
  class RightBoundary : public Boundary<T, boundaryType, algorithmType, Dimension> {
  private:
    using Base = Boundary<T, boundaryType, algorithmType, Dimension>;

  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {
      Base::applyXRight(iP, haloDistributionDestinationPtr, haloDistributionSourcePtr);
    }
  };

template <class T, BoundaryType boundaryType, AlgorithmType algorithmType, unsigned int Dimension>
  class BottomBoundary : public Boundary<T, boundaryType, algorithmType,Dimension> {
  private:
    using Base = Boundary<T, boundaryType, algorithmType, Dimension>;

  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {
      Base::applyYBottom(iP, haloDistributionDestinationPtr, haloDistributionSourcePtr);
    }
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType, unsigned int Dimension>
  class TopBoundary : public Boundary<T, boundaryType, algorithmType, Dimension> {
  private:
    using Base = Boundary<T, boundaryType, algorithmType, Dimension>;

  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {
      Base::applyYTop(iP, haloDistributionDestinationPtr, haloDistributionSourcePtr);
    }
  };

template <class T, BoundaryType boundaryType, AlgorithmType algorithmType>
  class BottomBoundary<T, boundaryType, algorithmType, 1> {
  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {}
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType>
  class TopBoundary<T, boundaryType, algorithmType, 1> {
  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {}
  };


  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType, unsigned int Dimension>
  class FrontBoundary : public Boundary<T, boundaryType, algorithmType, Dimension> {
  private:
    using Base = Boundary<T, boundaryType, algorithmType, Dimension>;

  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {
      Base::applyZFront(iP, haloDistributionDestinationPtr, haloDistributionSourcePtr);
    }
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType, unsigned int Dimension>
  class BackBoundary : public Boundary<T, boundaryType, algorithmType, Dimension> {
  private:
    using Base = Boundary<T, boundaryType, algorithmType, Dimension>;

  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {
      Base::applyZBack(iP, haloDistributionDestinationPtr, haloDistributionSourcePtr);
    }
  };

template <class T, BoundaryType boundaryType, AlgorithmType algorithmType>
  class FrontBoundary<T, boundaryType, algorithmType, 1> {
  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {}
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType>
  class BackBoundary<T, boundaryType, algorithmType, 1> {
  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {}
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType>
  class FrontBoundary<T, boundaryType, algorithmType, 2> {
  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {}
  };

  template <class T, BoundaryType boundaryType, AlgorithmType algorithmType>
  class BackBoundary<T, boundaryType, algorithmType, 2> {
  public:
    LBM_HOST LBM_DEVICE void operator()(const Position& iP, T* haloDistributionDestinationPtr,
                                        T* haloDistributionSourcePtr) {}
  };


  template <class T, AlgorithmType algorithmType, unsigned int Dimension>
  class Boundary<T, BoundaryType::BounceBack_Halfway, algorithmType,
                 Dimension>
    : public Boundary<T, BoundaryType::Generic, algorithmType, Dimension> {
  public:
    void apply() {}
  };

  template <class T, AlgorithmType algorithmType, unsigned int Dimension>
  class Boundary<T, BoundaryType::Entropic, algorithmType, Dimension>
    : public Boundary<T, BoundaryType::Generic, algorithmType, Dimension> {
  public:
    void apply() {}
  };

}  // namespace lbm
