#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "StaticArray.h"
#include "Domain.h"

namespace lbm {

  template<class T>
  struct Packer {
  public:
    DEVICE HOST
    inline void operator()(const Position& iP,
                           T * const local, T * halo) {
      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        local[hSD::getIndexLocal(iP, iQ)] = halo[hSD::getIndex(iP, iQ)];
      }
    }
  };

  template<class T>
  struct Unpacker {
  public:
    DEVICE HOST
    inline void operator()(const Position& iP,
                           T * halo, T * const local) {
      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        halo[hSD::getIndex(iP, iQ)] = local[hSD::getIndexLocal(iP, iQ)];
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
    DEVICE HOST
    inline void applyX(const Position& iP,
                              T * f) {
      { INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyX",5) }

      Position iP_Origin = {L::halo()[d::X], iP[d::Y], iP[d::Z]};
      Position iP_Destination = {L::halo()[d::X] + lSD::sLength()[d::X],
                                                    iP[d::Y], iP[d::Z]};

      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
      }

      iP_Origin =  Position({L::halo()[d::X]+ lSD::sLength()[d::X] -1,
            iP[d::Y], iP[d::Z]});
      iP_Destination =  Position({0, iP[d::Y], iP[d::Z]});

        for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
          f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
        }
    }

    DEVICE HOST
    inline void applyY(const Position& iP,
                              T * f) {
      INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyY",5)

      Position iP_Origin = {iP[d::X], L::halo()[d::Y], iP[d::Z]};
      Position iP_Destination = {iP[d::X],
                                                    L::halo()[d::Y] + lSD::sLength()[d::Y],
                                                    iP[d::Z]};

      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
      }

      iP_Origin =  Position({iP[d::X],
            L::halo()[d::Y]+ lSD::sLength()[d::Y] -1, iP[d::Z]});
      iP_Destination =  Position({iP[d::X], 0, iP[d::Z]});

      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
      }
    }

    DEVICE HOST
    inline void applyZ(const Position& iP,
                       T * f) {
      { INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyZ",5) }

      Position iP_Origin = {iP[d::X], iP[d::Y], L::halo()[d::Z]};
      Position iP_Destination = {iP[d::X], iP[d::Y],
            L::halo()[d::Z] + lSD::sLength()[d::Z]};

      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
      }

      iP_Origin =  Position({iP[d::X], iP[d::Y],
            L::halo()[d::Z] + lSD::sLength()[d::Z] - 1});
      iP_Destination =  Position({iP[d::X], iP[d::Y], 0});

      for(auto iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
      }
    }
  };

  // template<class T, PartitionningType partitionningType>
  // class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
  //                partitionningType, 1> {};

  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::OneD, Implementation::MPI, 1>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 1> {
  public:
    DEVICE HOST
    inline void operator()(const Position& iP,
                             T * f) {
    }
  };

  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::OneD, Implementation::MPI, 2>
 : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 2> {
  private:
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                   PartitionningType::Generic, Implementation::Generic, 2>::applyY;

  public:
    DEVICE HOST
    inline void operator()(const Position& iP,
                           T * f) {
      applyY(iP, f);
    }
  };

  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::TwoD, Implementation::MPI, 2>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 2> {
  public:
    DEVICE HOST
    inline void operator()(const Position& iP,
                             T * f) {
    }
  };


  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::OneD, Implementation::MPI, 3>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 3> {
  private:
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                   PartitionningType::Generic, Implementation::Generic, 3>::applyY;
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                   PartitionningType::Generic, Implementation::Generic, 3>::applyZ;

  public:
    DEVICE HOST
    inline void operator()(const Position& iP,
                             T * f) {
      applyY(iP, f);
      applyZ(iP, f);
    }
  };

  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::TwoD, Implementation::MPI, 3>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 3> {
  private:
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                   PartitionningType::Generic, Implementation::Generic, 3>::applyZ;

  public:
    DEVICE HOST
    inline void operator()(const Position& iP,
                             T * f) {
      applyZ(iP, f);
    }
  };


  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::ThreeD, Implementation::MPI, 3>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 3> {
  public:
    DEVICE HOST
    inline void operator()(const Position& iP,
                             T * f) {
    }
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

#endif // BOUNDARY_H
