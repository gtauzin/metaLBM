#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "StaticArray.h"
#include "Domain.h"

namespace lbm {

  template<class T, BoundaryType boundaryType, AlgorithmType algorithmType,
    PartitionningType partitionningType, Implementation implementation,
    unsigned int Dimension>
  class Boundary {};

  template<class T, unsigned int Dimension>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::Generic, Implementation::Generic, Dimension> {
  public:
    DEVICE HOST
    inline void applyX(const MathVector<unsigned int, 3>& iP,
                              T * RESTRICT f) {
      INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyX",5)

        MathVector<unsigned int, 3> iP_Origin = {L::halo()[d::X], iP[d::Y], iP[d::Z]};
      MathVector<unsigned int, 3> iP_Destination = {L::halo()[d::X] + lSD::length()[d::X],
                                                    iP[d::Y], iP[d::Z]};

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
      }

      iP_Origin =  MathVector<unsigned int, 3>({L::halo()[d::X]+ lSD::length()[d::X] -1,
            iP[d::Y], iP[d::Z]});
      iP_Destination =  MathVector<unsigned int, 3>({0, iP[d::Y], iP[d::Z]});

        for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
          f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
        }
    }

    DEVICE HOST
    inline void applyY(const MathVector<unsigned int, 3>& iP,
                              T * RESTRICT f) {
      INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyY",5)

      MathVector<unsigned int, 3> iP_Origin = {iP[d::X], L::halo()[d::Y], iP[d::Z]};
      MathVector<unsigned int, 3> iP_Destination = {iP[d::X],
                                                    L::halo()[d::Y] + lSD::length()[d::Y],
                                                    iP[d::Z]};

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
      }

      iP_Origin =  MathVector<unsigned int, 3>({iP[d::X],
            L::halo()[d::Y]+ lSD::length()[d::Y] -1, iP[d::Z]});
      iP_Destination =  MathVector<unsigned int, 3>({iP[d::X], 0, iP[d::Z]});

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
      }
    }

    DEVICE HOST
    inline void applyZ(const MathVector<unsigned int, 3>& iP,
                              T * RESTRICT f) {
      INSTRUMENT_OFF("Boundary<T, boundaryType, algorithmType>::applyZ",5)

        MathVector<unsigned int, 3> iP_Origin = {iP[d::X], iP[d::Y], L::halo()[d::Z]};
      MathVector<unsigned int, 3>iP_Destination = {iP[d::X], iP[d::Y],
            L::halo()[d::Z] + lSD::length()[d::Z]};

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
      }

      iP_Origin =  MathVector<unsigned int, 3>({iP[d::X], iP[d::Y],
            L::halo()[d::Z] + lSD::length()[d::Z] - 1});
      iP_Destination =  MathVector<unsigned int, 3>({iP[d::X], iP[d::Y], 0});

      for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
        f[hSD::getIndex(iP_Destination, iQ)] = f[hSD::getIndex(iP_Origin, iQ)];
      }
    }
  };

  // template<class T, PartitionningType partitionningType>
  // class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
  //                partitionningType, 1> {};

  template<class T, PartitionningType partitionningType>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 partitionningType, Implementation::Serial, 1>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 1> {
  private:
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
      PartitionningType::Generic, Implementation::Generic, 1>::applyX;
  public:
    DEVICE HOST
    inline void operator()(const MathVector<unsigned int, 3>& iP,
                             T * RESTRICT f) {
      applyX(iP, f);
    }
  };

  template<class T>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 PartitionningType::OneD, Implementation::MPI, 1>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 1> {
  public:
    DEVICE HOST
    inline void operator()(const MathVector<unsigned int, 3>& iP,
                             T * RESTRICT f) {
    }
  };


  // template<class T, PartitionningType partitionningType>
  // class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
  //                partitionningType, 2>
  //   : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
  //                     PartitionningType::Generic, 0> {
  // };

  template<class T, PartitionningType partitionningType>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 partitionningType, Implementation::Serial, 2>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 2> {
  private:
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                   PartitionningType::Generic, Implementation::Generic, 2>::applyX;
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                   PartitionningType::Generic, Implementation::Generic, 2>::applyY;

  public:
    DEVICE HOST
    inline void operator()(const MathVector<unsigned int, 3>& iP,
                             T * RESTRICT f) {
      applyY(iP, f);
      applyX(iP, f);
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
    inline void operator()(const MathVector<unsigned int, 3>& iP,
                             T * RESTRICT f) {
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
    inline void operator()(const MathVector<unsigned int, 3>& iP,
                             T * RESTRICT f) {
    }
  };


  template<class T, PartitionningType partitionningType>
  class Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                 partitionningType, Implementation::Serial, 3>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                      PartitionningType::Generic, Implementation::Generic, 3> {
  private:
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                   PartitionningType::Generic, Implementation::Generic, 3>::applyX;
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                   PartitionningType::Generic, Implementation::Generic, 3>::applyY;
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull,
                   PartitionningType::Generic, Implementation::Generic, 3>::applyZ;

  public:
    DEVICE HOST
    inline void operator()(const MathVector<unsigned int, 3>& iP,
                             T * RESTRICT f) {
      applyX(iP, f);
      applyY(iP, f);
      applyZ(iP, f);
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
    inline void operator()(const MathVector<unsigned int, 3>& iP,
                             T * RESTRICT f) {
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
    inline void operator()(const MathVector<unsigned int, 3>& iP,
                             T * RESTRICT f) {
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
    inline void operator()(const MathVector<unsigned int, 3>& iP,
                             T * RESTRICT f) {
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
