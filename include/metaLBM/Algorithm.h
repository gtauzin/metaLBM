#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <chrono>

#include "Commons.h"
#include "Options.h"
#include "Domain.h"
#include "Lattice.h"
#include "Field.h"
#include "Distribution.h"
#include "Moment.h"
#include "Collision.h"
#include "Boundary.h"
#include "Communication.h"
#include "Computation.h"

namespace lbm {

  template<class T, AlgorithmType, Architecture architecture>
  class Algorithm {};

  template<class T, Architecture architecture>
    class Algorithm<T, AlgorithmType::Generic, architecture> {
  protected:
    Communication_& communication;

    Field<T, 1, architecture, writeDensity>& densityField;
    Field<T, L::dimD, architecture, writeVelocity>& velocityField;
    Field<T, L::dimD, architecture, writeDensity>& forceField;
    Field<T, 1, architecture, writeDensity>& alphaField;
    Distribution<T, architecture>& f_Previous;
    Distribution<T, architecture>& f_Next;

    Collision_ collision;
    Moment<T> moment;
    Boundary_ boundary;

    std::chrono::duration<double> dtComputation;
    std::chrono::duration<double> dtCommunication;
    std::chrono::duration<double> dtTotal;

    Algorithm(Communication_& communication_in,
              Field<T, 1, architecture, writeDensity>& densityField_in,
              Field<T, L::dimD, architecture, writeVelocity>& velocityField_in,
              Field<T, L::dimD, architecture, writeDensity>& forceField_in,
              Field<T, 1, architecture, writeDensity>& alphaField_in,
              Distribution<T, architecture>& f_Previous_in,
              Distribution<T, architecture>& f_Next_in)
      : communication(communication_in)
      , densityField(densityField_in)
      , velocityField(velocityField_in)
      , forceField(forceField_in)
      , alphaField(alphaField_in)
      , f_Previous(f_Previous_in)
      , f_Next(f_Next_in)
      , collision(relaxationTime, forceAmplitude, forceWaveLength)
      , moment()
      , boundary()
    {}

    void storeLocalFields(const MathVector<unsigned int, 3>& iP,
                          const unsigned int iteration) {
      SCOREP_INSTRUMENT_OFF("Algorithm<T, AlgorithmType::Pull>::storeLocalFields")

      if(iteration%writeStep == 0) {
        const unsigned int indexLocal = hD::getIndexLocal(iP);

        densityField.setLocalValue(indexLocal, moment.getDensity());
        velocityField.setLocalVector(indexLocal, collision.getHydrodynamicVelocity());
        alphaField.setLocalValue(indexLocal, collision.getAlpha());
        forceField.setLocalVector(indexLocal, collision.getForce());
      }
    }

  public:
    double getCommunicationTime() {
      return dtCommunication.count();
    }

    double getComputationTime() {
      return dtComputation.count();
    }

    double getTotalTime() {
      return dtTotal.count();
    }
  };

  template<class T, Architecture architecture>
    class Algorithm<T, AlgorithmType::Pull, architecture>
    : public Algorithm<T, AlgorithmType::Generic, architecture> {
  private:
    using Algorithm<T, AlgorithmType::Generic, architecture>::communication;

    using Algorithm<T, AlgorithmType::Generic, architecture>::densityField;
    using Algorithm<T, AlgorithmType::Generic, architecture>::velocityField;
    using Algorithm<T, AlgorithmType::Generic, architecture>::forceField;
    using Algorithm<T, AlgorithmType::Generic, architecture>::alphaField;
    using Algorithm<T, AlgorithmType::Generic, architecture>::f_Previous;
    using Algorithm<T, AlgorithmType::Generic, architecture>::f_Next;

    using Algorithm<T, AlgorithmType::Generic, architecture>::collision;
    using Algorithm<T, AlgorithmType::Generic, architecture>::moment;
    using Algorithm<T, AlgorithmType::Generic, architecture>::boundary;

    using Algorithm<T, AlgorithmType::Generic, architecture>::dtComputation;
    using Algorithm<T, AlgorithmType::Generic, architecture>::dtCommunication;
    using Algorithm<T, AlgorithmType::Generic, architecture>::dtTotal;

    using Algorithm<T, AlgorithmType::Generic, architecture>::storeLocalFields;

  public:
    Algorithm(Communication_& communication_in,
              Field<T, 1, architecture, writeDensity>& densityField_in,
              Field<T, L::dimD, architecture, writeVelocity>& velocityField_in,
              Field<T, L::dimD, architecture, writeDensity>& forceField_in,
              Field<T, 1, architecture, writeDensity>& alphaField_in,
              Distribution<T, architecture>& f_Previous_in,
              Distribution<T, architecture>& f_Next_in)
      : Algorithm<T, AlgorithmType::Generic, architecture>(communication_in,
                                                           densityField_in, velocityField_in,
                                                           forceField_in, alphaField_in,
                                                           f_Previous_in, f_Next_in)
    {}

    void iterate(const unsigned int iteration) {
      SCOREP_INSTRUMENT_ON("Algorithm<T, AlgorithmType::Pull>::iterate")

      f_Previous.swapHalo(f_Next);

      //force.update(iteration);

      auto t0 = std::chrono::high_resolution_clock::now();
      communication.periodic(f_Previous.haloDeviceArray(),
                             f_Previous.haloHostArray(),
                             f_Previous.haloComputedData());

      //boundary.apply(f_Previous.haloData());

      auto t1 = std::chrono::high_resolution_clock::now();

      Computation_::Do(lD::start()+L::halo(), lD::end()+L::halo(),
                      [&] HOST DEVICE (MathVector<unsigned int, 3>& iP) {
          SCOREP_INSTRUMENT_ON("Algorithn<T, AlgorithmType::Pull>::lambda[fused_collide_and_push]")

          moment.calculateMoments(f_Previous.haloComputedData(), iP);

          collision.setForce(iP+gD::offset(communication.getRankMPI()));
          collision.setVariables(f_Previous.haloComputedData(), iP,
                                 moment.getDensity(), moment.getVelocity());

          UnrolledFor<0, L::dimQ>::Do([&] HOST DEVICE (unsigned int iQ) {
              f_Next.setHaloField(hD::getIndex(iP, iQ),
                                  collision.calculate(f_Previous.haloComputedData(),
                                                      iP-uiL::celerity()[iQ], iQ));
            });

          storeLocalFields(iP, iteration);
        });

      auto t2 = std::chrono::high_resolution_clock::now();

      dtCommunication = (t1 - t0);
      dtComputation = (t2 - t1);
      dtTotal = (t2 - t0);
    }

    using Algorithm<T, AlgorithmType::Generic, architecture>::getCommunicationTime;
    using Algorithm<T, AlgorithmType::Generic, architecture>::getComputationTime;
    using Algorithm<T, AlgorithmType::Generic, architecture>::getTotalTime;
  };


  typedef Algorithm<dataT, algorithmT, architecture> Algorithm_;

}




#endif // ALGORITHM_H
