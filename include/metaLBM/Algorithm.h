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

  template<class T, AlgorithmType>
  class Algorithm {};

  template<class T>
  class Algorithm<T, AlgorithmType::Generic> {
  protected:
    Communication_& communication;

    Field<T, 1, writeDensity>& densityField;
    Field<T, L::dimD, writeVelocity>& velocityField;
    Field<T, L::dimD, writeDensity>& forceField;
    Field<T, 1, writeDensity>& alphaField;
    Distribution<T>& f_Previous;
    Distribution<T>& f_Next;

    Collision_ collision;
    Moment<T> moment;
    Boundary_ boundary;

    std::chrono::duration<double> dtComputation;
    std::chrono::duration<double> dtCommunication;
    std::chrono::duration<double> dtTotal;

    Algorithm(Communication_& communication_in,
              Field<T, 1, writeDensity>& densityField_in,
              Field<T, L::dimD, writeVelocity>& velocityField_in,
              Field<T, L::dimD, writeDensity>& forceField_in,
              Field<T, 1, writeDensity>& alphaField_in,
              Distribution<T>& f_Previous_in,
              Distribution<T>& f_Next_in)
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

  template<class T>
  class Algorithm<T, AlgorithmType::Pull>
    : public Algorithm<T, AlgorithmType::Generic> {
  private:
    using Algorithm<T, AlgorithmType::Generic>::communication;

    using Algorithm<T, AlgorithmType::Generic>::densityField;
    using Algorithm<T, AlgorithmType::Generic>::velocityField;
    using Algorithm<T, AlgorithmType::Generic>::forceField;
    using Algorithm<T, AlgorithmType::Generic>::alphaField;
    using Algorithm<T, AlgorithmType::Generic>::f_Previous;
    using Algorithm<T, AlgorithmType::Generic>::f_Next;

    using Algorithm<T, AlgorithmType::Generic>::collision;
    using Algorithm<T, AlgorithmType::Generic>::moment;
    using Algorithm<T, AlgorithmType::Generic>::boundary;

    using Algorithm<T, AlgorithmType::Generic>::dtComputation;
    using Algorithm<T, AlgorithmType::Generic>::dtCommunication;
    using Algorithm<T, AlgorithmType::Generic>::dtTotal;

    using Algorithm<T, AlgorithmType::Generic>::storeLocalFields;

  public:
    Algorithm(Communication_& communication_in,
              Field<T, 1, writeDensity>& densityField_in,
              Field<T, L::dimD, writeVelocity>& velocityField_in,
              Field<T, L::dimD, writeDensity>& forceField_in,
              Field<T, 1, writeDensity>& alphaField_in,
              Distribution<T>& f_Previous_in,
              Distribution<T>& f_Next_in)
      : Algorithm<T, AlgorithmType::Generic>(communication_in,
                                             densityField_in, velocityField_in,
                                             forceField_in, alphaField_in,
                                             f_Previous_in, f_Next_in)
    {}

    void iterate(const unsigned int iteration) {
      SCOREP_INSTRUMENT_ON("Algorithm<T, AlgorithmType::Pull>::iterate")

      f_Previous.swapHalo(f_Next);

      //force.update(iteration);

      auto t0 = std::chrono::high_resolution_clock::now();
      communication.periodic(f_Previous.haloData());

      //boundary.apply(f_Previous.haloData());

      auto t1 = std::chrono::high_resolution_clock::now();

      Computation::Do(lD::start()+L::halo(), lD::end()+L::halo(),
                      [&] (MathVector<unsigned int, 3>& iP) {
          SCOREP_INSTRUMENT_ON("Algorithn<T, AlgorithmType::Pull>::lambda[fused_collide_and_push]")

          moment.calculateMoments(f_Previous.haloData(), iP);

          collision.setForce(iP+gD::offset(communication.getRankMPI()));
          collision.setVariables(f_Previous.haloData(), iP,
                                 moment.getDensity(), moment.getVelocity());

          UnrolledFor<0, L::dimQ>::Do([&] (unsigned int iQ) {
              f_Next.setHaloField(hD::getIndex(iP, iQ),
                                  collision.calculate(f_Previous.haloData(),
                                                      iP-uiL::celerity()[iQ], iQ));
            });

          storeLocalFields(iP, iteration);
        });

      auto t2 = std::chrono::high_resolution_clock::now();

      dtCommunication = (t1 - t0);
      dtComputation = (t2 - t1);
      dtTotal = (t2 - t0);
    }

    using Algorithm<T, AlgorithmType::Generic>::getCommunicationTime;
    using Algorithm<T, AlgorithmType::Generic>::getComputationTime;
    using Algorithm<T, AlgorithmType::Generic>::getTotalTime;
  };


  typedef Algorithm<dataT, algorithmT> Algorithm_;

}




#endif // ALGORITHM_H
