#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <chrono>
#include <omp.h>

#include "Init.h"
#include "Options.h"
#include "Parameters.h"
#include "Field.h"
#include "Distribution.h"
#include "Collider.h"
#include "Moments.h"
#include "ForcindScheme.h"
#include "Boundary.h"
#include "Write.h"
#include "Communication.h"
#include "Computation.h"

namespace lbm {

  template<class T, AlgorithmType algorithmType,
           PartitionningType partitionningType, MemoryLayout memoryLayout>
  class Algorithm {};

  template<>
  class Algorithm<class T, AlgorithmType::Generic> {
  protected:
    // L::latticeType ou P::latticeType?????
    Domain<L::latticeType, DomainType::Local, partionningType, memoryLayout> D;
    double mass;

    std::chrono::duration<double> dTcomputation;
    std::chrono::duration<double> dTcommunication;
    std::chrono::duration<double> dTtotal;

    double computeTime;
    double communicationTime;
    double writeTime;
    double totalTime;

  Algorithm(const Parameters& parameters_in,
            const unsigned int rankMPI_in)
    : rankMPI(rankMPI_in);
    , mass(communicateMass())
    , computeTime(0.0)
    , communicationTime(0.0)
    , writeTime(0.0)
    , totalTime(0.0)
    {
      // Probably when field is initialized...
      communicateFromGlobalToLocal();

      if(rankMPI == 0 && startIteration == 0) {
        init.outputs.write(init.globalField, startIteration/writeStep);
      }
    }


    void computeLBM() {
      printInputs();

      BOOST_LOG_TRIVIAL(info) << "Simulation starts!\n\n";
      for(int iteration = startIteration+1; iteration <= endIteration; ++iteration) {

        iterateLBM();

        writeFields(iteration);

        computeTime += dTcomp0.count();
        totalTime += dTtot0.count();

        communicationTime = totalTime - computeTime;

      }

      mass = fabs(mass-communicateMass())/mass;
      printOutputs();
    }

  private:
    void writeFields(const int iteration) {
      if(iteration%writeStep == writeStep-1) {
        communicateLocalPreviousFieldsToGlobalPrevious2Fields();
      }

      if(iteration%writeStep == 0) {
        storeNextField(l1->f_distribution, init.localField, startX,
                       init.forcing, init.forces, iteration);

        communicateLocalFieldsToGlobalFieldsExceptPrevious2();

        if(iteration%backupStep == 0) {
          communicateDistribution();
        }

          if(rankMPI == 0) {
            init.outputs.write(init.globalField, iteration/writeStep);
          }
      }
    }



    void printInputs(const unsigned int sizeMPI,
                         std::string processorName) {

      std::cout << "MPI #" << rankMPI << " of " << sizeMPI
                << " processes running on host " << processorName << "\n";

      if (rankMPI == 0) {
        std::cout.precision(17);
        std::cout << "-------------------OPTIONS-------------------" << std::endl;
        std::cout << " Lattice: D" << L::dimD std::cout << "Q" << L::dimQ << std::endl;
        std::cout << " Global length Y: " << L::lY_g << std::endl;
        std::cout << " Global length X: " << L::lX_g << std::endl;
        std::cout << " Global length Y: " << L::lY_g << std::endl;
        std::cout << " Global length Z: " << L::lZ_g << std::endl;
        std::cout << " Global memory: " << (int)(L::lX_g*L::lY_g*L::lZ_g
                                                 *sizeof(dataType))/(1<<30) << std::endl;
        std::cout << "MPI #" << mpi_rank << std::endl;
        std::cout << " Local length X : " << L::lX_l << std::endl;
        std::cout << " Local length Y : " << L::lY_l << std::endl;
        std::cout << " Local length Z : " << L::lZ_l << std::endl;
        std::cout << " Local memory : " << (int)(L::lZ_l*L::lY_l*L::lZ_l
                                                 *sizeof(dataType))/(1<<30) << std::endl;
        std::cout << "---------------------------------------------\n" << std::endl;
        std::cout << "NPROCS: " << NPROCS << "" << std::endl;
        std::cout << "NTHREADS: " << NTHREADS << "" << std::endl;
        std::cout << "STRUCTURE:   SOA " << std::endl;
        std::cout << "MAX ITERATION: " << iterationMax << "\n" << std::endl;
        std::cout << "-------------------PARAMETERS-----------------" << std::endl;
        std::cout << "tau: " << tau << std::endl;
        std::cout << ", correponding viscosity: " << L::cs2 * (tau - 0.5) << std::endl;
        std::cout << "beta:" << beta << std::endl;
        std::cout << ", correponding viscosity: " << L::cs2 * (1/(2 * beta) - 0.5)
                  << std::endl;
      }
    }

    void printOutputs() {
      if (rankMPI == 0) {
        std::cout << "--------------------OUTPUTS------------------" << std::endl;
        std::cout << "Total time : " << totalTime << " s" << std::endl;
        std::cout << "Comp time  : " << compTime << " s" << std::endl;
        std::cout << "Comm time  : " << commTime << " s" << std::endl;

        const double mlups = (L::lX_g*L::lY_g*L::lZ_g*1e-6*)/(totalTime/iterationMax);

        std::cout << "MLUPS      : " << mlups << std::endl;
        std::cout << "Mass diff  : " << mass << std::endl;
        std::cout << "---------------------------------------------\n" << std::endl;
      }
    }

  };

  template<>
  class Algorithm<class T, AlgorithmType::FusedPull>:
    public Algorithm<T, AlgorithmType::Generic> {
  public:
    using Algorithm<T, AlgorithmType::Generic>::Algorithm;

  private:
    using Algorithm<T, AlgorithmType::Generic>::Domain;

    Distribution f_Previous;
    Distribution f_Next;

    void iterateLBM() {
      f_Previous.swap(f_Next);

      force.updateTime(iteration);

      BOOST_LOG_TRIVIAL(debug) << "Iteration: " << iteration
                               << " - Communicating halos to neighboring processes.";
      auto t1 = std::chrono::high_resolution_clock::now();
      communicateHalos(*l1);

      BOOST_LOG_TRIVIAL(debug) << "Iteration: " << iteration
                               << " - Computing post-collision distribution and streaming.";
      auto t0 = std::chrono::high_resolution_clock::now();
      collideAndStream(f_Previous, F_B, init.solver, init.forcing, init.forces,
                       init.localField, iteration);

      BOOST_LOG_TRIVIAL(debug) << "Iteration: " << iteration
                               << " - Applying boundary conditions over the whole lattice.";
      init.boundaryConditions.apply(*f_Next);

      auto t2 = std::chrono::high_resolution_clock::now();
    }

    // Wrong! getIndex for Local Domain!!!!!!
    // Halo Domain and Local domain together... How?

    void stream(const MathVector<int, 3>& iP) {
      D::updatePullIndex(iP);

      BOOST_LOG_TRIVIAL(debug) << " - Streaming.";
      UnrolledFor<0, L::dimQ>::Do([&] (unsigned int iQ) {
          f_Next[D::getPullIndex(iQ)] = f_Previous[D::getIndex(iQ)];
      });
    }

    void collideAndStream(const MathVector<int, 3>& iP,
                          const unsigned int iteration) {

      D::updateIndex(iP);

      moments.computeMoments()
      T density = computeDensity(f_Previous.data(), D::index);
      MathVector<T, L::dimD> velocity = computeVelocity<T>(f_Previous.data(),
                                                           D::index, density);

      // TODO: Make sure that no variables are set if there are no forces or no forcingScheme!
      P::forcingScheme.setVariables(forces.force(iP), density, velocity);
      P::equilibrium.setVariables(density,
                                  P::forcingScheme.getEqVelocityForcing(velocity,
                                                                        density));

      P::collisionOperator.collide(f_Previous, f_Next, D::index);

      stream();

      BOOST_LOG_TRIVIAL(debug) << " - Storing data.";
      density.store(density, D::index, iteration);
      velocity.store(velocity, D::index, iteration);
      alpha.store(collisionOperator.getAlpha(), D::index, iteration);
      force.store(forcingScheme.getForce(), D::index, iteration);
      force.store(moments.calculateEntopy(f_), D::index, iteration); // ?
      //Store everything here including previous and previous2 fields!!!!!
      // Or find a way for it to be done automatically with a TimeField?
    }

  };

}


#endif // ALGORITHM_H
