#ifndef COMPUTE_SERIAL_H
#define COMPUTE_SERIAL_H

#include <array>
#include <vector>
#include <memory>
#include <chrono>
#include <omp.h>

#include "init.h"
#include "structure.h"
#include "commons.h"
#include "lattice.h"
#include "solver.h"
#include "calculate.h"
#include "force.h"
#include "boundary.h"
#include "output.h"
#include "communication.h"

namespace lbm {

  template<class T, LatticeType L>
    void push_fusedCollideAndStream(Lattice<T, L>& l_previous, Lattice<T, L>& l_next,
                                    std::shared_ptr<Solver<T, L>> solver,
                                    std::shared_ptr<Forcing<T, L>> forcing,
                                    Forces<T, L>& forces,
                                    LocalField<T, L>& field,
                                    const int startX, const int iteration) {

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(int iZ = hZ<T, L>(); iZ < hZ<T, L>()+lZ_l<T, L>(); ++iZ) {
      for(int iY = hY<T, L>(); iY < hY<T, L>()+lY_l<T, L>(); ++iY) {
#pragma omp simd
        for(int iX = hX<T, L>(); iX < hX<T, L>()+lX_l<T, L>(); ++iX) {
          BOOST_LOG_TRIVIAL(debug) << " - (x, y, z) = "
                                   << "(" << iX
                                   << ", " << iY
                                   << ", " << iZ << ")";
          int idx_lattice = idxL<T, L>(iX, iY, iZ);
          int idx_field = idx_inF<T, L>(iX, iY, iZ);

          T previousDensity;
          MathVector<T, dimD<T, L>()> previousVelocity;
          calculateMoments<T, L>(l_previous.f_distribution.data(), idx_lattice,
                                 previousDensity, previousVelocity);


          BOOST_LOG_TRIVIAL(debug) << " - Computing force.";
          forcing->force = forces.force(iX-hX<T, L>()+startX, iY-hY<T, L>(), iZ-hZ<T, L>());
          //MathVector<T, dimD<T, L>()> nextForce = forcing->force;

          T previousVelocity2 = previousVelocity.norm2();

          MathVector<T, dimD<T, L>()> previousEqVelocity = previousVelocity
            + forcing->getEqVelocityForcing(previousDensity);

          T previousEqVelocity2 = previousEqVelocity.norm2();

          MathVector<T, dimQ<T, L>()> fNeq;
          MathVector<T, dimQ<T, L>()> fForced;

          UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {
              fForced[iQ] = l_previous.f_distribution[idxPop<T, L>(idx_lattice, iQ)]
                + forcing->getCollisionForcing(iQ, previousDensity,
                                               previousVelocity,
                                               previousVelocity2);

              fNeq[iQ] = l_previous.f_distribution[idxPop<T, L>(idx_lattice, iQ)]
                - computeEquilibrium<T, L>(iQ, previousDensity,
                                           previousEqVelocity,
                                           previousEqVelocity2);
            });

          T nextAlpha = solver->computeAlpha(fForced, fNeq, 2.0);

          BOOST_LOG_TRIVIAL(debug) << " - Colliding and Streaming.";
          UnrolledFor<0, dimQ<T, L>()>::Do([&] (int iQ) {

              int idx_lattice_next = idxL<T, L>(iX + celerity<T, L>(d::X, iQ),
                                                iY + celerity<T, L>(d::Y, iQ),
                                                iZ + celerity<T, L>(d::Z, iQ));
              l_next.f_distribution[idxPop<T, L>(idx_lattice_next, iQ)] = fForced[iQ]
                - nextAlpha*beta*fNeq[iQ];
            });

          //if(iteration%writeStep == 0) {
          BOOST_LOG_TRIVIAL(debug) << " - (x, y, z) = "
                                   << "(" << iX
                                   << ", " << iY
                                   << ", " << iZ << ")";

          BOOST_LOG_TRIVIAL(debug) << " - Storing bulk density, velocities, forces, and alpha.";
          field.nextAlpha[idx_field] = nextAlpha;
          field.previousDensity[idx_field] = previousDensity;
          previousVelocity += forcing->getHydroVelocityForcing(previousDensity);
          field.previousVelocity[idx_field] = previousVelocity;
          // field.nextForce[idx_field] = nextForce;
          //}
        }
      }
    }
  }

  template<class T, LatticeType L>
    void storeNextField(const vector<T, CACHE_LINE>& f_next,
                        Field<T, L>& field, const int startX,
                        std::shared_ptr<Forcing<T, L>>& forcing,
                        Forces<T, L>& forces,
                        const int iteration) {

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(int iZ = hZ<T, L>(); iZ < hZ<T, L>()+lZ_l<T, L>(); ++iZ) {
      for(int iY = hY<T, L>(); iY < hY<T, L>()+lY_l<T, L>(); ++iY) {
#pragma omp simd
        for(int iX = hX<T, L>(); iX < hX<T, L>()+lX_l<T, L>(); ++iX) {

          int idx_lattice = idxL<T, L>(iX, iY, iZ);
          int idx_field = idx_inF<T, L>(iX, iY, iZ);

          T nextDensity;
          MathVector<T, dimD<T, L>()> nextVelocity;
          calculateMoments<T, L>(f_next.data(), idx_lattice,
                                 nextDensity, nextVelocity);

          BOOST_LOG_TRIVIAL(debug) << " - Computing force.";
          forcing->force = forces.force(iX-hX<T, L>()+startX, iY-hY<T, L>(), iZ-hZ<T, L>());

          field.nextDensity[idx_field] = nextDensity;
          nextVelocity += forcing->getHydroVelocityForcing(nextDensity);
          field.nextVelocity[idx_field] = nextVelocity;
        }
      }
    }
  }

  template<class T, LatticeType L>
    void compute(Init<T, L>& init) {

    int startX = 0;
    init.globalField.nextDensity = init.localField.nextDensity;
    init.globalField.nextVelocity = init.localField.nextVelocity;
    init.globalField.nextAlpha = init.localField.nextAlpha;
    init.globalField.nextDistribution = init.localField.nextDistribution;
    init.globalField.previousDensity = init.localField.previousDensity;
    init.globalField.previousVelocity = init.localField.previousVelocity;

    std::shared_ptr<Lattice<T, L> > l0 = std::shared_ptr<Lattice<T, L> >(new Lattice<T, L>(init.localField.unpackDistribution()));
    std::shared_ptr<Lattice<T, L> > l1 = std::shared_ptr<Lattice<T, L> >(new Lattice<T, L>(init.localField.unpackDistribution()));

    double local_startMass = init.localField.getMass();
    double global_startMass = local_startMass;

    double totalTime = 0.;
    double compTime = 0.;
    BOOST_LOG_TRIVIAL(info) << "Simulation starts!";

    BOOST_LOG_TRIVIAL(info) << "Iteratilaon: " << 0
                            << " - Writing outputs.";
    if(startIteration == 0) {
      init.outputs.write(init.globalField, startIteration/writeStep);
    }

    for(int iteration = startIteration+1; iteration <= iterationMax; ++iteration) {

      BOOST_LOG_TRIVIAL(debug) << "Iteration: " << iteration
                               << " - Looping over lattice nodes.";

      auto t0 = std::chrono::high_resolution_clock::now();
      std::swap(l0, l1);

      init.localField.previousVelocity = init.globalField.previousVelocity;

      push_fusedCollideAndStream<T, L>(*l0, *l1, init.solver, init.forcing, init.forces,
                                       init.localField, startX, iteration);
      push_handleHalos<T, L>(*l1);

      auto t1 = std::chrono::high_resolution_clock::now();

      BOOST_LOG_TRIVIAL(debug) << "Iteration: " << iteration
                               << " - Applying boundary conditions over the whole lattice.";
      init.boundaryConditions.apply(*l1);

      auto t2 = std::chrono::high_resolution_clock::now();
      auto dTcomp0 = t1 - t0;
      auto dTtot0 = t2 - t0;
      compTime += dTcomp0.count();
      totalTime += dTtot0.count();

      if(iteration%writeStep == 0) {
        storeNextField<T, L>(l1->f_distribution, init.localField, startX,
                             init.forcing, init.forces, iteration);

        BOOST_LOG_TRIVIAL(info) << "Iteration: " << iteration
                                << " - Writing outputs.";
        init.localField.nextDensity = init.globalField.nextDensity;
        init.localField.nextVelocity = init.globalField.nextVelocity;
        init.localField.nextAlpha = init.globalField.nextAlpha;
        init.localField.previousDensity = init.globalField.previousDensity;
        init.localField.previousVelocity = init.globalField.previousVelocity;

        if(iteration%backupStep == 0) {
          BOOST_LOG_TRIVIAL(info) << "Iteration: " << iteration
                                  << " - Writing backup.";
          init.localField.packDistribution(*l1);
          init.localField.nextDistribution = init.globalField.nextDistribution;
        }

        init.outputs.write(init.globalField, iteration/writeStep);
      }

    }

    double local_endMass = init.localField.getMass();
    double global_endMass = local_endMass;

    const double commTime = totalTime - compTime;

    const long nCells = size_g;
    const double mlups = (1e-6*nCells) / (totalTime / iterationMax);

    const double relative_MassDiff = fabs(global_endMass-global_startMass)/global_startMass;

    std::cout << "Prefix     : " << prefix << std::endl;
    std::cout << "NPROCS     : " << NPROCS << std::endl;
    std::cout << "NTHREADS   : " << NTHREADS << std::endl;
    std::cout << "Total time : " << totalTime << " s" << std::endl;
    std::cout << "Comp time  : " << compTime << " s" << std::endl;
    std::cout << "Comm time  : " << commTime << " s" << std::endl;
    std::cout << "MLUPS      : " << mlups << std::endl;
    std::cout << "Mass diff  : " << relative_MassDiff << std::endl;

  }
}

#endif // COMPUTE_SERIAL_H
