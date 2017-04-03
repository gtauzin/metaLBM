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
#include "field.h"
#include "distribution.h"
#include "solver.h"
#include "calculate.h"
#include "force.h"
#include "boundary.h"
#include "output.h"
#include "communication.h"

namespace lbm {

  template<class T>
    void pull_fusedCollideAndStream(Distribution<T>& f_previous, Distribution<T>& f_next,
                                    std::shared_ptr<Solver<T>> solver,
                                    std::shared_ptr<ForcingScheme<T>> forcingScheme,
                                    Forces<T>& forces,
                                    LocalField<T>& field,
                                    const int startX, const int iteration) {

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(int iZ = L::hZ; iZ < L::hZ + L::lZ_l; ++iZ) {
      for(int iY = L::hY; iY < L::hY + L::lY_l; ++iY) {
#pragma omp simd
        for(int iX = L::hX; iX < L::hX + L::lX_l; ++iX) {
          BOOST_LOG_TRIVIAL(debug) << " - (x, y, z) = "
                                   << "(" << iX
                                   << ", " << iY
                                   << ", " << iZ << ")";
            MathVector<int, 3> iP{{iX, iY, iZ}};
            int idx_lattice = idxL(iX, iY, iZ);
            int idx_field = idx_inF(iX, iY, iZ);

          T previousDensity;
          MathVector<T, L::dimD> previousVelocity;
          calculateMoments<T>(f_previous.fPop.data(), idx_lattice,
                                 previousDensity, previousVelocity);

          MathVector<int, 3> iP_lF{iX-L::hX+startX, iY-L::hY, iZ-L::hZ};
          BOOST_LOG_TRIVIAL(debug) << " - Computing force.";
          forcingScheme->force = forces.force(iP_lF);
          //MathVector<T, L::dimD> nextForce = forcing->force;

          T previousVelocity2 = previousVelocity.norm2();

          MathVector<T, L::dimD> previousEqVelocity = previousVelocity
            + forcingScheme->getEqVelocityForcing(previousDensity);

          T previousEqVelocity2 = previousEqVelocity.norm2();

          MathVector<T, L::dimQ> fNeq;
          MathVector<T, L::dimQ> fForced;

          UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {
              int idx_lattice_previous = idxL(iP - L::celerity()[iQ]);
              fForced[iQ] = f_previous[idxPop(idx_lattice_previous, iQ)]
                + forcingScheme->getCollisionForcing(iQ, previousDensity,
                                               previousVelocity,
                                               previousVelocity2);

              fNeq[iQ] = f_previous[idxPop(idx_lattice_previous, iQ)]
                - computeEquilibrium<T>(iQ, previousDensity,
                                           previousEqVelocity,
                                           previousEqVelocity2);
            });

          T nextAlpha = solver->computeAlpha(fForced, fNeq, 2.0);

          BOOST_LOG_TRIVIAL(debug) << " - Colliding and Streaming.";
          UnrolledFor<0, L::dimQ>::Do([&] (int iQ) {

              f_next[idxPop(idx_lattice, iQ)] = fForced[iQ]
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
          previousVelocity += forcingScheme->getHydroVelocityForcing(previousDensity);
          field.previousVelocity[idx_field] = previousVelocity;
          // field.nextForce[idx_field] = nextForce;
          //}
        }
      }
    }
  }

  template<class T>
    void storeNextField(const vector<T, CACHE_LINE>& f_next,
                        Field<T>& field, const int startX,
                        std::shared_ptr<ForcingScheme<T>>& forcingScheme,
                        Forces<T>& forces,
                        const int iteration) {

#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(int iZ = L::hZ; iZ < L::hZ + L::lZ_l; ++iZ) {
      for(int iY = L::hY; iY < L::hY + L::lY_l; ++iY) {
#pragma omp simd
        for(int iX = L::hX; iX < L::hX + L::lX_l; ++iX) {

          int idx_lattice = idxL(iX, iY, iZ);
          int idx_field = idx_inF(iX, iY, iZ);

          T nextDensity;
          MathVector<T, L::dimD> nextVelocity;
          calculateMoments<T>(f_next.data(), idx_lattice,
                                 nextDensity, nextVelocity);

          MathVector<int, 3> iP_lF{iX-L::hX+startX, iY-L::hY, iZ-L::hZ};
          BOOST_LOG_TRIVIAL(debug) << " - Computing force.";
          forcingScheme->force = forces.force(iP_lF);

          field.nextDensity[idx_field] = nextDensity;
          nextVelocity += forcingScheme->getHydroVelocityForcing(nextDensity);
          field.nextVelocity[idx_field] = nextVelocity;
        }
      }
    }
  }

  template<class T>
    void compute(Init<T>& init) {

    int startX = 0;
    init.globalField.nextDensity = init.localField.nextDensity;
    init.globalField.nextVelocity = init.localField.nextVelocity;
    init.globalField.nextAlpha = init.localField.nextAlpha;
    init.globalField.nextDistribution = init.localField.nextDistribution;
    init.globalField.previousDensity = init.localField.previousDensity;
    init.globalField.previousVelocity = init.localField.previousVelocity;

    std::shared_ptr<Distribution<T> > f0 = std::shared_ptr<Distribution<T> >(new Distribution<T>(init.localField.unpackDistribution()));
    std::shared_ptr<Distribution<T> > f1 = std::shared_ptr<Distribution<T> >(new Distribution<T>(init.localField.unpackDistribution()));

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
      std::swap(f0, f1);

      pull_handleHalos<T>(*f0);

      pull_fusedCollideAndStream<T>(*f0, *f1, init.solver, init.forcingScheme, init.forces,
                                       init.localField, startX, iteration);
      auto t1 = std::chrono::high_resolution_clock::now();

      BOOST_LOG_TRIVIAL(debug) << "Iteration: " << iteration
                               << " - Applying boundary conditions over the whole lattice.";
      init.boundaryConditions.apply(*f1);

      auto t2 = std::chrono::high_resolution_clock::now();
      auto dTcomp0 = t1 - t0;
      auto dTtot0 = t2 - t0;
      compTime += dTcomp0.count();
      totalTime += dTtot0.count();

      if(iteration%writeStep == 0) {
        storeNextField<T>(f1->fPop, init.localField, startX,
                             init.forcingScheme, init.forces, iteration);

        BOOST_LOG_TRIVIAL(info) << "Iteration: " << iteration
                                << " - Writing outputs.";
        init.globalField.nextDensity = init.localField.nextDensity;
        init.globalField.nextVelocity = init.localField.nextVelocity;
        init.globalField.nextAlpha = init.localField.nextAlpha;
        init.globalField.nextDistribution = init.localField.nextDistribution;
        init.globalField.previousDensity = init.localField.previousDensity;
        init.globalField.previousVelocity = init.localField.previousVelocity;

        if(iteration%backupStep == 0) {
          BOOST_LOG_TRIVIAL(info) << "Iteration: " << iteration
                                  << " - Writing backup.";
          init.localField.packDistribution((*f1).fPop);
          init.localField.nextDistribution = init.globalField.nextDistribution;
        }

        init.outputs.write(init.globalField, iteration/writeStep);
      }

    }

    double local_endMass = init.localField.getMass();
    double global_endMass = local_endMass;

    const double commTime = totalTime - compTime;

    const long nCells = s_g();
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
