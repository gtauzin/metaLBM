#include "compute.h"

#include "input.h"
#include "commons.h"
#include "helpers.h"
#include "init.h"
#include "communication.h"
#include "lattice.h"
#include "boundary.h"
#include "force.h"
#include "output.h"

extern "C" {
#include "ex1.h"
#include "kernels.h"
}

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <array>
#include <memory>
#include <math.h>
#include <chrono>
#include <algorithm>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

namespace lbm {

  void push_fusedCollideAndStream(Lattice& l_previous, Lattice& l_next,
                                  std::shared_ptr<Solver> solver,
                                  std::shared_ptr<Forcing>& forcing, Forces& forces,
                                  LocalField& field, const int startX, const int iteration) {
#pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(int iX = haloX; iX < haloX+lengthX_l; ++iX) {
      for(int iY = haloY; iY < haloY+lengthY_l; ++iY) {
        BOOST_LOG_TRIVIAL(debug) << " - (x, y) = (" << iX << ", " << iY << ")";
        int idx_lattice = idxL(iX, iY);
        int idx_field = idx_inF(iX, iY);

        double previousDensity, previousVelocityX, previousVelocityY;
        calculateMoments(l_previous.f_distribution, idx_lattice,
                         previousDensity, previousVelocityX, previousVelocityY);

        BOOST_LOG_TRIVIAL(debug) << " - Computing force.";
        forcing->forceX = forces.forceX(iX-haloX+startX, iY-haloY, iteration);
        forcing->forceY = forces.forceY(iX-haloX+startX, iY-haloY, iteration);

        previousVelocityX += forcing->getVelocityXForcing(previousDensity);
        previousVelocityY += forcing->getVelocityYForcing(previousDensity);
        double previousVelocity2 = previousVelocityX * previousVelocityX
                                 + previousVelocityY * previousVelocityY;

        std::array<double, 9> fNeq;
        std::array<double, 9> f;
        for(int i = 0; i < 9; ++i) {
          f[i] = l_previous.f_distribution[idxPop(idx_lattice, i)]
            + forcing->getCollisionForcing(i, previousDensity,
                                           previousVelocityX, previousVelocityY, previousVelocity2);
          fNeq[i] = f[i] - computeEquilibrium(i, previousDensity,
                                              previousVelocityX, previousVelocityY, previousVelocity2);
        }

        double nextAlpha = solver->computeAlpha(f, fNeq, 2.0);

        BOOST_LOG_TRIVIAL(debug) << " - Colliding and Streaming.";
        for(int i = 0; i < 9; ++i) {
          int idx_lattice_next = idxL(iX + celerityX[i], iY + celerityY[i]);
          l_next.f_distribution[idxPop(idx_lattice_next, i)] = f[i] - nextAlpha*beta*fNeq[i];
        }

        if(iteration%writeStep == 0) {
          BOOST_LOG_TRIVIAL(debug) << " - (x, y) = (" << iX << ", " << iY << ")"
                                   << " - Storing bulk density, velocity and alpha.";
          field.nextAlpha[idx_field] = nextAlpha;
          field.previousDensity[idx_field] = previousDensity;
          field.previousVelocityX[idx_field] = previousVelocityX;
          field.previousVelocityY[idx_field] = previousVelocityY;

        }
      }
    }

  }

  void storeNextField(const std::vector<double>& f_next, Field& field, const int startX,
                      std::shared_ptr<Forcing>& forcing, Forces& forces, const int iteration) {
    for(int iX = haloX; iX < haloX+lengthX_l; ++iX) {
      for(int iY = haloY; iY < haloY+lengthY_l; ++iY) {

        int idx_lattice = idxL(iX, iY);
        int idx_field = idx_inF(iX, iY);

        double nextDensity, nextVelocityX, nextVelocityY;
        calculateMoments(f_next, idx_lattice,
                         nextDensity, nextVelocityX, nextVelocityY);

        BOOST_LOG_TRIVIAL(debug) << " - Computing force.";
        forcing->forceX = forces.forceX(iX-haloX+startX, iY-haloY, iteration+1);
        forcing->forceY = forces.forceY(iX-haloX+startX, iY-haloY, iteration+1);

        field.nextDensity[idx_field] = nextDensity;
        field.nextVelocityX[idx_field] = nextVelocityX
                                       + forcing->getVelocityXForcing(nextDensity);
        field.nextVelocityY[idx_field] = nextVelocityY
                                       + forcing->getVelocityYForcing(nextDensity);
      }
    }
  }

  void compute(Init& init, const int mpi_rank) {
    int mpi_left  = ( mpi_rank + NPROCS - 1 ) % NPROCS;
    int mpi_right = ( mpi_rank          + 1 ) % NPROCS;

    std::array<double, size_buf> sendLeft_buffer;
    std::array<double, size_buf> sendRight_buffer;
    std::array<double, size_buf> receiveLeft_buffer;
    std::array<double, size_buf> receiveRight_buffer;
    MPI_Status status;

    MPI_Scatter(init.globalField.nextDensity.data(), size_l, MPI_DOUBLE,
                init.localField.nextDensity.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.nextVelocityX.data(), size_l, MPI_DOUBLE,
                init.localField.nextVelocityX.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.nextVelocityY.data(), size_l, MPI_DOUBLE,
                init.localField.nextVelocityY.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.nextAlpha.data(), size_l, MPI_DOUBLE,
                init.localField.nextAlpha.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.nextDistribution.data(), 9*size_l, MPI_DOUBLE,
                init.localField.nextDistribution.data(), 9*size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.previousDensity.data(), size_l, MPI_DOUBLE,
                init.localField.previousDensity.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.previousVelocityX.data(), size_l, MPI_DOUBLE,
                init.localField.previousVelocityX.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.previousVelocityY.data(), size_l, MPI_DOUBLE,
                init.localField.previousVelocityY.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    std::shared_ptr<Lattice> l0 = std::shared_ptr<Lattice>(new Lattice(init.localField.unpackDistribution()));
    std::shared_ptr<Lattice> l1 = std::shared_ptr<Lattice>(new Lattice(init.localField.unpackDistribution()));

    const int startX = lengthX_l*mpi_rank;

    double local_startMass = init.localField.getMass();
    MPI_Barrier(MPI_COMM_WORLD);
    double global_startMass;
    MPI_Reduce(&local_startMass, &global_startMass, 1, MPI_DOUBLE,
               MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    double totalTime = 0.;
    double compTime = 0.;
    BOOST_LOG_TRIVIAL(info) << "Simulation starts!";

    BOOST_LOG_TRIVIAL(info) << "Iteratilaon: " << 0
                            << " - Writing outputs.";
    if(mpi_rank == 0) init.outputs.write(init.globalField, 0);

    for(int iteration = startIteration; iteration <= iterationMax; ++iteration) {

      BOOST_LOG_TRIVIAL(debug) << "Iteration: " << iteration
                               << " - Looping over lattice nodes.";

      auto t0 = std::chrono::high_resolution_clock::now();
      std::swap(l0, l1);

      push_fusedCollideAndStream(*l0, *l1, init.solver, init.forcing, init.forces,
                                 init.localField, startX, iteration);

      dummyPropagateDevice((*l1).f_distribution.data(), (*l0).f_distribution.data(),
                           haloX, haloX+lengthX_g, haloY, haloY+lengthY_g);
      sayHi(mpi_rank, 2);

      push_handleHalos(*l1);

      auto t1 = std::chrono::high_resolution_clock::now();

      push_packToSend_mpi(*l1, sendLeft_buffer, sendRight_buffer);
      MPI_Sendrecv(sendRight_buffer.data(), size_buf, MPI_DOUBLE, mpi_right, 17,
                   receiveLeft_buffer.data(), size_buf, MPI_DOUBLE, mpi_left, 17,
                   MPI_COMM_WORLD, &status);
      MPI_Sendrecv(sendLeft_buffer.data(), size_buf, MPI_DOUBLE, mpi_left, 23,
                   receiveRight_buffer.data(), size_buf, MPI_DOUBLE, mpi_right, 23,
                   MPI_COMM_WORLD, &status);
      push_unpackReceived_mpi(*l1, receiveLeft_buffer, receiveRight_buffer);

      BOOST_LOG_TRIVIAL(debug) << "Iteration: " << iteration
                               << " - Applying boundary conditions over the whole lattice.";
      init.boundaryConditions.apply(*l1);

      auto t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> dTcomp0 = t1 - t0;
      std::chrono::duration<double> dTtot0 = t2 - t0;
      compTime += dTcomp0.count();
      totalTime += dTtot0.count();

      if(iteration%writeStep == 0) {
        storeNextField(l1->f_distribution, init.localField, startX,
                       init.forcing, init.forces, iteration);

        BOOST_LOG_TRIVIAL(info) << "Iteration: " << iteration
                                << " - Writing outputs.";
        MPI_Gather(init.localField.nextDensity.data(), size_l, MPI_DOUBLE,
                   init.globalField.nextDensity.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        MPI_Gather(init.localField.nextVelocityX.data(), size_l, MPI_DOUBLE,
                   init.globalField.nextVelocityX.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        MPI_Gather(init.localField.nextVelocityY.data(), size_l, MPI_DOUBLE,
                   init.globalField.nextVelocityY.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        MPI_Gather(init.localField.nextAlpha.data(), size_l, MPI_DOUBLE,
                   init.globalField.nextAlpha.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        MPI_Gather(init.localField.previousDensity.data(), size_l, MPI_DOUBLE,
                   init.globalField.previousDensity.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        MPI_Gather(init.localField.previousVelocityX.data(), size_l, MPI_DOUBLE,
                   init.globalField.previousVelocityX.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        MPI_Gather(init.localField.previousVelocityY.data(), size_l, MPI_DOUBLE,
                   init.globalField.previousVelocityY.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);

        if(iteration%backupStep == 0) {
          BOOST_LOG_TRIVIAL(info) << "Iteration: " << iteration
                                  << " - Writing backup.";
          init.localField.packDistribution(*l1);
          MPI_Gather(init.localField.nextDistribution.data(), 9*size_l, MPI_DOUBLE,
                     init.globalField.nextDistribution.data(), 9*size_l, MPI_DOUBLE,
                     0, MPI_COMM_WORLD);
        }

        if(mpi_rank == 0) init.outputs.write(init.globalField, iteration/writeStep);
      }

    }

    double local_endMass = init.localField.getMass();
    MPI_Barrier(MPI_COMM_WORLD);
    double global_endMass;
    MPI_Reduce(&local_endMass, &global_endMass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if(mpi_rank == 0) {
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


}
