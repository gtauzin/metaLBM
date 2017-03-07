#include "compute.h"

#include "input.h"
#include "commons.h"
#include "helpers.h"
#include "analysis.h"
#include "init.h"
#include "communication.h"
#include "lattice.h"
#include "boundary.h"
#include "force.h"
#include "output.h"

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <array>
#include <memory>
#include <math.h>
#include <chrono>
#include <algorithm>
#include <complex>

#include <boost/align/aligned_allocator.hpp>
#include <vector>
template<class T, std::size_t Alignment = 1>
#ifdef ALIGN_MEMORY
  using vector = std::vector<T,
  boost::alignment::aligned_allocator<T, Alignment> >;
#else
using vector = std::vector<T>;
#endif

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

namespace lbm {

  void calculateMomentsField(Lattice& l_previous, LocalField& field) {
    #pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(int iX = haloX; iX < haloX+lengthX_l; ++iX) {
      #pragma omp simd
      for(int iY = haloY; iY < haloY+lengthY_l; ++iY) {
        BOOST_LOG_TRIVIAL(debug) << " - (x, y) = (" << iX << ", " << iY << ")";
        int idx_lattice = idxL(iX, iY);
        int idx_field = idx_inF(iX, iY);

        double previousDensity, previousVelocityX, previousVelocityY;
        calculateMoments(l_previous.f_distribution.data(), idx_lattice,
                         previousDensity, previousVelocityX, previousVelocityY);

        field.nextVelocityX[idx_field] = previousVelocityX;
        field.nextVelocityY[idx_field] = previousVelocityY;
      }
    }

  }


  void push_fusedCollideAndStream(Lattice& l_previous, Lattice& l_next,
                                  std::shared_ptr<Solver> solver,
                                  std::shared_ptr<Forcing>& forcing, Forces& forces,
                                  LocalField& field, const int startX, const int iteration) {

    forces.calculateTimeFactor(iteration);

    #pragma omp parallel for schedule(static) num_threads(NTHREADS)
    for(int iX = haloX; iX < haloX+lengthX_l; ++iX) {
      #pragma omp simd
      for(int iY = haloY; iY < haloY+lengthY_l; ++iY) {
        BOOST_LOG_TRIVIAL(debug) << " - (x, y) = (" << iX << ", " << iY << ")";
        int idx_lattice = idxL(iX, iY);
        int idx_field = idx_inF(iX, iY);

        double previousDensity, previousVelocityX, previousVelocityY;
        calculateMoments(l_previous.f_distribution.data(), idx_lattice,
                         previousDensity, previousVelocityX, previousVelocityY);


        BOOST_LOG_TRIVIAL(debug) << " - Computing force.";
        forcing->forceX = forces.forceX(iX-haloX+startX, iY-haloY, iteration);
        forcing->forceY = forces.forceY(iX-haloX+startX, iY-haloY, iteration);
        double nextForceX = forcing->forceX;
        double nextForceY = forcing->forceY;

        double previousVelocity2 = previousVelocityX * previousVelocityX
          + previousVelocityY * previousVelocityY;

        double previousEqVelocityX = previousVelocityX + forcing->getEqVelocityXForcing(previousDensity);
        double previousEqVelocityY = previousVelocityY + forcing->getEqVelocityYForcing(previousDensity);

        double previousEqVelocity2 = previousEqVelocityX * previousEqVelocityX
          + previousEqVelocityY * previousEqVelocityY;

        std::array<double, 9> fNeq;
        std::array<double, 9> fForced;
        for(int i = 0; i < 9; ++i) {
          fForced[i] = l_previous.f_distribution[idxPop(idx_lattice, i)]
            + forcing->getCollisionForcing(i, previousDensity,
                                           previousVelocityX, previousVelocityY,
                                           previousVelocity2);

          fNeq[i] = l_previous.f_distribution[idxPop(idx_lattice, i)]
            - computeEquilibrium(i, previousDensity,
                                 previousEqVelocityX, previousEqVelocityY,
                                 previousEqVelocity2);
        }

        double nextAlpha = solver->computeAlpha(fForced, fNeq, 2.0);

        BOOST_LOG_TRIVIAL(debug) << " - Colliding and Streaming.";
        for(int i = 0; i < 9; ++i) {
          int idx_lattice_next = idxL(iX + celerityX[i], iY + celerityY[i]);
          l_next.f_distribution[idxPop(idx_lattice_next, i)] = fForced[i] - nextAlpha*beta*fNeq[i];
        }

        //if(iteration%writeStep == 0) {
          BOOST_LOG_TRIVIAL(debug) << " - (x, y) = (" << iX << ", " << iY << ")"
                                   << " - Storing bulk density, velocities, forces, and alpha.";
          field.nextAlpha[idx_field] = nextAlpha;
          field.nextForceX[idx_field] = nextForceX;
          field.nextForceY[idx_field] = nextForceY;
          field.nextEntropy[idx_field] = calculateEntropy(l_next, iX, iY);
          field.previousDensity[idx_field] = previousDensity;
          field.previousVelocityX[idx_field] = previousVelocityX
             + forcing->getHydroVelocityXForcing(previousDensity);
          field.previousVelocityY[idx_field] = previousVelocityY
             + forcing->getHydroVelocityYForcing(previousDensity);
          field.previousEntropy[idx_field] = calculateEntropy(l_previous, iX, iY);
          //}
      }
    }

  }

  void storeNextField(const vector<double, CACHE_LINE>& f_next,
                      Field& field, const int startX,
                      std::shared_ptr<Forcing>& forcing, Forces& forces,
                      const int iteration)
  {

    forces.calculateTimeFactor(iteration+1);

    for(int iX = haloX; iX < haloX+lengthX_l; ++iX) {
      for(int iY = haloY; iY < haloY+lengthY_l; ++iY) {

        int idx_lattice = idxL(iX, iY);
        int idx_field = idx_inF(iX, iY);

        double nextDensity, nextVelocityX, nextVelocityY;
        calculateMoments(f_next.data(), idx_lattice,
                         nextDensity, nextVelocityX, nextVelocityY);

        BOOST_LOG_TRIVIAL(debug) << " - Computing force.";
        forcing->forceX = forces.forceX(iX-haloX+startX, iY-haloY, iteration+1);
        forcing->forceY = forces.forceY(iX-haloX+startX, iY-haloY, iteration+1);

        field.nextDensity[idx_field] = nextDensity;
        field.nextVelocityX[idx_field] = nextVelocityX
                                       + forcing->getHydroVelocityXForcing(nextDensity);
        field.nextVelocityY[idx_field] = nextVelocityY
                                       + forcing->getHydroVelocityYForcing(nextDensity);
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
    MPI_Scatter(init.globalField.nextForceX.data(), size_l, MPI_DOUBLE,
                init.localField.nextForceX.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.nextForceY.data(), size_l, MPI_DOUBLE,
                init.localField.nextForceY.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.nextEntropy.data(), size_l, MPI_DOUBLE,
                init.localField.nextEntropy.data(), size_l, MPI_DOUBLE,
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
    MPI_Scatter(init.globalField.previousEntropy.data(), size_l, MPI_DOUBLE,
                init.localField.previousEntropy.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.previous2Density.data(), size_l, MPI_DOUBLE,
                init.localField.previous2Density.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.previous2VelocityX.data(), size_l, MPI_DOUBLE,
                init.localField.previous2VelocityX.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(init.globalField.previous2VelocityY.data(), size_l, MPI_DOUBLE,
                init.localField.previous2VelocityY.data(), size_l, MPI_DOUBLE,
                0, MPI_COMM_WORLD);


    std::shared_ptr<Lattice> l0 = std::shared_ptr<Lattice>(new Lattice(init.localField.unpackNextDistribution()));
    std::shared_ptr<Lattice> l1 = std::shared_ptr<Lattice>(new Lattice(init.localField.unpackNextDistribution()));

    auto velocityXFFTPtr = std::make_shared<vector<std::complex<double> > >(lengthX_g*(lengthY_g/2 +1));
    auto velocityYFFTPtr = std::make_shared<vector<std::complex<double> > >(lengthX_g*(lengthY_g/2 +1));
    init.forces.setVelocityFFTPtr(velocityXFFTPtr, velocityYFFTPtr);

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
    if(mpi_rank == 0 && startIteration == 0) {
      init.outputs.write(init.globalField, startIteration/writeStep);
    }

    for(int iteration = startIteration+1; iteration <= iterationMax; ++iteration) {

      BOOST_LOG_TRIVIAL(debug) << "Iteration: " << iteration
                               << " - Looping over lattice nodes.";

      auto t0 = std::chrono::high_resolution_clock::now();
      std::swap(l0, l1);

      //calculateMomentsField(*l0, init.localField);
      MPI_Gather(init.localField.previousVelocityX.data(), size_l, MPI_DOUBLE,
                 init.globalField.previousVelocityX.data(), size_l, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
      MPI_Gather(init.localField.previousVelocityY.data(), size_l, MPI_DOUBLE,
                 init.globalField.previousVelocityY.data(), size_l, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

      if(mpi_rank == 0) {
        calculate2DField_FFT(init.globalField.previousVelocityX, *velocityXFFTPtr);
        calculate2DField_FFT(init.globalField.previousVelocityY, *velocityYFFTPtr);
      }

      MPI_Bcast((*velocityXFFTPtr).data(), lengthX_g*(lengthY_g/2 +1), MPI_COMPLEX,
                0, MPI_COMM_WORLD);
      MPI_Bcast((*velocityYFFTPtr).data(), lengthX_g*(lengthY_g/2 +1), MPI_COMPLEX,
                0, MPI_COMM_WORLD);

      push_fusedCollideAndStream(*l0, *l1, init.solver, init.forcing, init.forces,
                                 init.localField, startX, iteration);
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

      if(iteration%writeStep == writeStep-1) {
        BOOST_LOG_TRIVIAL(info) << "Iteration: " << iteration
                                << " - Writing outputs.";
        MPI_Gather(init.localField.previousDensity.data(), size_l, MPI_DOUBLE,
                   init.globalField.previous2Density.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        MPI_Gather(init.localField.previousVelocityX.data(), size_l, MPI_DOUBLE,
                   init.globalField.previous2VelocityX.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        MPI_Gather(init.localField.previousVelocityY.data(), size_l, MPI_DOUBLE,
                   init.globalField.previous2VelocityY.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        }


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
        MPI_Gather(init.localField.nextForceX.data(), size_l, MPI_DOUBLE,
                   init.globalField.nextForceX.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        MPI_Gather(init.localField.nextForceY.data(), size_l, MPI_DOUBLE,
                   init.globalField.nextForceY.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        MPI_Gather(init.localField.nextEntropy.data(), size_l, MPI_DOUBLE,
                   init.globalField.nextEntropy.data(), size_l, MPI_DOUBLE,
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
        MPI_Gather(init.localField.previousEntropy.data(), size_l, MPI_DOUBLE,
                   init.globalField.previousEntropy.data(), size_l, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);

        if(iteration%backupStep == 0) {
          BOOST_LOG_TRIVIAL(info) << "Iteration: " << iteration
                                  << " - Writing backup.";
          init.localField.packNextDistribution(*l1);
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
