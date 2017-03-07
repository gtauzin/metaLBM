#ifndef COMPUTE_H
#define COMPUTE_H

#include "init.h"
#include "lattice.h"
#include "force.h"

#include <vector>
#include <memory>

namespace lbm {

  std::shared_ptr<Lattice> initLattice(LocalField& field);

  #pragma omp declare simd
  inline void calculateMoments(const double * f, const int idx_lattice,
                               double& density, double& velocityX, double& velocityY) {

    BOOST_LOG_TRIVIAL(debug) << " - Computing density.";
    density = computeDensity(f, idx_lattice);

    BOOST_LOG_TRIVIAL(debug) << " - Computing velocity.";
    velocityX = computeVelocityX(f, idx_lattice, density);
    velocityY = computeVelocityY(f, idx_lattice, density);
  }

  void calculateMomentsField(Lattice& l_previous, LocalField& field);

  void storeNextField(const std::vector<double>& f_next, Field& field, const int startX,
                      std::shared_ptr<Forcing>& forcing, Forces& forces, const int iteration);

  void push_fusedCollideAndStream(Lattice& l_previous, Lattice& l_next,
                                  Solver& solver, std::shared_ptr<Forcing>& forcing, Forces& forces,
                                  LocalField& field, const int startX, const int iteration);

  void compute(Init& init, const int mpi_rank);

}

#endif // COMPUTE_H
