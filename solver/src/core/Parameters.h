#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include "Input.h"
#include "structure.h"
#include "commons.h"
#include "solver.h"
#include "field.h"
#include "distribution.h"
#include "boundary.h"
#include "forcingScheme.h"
#include "force.h"
#include "output.h"


namespace lbm {

  typedef lbm::Lattice<lbm::dataType, lbm::latticeType> L;

  template<class T>
  class Parameters {
  public:
    int iterationStart;
    LocalField<T> localField;
    GlobalField<T> globalField;
    std::shared_ptr<Solver<T>> solver;
    std::shared_ptr<ForcingScheme<T>> forcingScheme;
    BoundaryConditions<T> boundaryConditions;
    Forces<T> forces;
    Outputs<T> outputs;

    Parameters(int iterationStart_in,
         LocalField<T> localField_in,
         GlobalField<T> globalField_in,
         std::shared_ptr<Solver<T>> solver_in,
         std::shared_ptr<ForcingScheme<T>> forcingScheme_in,
         BoundaryConditions<T> boundaryConditions_in,
         Forces<T> forces_in,
         Outputs<T> outputs_in)
      : iterationStart(iterationStart_in)
      , localField(localField_in)
      , globalField(globalField_in)
      , solver(solver_in)
      , forcingScheme(forcingScheme_in)
      , boundaryConditions(boundaryConditions_in)
      , forces(forces_in)
      , outputs(outputs_in)
    {}
  };

  template <class T>
  Parameters<T> init_Simulation(const int mpi_rank) {
    constexpr int iterationStart = 1;

    vector<T, CACHE_LINE> initNextDensity = generate_initDensity<T>();
    vector<MathVector<T, L::dimD>, CACHE_LINE> initNextVelocity = generate_initVelocity<T>();
    vector<T, CACHE_LINE> initNextAlpha = generate_initAlpha<T>();

    vector<T, CACHE_LINE> initNextDistribution = generate_initDistribution<T>(initNextDensity,
                                                                              initNextVelocity);

    vector<T, CACHE_LINE> initPreviousDensity = initNextDensity;
    vector<MathVector<T, L::dimD>, CACHE_LINE> initPreviousVelocity = initNextVelocity;

    LocalField<T> localField = LocalField<T>(initNextDensity,
                                             initNextVelocity,
                                             initNextAlpha,
                                             initNextDistribution,
                                             initPreviousDensity,
                                             initPreviousVelocity);


    GlobalField<T> globalField = GlobalField<T>(initNextDensity,
                                                initNextVelocity,
                                                initNextAlpha,
                                                initNextDistribution,
                                                initPreviousDensity,
                                                initPreviousVelocity);

    std::shared_ptr<Solver<T>> solver = Create<T>(solverMethod);

    std::shared_ptr<ForcingScheme<T>> forcingScheme = Create<T>(forcingSchemeMethod);

    std::vector<std::shared_ptr<BoundaryCondition<T>>> boundaryConditionsVector = convert_boundaryConditionsVector<T>();
    BoundaryConditions<T> boundaryConditions = BoundaryConditions<T>(boundaryConditionsVector);

    Forces<T> forces = Forces<T>(convert_forcesArray<T>());

    Outputs<T> outputs = Outputs<T>(convert_outputsArray<T>());

    return Parameters<T>(iterationStart,
                   localField,
                   globalField,
                   solver,
                   forcingScheme,
                   boundaryConditions,
                   forces,
                   outputs);
  }
}

#endif // PARAMETERS_H
