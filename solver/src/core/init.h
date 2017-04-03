#ifndef INIT_H
#define INIT_H

#include <iostream>
#include <string>
#include <array>
#include <memory>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include <boost/align/aligned_allocator.hpp>
#include <vector>
template<class T, std::size_t Alignment = 1>
#ifdef ALIGN_MEMORY
  using vector = std::vector<T,
  boost::alignment::aligned_allocator<T, Alignment> >;
#else
using vector = std::vector<T>;
#endif

#include "input.h"
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

  void initLogging(const int mpi_rank) {
#ifdef ENABLE_LOG
    std::cout << "Logging enabled\n";
    BOOST_LOG_TRIVIAL(debug) << "Logging for debug starts.";

    std::ostringstream rank;
    rank << mpi_rank;
    logging::add_file_log(//"../log/log_rank-" + rank.str() + ".log");
                          keywords::file_name = "../log/log_rank-" + rank.str() + ".log",
                          keywords::rotation_size = 10 * 1024 * 1024,
                          keywords::time_based_rotation = sinks::file::rotation_at_time_point(0, 0, 0),
                          keywords::format = "[%TimeStamp%]: %Message%"
                          );

    logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::debug);

#else
    logging::core::get()->set_logging_enabled(false);

#endif
  }

  template<class T>
    void initPrint(const int mpi_rank, const int mpi_size, std::string processorName) {

    std::cout << "MPI #" << mpi_rank << " of " << mpi_size
              << " processes running on host " << processorName << "\n";

    if ( mpi_rank == 0 ) {
      std::cout.precision(17);
      std::cout << "-------------------OPTIONS-------------------\n"
                << "MPI #" << mpi_rank
                << " Lattice: D" << L::dimD << "Q" << L::dimQ
                << " Global length Y: " << L::lY_g
                << " Global length X: " << L::lX_g
                << " Global length Y: " << L::lY_g
                << " Global length Z: " << L::lY_g
                /* << " Global memory: " << (int)(s_g*sizeof(valueType)) / (1<<30) << "\n" */
                << "MPI #" << mpi_rank
                << " Local length X : " << L::lX_l
                << " Local length Y : " << L::lY_l
                << " Local length Z : " << L::lZ_l
        //      << " Local memory : " << (int)(s_l*sizeof(valueType)) / (1<<30) << "\n"
                << "---------------------------------------------\n\n"
        //   << "[#MPI%02d] ARCH:     CPUs           \n", mpi_rank)
                << "NPROCS: " << NPROCS << "\n"
                << "NTHREADS: " << NTHREADS << "\n"
                << "STRUCTURE:   SOA \n"
                << "MAX ITERATION: " << iterationMax << "\n\n"
                << "-------------------PARAMETERS----------------\n\n"
                << "tau: " << tau
                << ", correponding viscosity: " << L::cs2 * (tau - 0.5) << "\n"
                << "beta:" << beta
                << ", correponding viscosity: " << L::cs2 * (1/(2 * beta) - 0.5) << "\n";
    }
  }

  template<class T>
    vector<T, CACHE_LINE> generate_initDensity() {
    vector<T, CACHE_LINE> initDensityR(s_g(), initDensityValue);

    switch(initDensityType){
    case InitDensityType::homogeneous:{
      break;
    }
    case InitDensityType::peak:{
      const T densityPeakValue = 3.0 * initDensityValue;
      int centerX = static_cast<int>((L::lX_g-1)*0.4);
      int centerY = static_cast<int>((L::lY_g-1)*0.3);
      int centerZ = static_cast<int>((L::lZ_g-1)*0.2);
      int idx = idx_gF(centerX, centerY, centerZ);
      initDensityR[idx] = densityPeakValue;
      break;
    }
    default:{
      BOOST_LOG_TRIVIAL(error) << "Wrong type of density initialization.";
    }
    }
    return initDensityR;
  }

  template<class T>
    vector<MathVector<T, L::dimD>, CACHE_LINE> generate_initVelocity() {
    vector<MathVector<T, L::dimD>, CACHE_LINE> initVelocityXR(s_g(),
                                                                   MathVector<T, L::dimD>{{initVelocityXValue}});


    switch(initVelocityType){
    case InitVelocityType::homogeneous:{
      break;
    }

    default:{
      BOOST_LOG_TRIVIAL(error) << "Wrong type of velocity initialization.";
    }
    }
    return initVelocityXR;
  }

  template<class T>
    vector<T, CACHE_LINE> generate_initAlpha() {
    vector<T, CACHE_LINE> initAlphaR(s_g(), 2.);

    return initAlphaR;
  }

  template<class T>
    vector<T, CACHE_LINE> generate_initDistributionStart(const vector<T, CACHE_LINE>& density,
                                                         const vector<MathVector<T, L::dimD>, CACHE_LINE>& velocity) {
    vector<T, CACHE_LINE> initDistributionR(L::dimQ*s_g(), 0.0);

    for(int iX = 0; iX < L::lX_g; ++iX) {
      for(int iY = 0; iY < L::lY_g; ++iY) {
        for(int iZ = 0; iZ < L::lZ_g; ++iZ) {
          int idx = idx_gF(iX, iY, iZ);
          T velocity2 = velocity[idx].norm2();
          for(int iQ = 0; iQ < L::dimQ; ++iQ) {
            //T density = density[idx];
            initDistributionR[idxPop_gF(iX, iY, iZ, iQ)] = computeEquilibrium<T>(iQ, density[idx], velocity[idx], velocity2);
          }
        }
      }
    }
    return initDistributionR;
  }

  template<class T>
    vector<T, CACHE_LINE> generate_initDistributionRestart() {
    vector<T, CACHE_LINE> initDistributionR(L::dimQ*s_g(), 0.0);
    std::ostringstream number;
    number << startIteration;

    std::string inputFilename = "../../output/outputBackup/" + std::string(prefix) + "-" + number.str() + ".vtr";
    vector<T, CACHE_LINE> distributionVTK = readVTK<T>(inputFilename, "Distribution");

    for(int iX = 0; iX < L::lX_g; ++iX) {
      for(int iY = 0; iY < L::lY_g; ++iY) {
        for(int iZ = 0; iZ < L::lZ_g; ++iZ) {
          for(int iQ = 0; iQ < L::dimQ; ++iQ) {
            initDistributionR[idxPop_gF(iX, iY, iZ, iQ)] = distributionVTK[L::dimQ*(L::lZ_g*(L::lY_g*iX + iY) + iZ)+ iQ];
          }
        }
      }
    }

    return initDistributionR;
  }

  template<class T>
    vector<T, CACHE_LINE> generate_initDistribution(const vector<T, CACHE_LINE>& density,
                                                    const vector<MathVector<T, L::dimD>, CACHE_LINE>& velocity) {
    if(!startIteration) {
      return generate_initDistributionStart<T>(density, velocity);
    }
    else {
      return generate_initDistributionRestart<T>();
    }
  }

  template<class T>
    std::array<std::shared_ptr<Force<T>>, numberForces> convert_forcesArray() {
    std::array<std::shared_ptr<Force<T>>, numberForces> forcesVectorR;
    for(int k = 0; k < numberForces; ++k) {
      BOOST_LOG_TRIVIAL(info) << "Initializing Force";
      forcesVectorR[k] = Create<T>(forceTypeArray[k],
                                      MathVector<T, 3>({forceAmplitudeXArray[k],
                                            forceAmplitudeYArray[k],
                                            forceAmplitudeZArray[k]}),
                                      MathVector<T, 3>({forceWaveLengthXArray[k],
                                            forceWaveLengthYArray[k],
                                            forceWaveLengthZArray[k]}));
    }
    return forcesVectorR;
  }

  template<class T>
    std::vector<std::shared_ptr<BoundaryCondition<T>> > convert_boundaryConditionsVector() {
    std::vector<std::shared_ptr<BoundaryCondition<T>> > boundaryConditionsVectorR;
    for(int k = 0; k < numberBCs; ++k) {
      BOOST_LOG_TRIVIAL(info) << "Initializing Boundary condition";
      boundaryConditionsVectorR.push_back(Create<T>(boundaryTypeArray[k],
                                                       boundaryPositionArray[k],
                                                       boundaryStartArray[k],
                                                       boundaryEndArray[k],
                                                       boundaryPressureArray[k],
                                                       MathVector<T, L::dimD>{{boundaryVelocityXArray[k],
                                                             boundaryVelocityYArray[k]}}));
    }
    return boundaryConditionsVectorR;
  }

  template<class T>
    std::vector<std::shared_ptr<BoundaryCondition<T>> > localize_BoundaryConditions(const std::vector<std::shared_ptr<BoundaryCondition<T>> >& boundaryConditionsArray_g,
                                                                                       const int  mpi_rank) {
    std::vector<std::shared_ptr<BoundaryCondition<T>> > boundaryConditionsVectorR;

    return boundaryConditionsVectorR;
  }

  template<class T>
    std::array<std::shared_ptr<Output<T> > , numberOutputs> convert_outputsArray() {
    std::array<std::shared_ptr<Output<T> > , numberOutputs> outputsArrayR;
    for(int k = 0; k < numberOutputs; ++k) {
      BOOST_LOG_TRIVIAL(info) << "Initializing Output";
      outputsArrayR[k] = Create<T>(prefix, outputTypeArray[k]);
    }
    return outputsArrayR;
  }

  template<class T>
    class Init {
  public:
    int iterationStart;
    LocalField<T> localField;
    GlobalField<T> globalField;
    std::shared_ptr<Solver<T>> solver;
    std::shared_ptr<ForcingScheme<T>> forcingScheme;
    BoundaryConditions<T> boundaryConditions;
    Forces<T> forces;
    Outputs<T> outputs;

  Init(int iterationStart_in,
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
    Init<T> init_Simulation(const int mpi_rank) {
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

    return Init<T>(iterationStart,
                   localField,
                   globalField,
                   solver,
                   forcingScheme,
                   boundaryConditions,
                   forces,
                   outputs);
  }
}

#endif // INIT_H
