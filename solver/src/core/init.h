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
#include "lattice.h"
#include "boundary.h"
#include "force.h"
#include "output.h"


namespace lbm {

  void initLogging(const int mpi_rank) {
#ifdef ENABLE_LOG
    std::cout << "Logging enabled\n";
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

  template<class T, LatticeType L>
    void initPrint(const int mpi_rank, const int mpi_size, std::string processorName) {

    std::cout << "MPI #" << mpi_rank << " of " << mpi_size
              << " processes running on host " << processorName << "\n";

    if ( mpi_rank == 0 ) {
      std::cout.precision(17);
      std::cout << "-------------------OPTIONS-------------------\n"
                << "MPI #" << mpi_rank
                << " Lattice: D" << dimD<T, L>() << "Q" << dimQ<T, L>()
                << " Global length Y: " << lY_g<T, L>()
                << " Global length X: " << lX_g<T, L>()
                << " Global length Y: " << lY_g<T, L>()
                << " Global length Z: " << lY_g<T, L>()
                << " Global memory: " << (T)(s_g<T, L>()*sizeof(T)) / (1<<30) << "\n"
                << "MPI #" << mpi_rank
                << " Local length X : " << lX_l<T, L>()
                << " Local length Y : " << lY_l<T, L>()
                << " Local length Z : " << lZ_l<T, L>()
                << " Local memory : " << (T)(s_l<T, L>()*sizeof(T)) / (1<<30) << "\n"
                << "---------------------------------------------\n\n"
        //   << "[#MPI%02d] ARCH:     CPUs           \n", mpi_rank)
                << "NPROCS: " << NPROCS << "\n"
                << "NTHREADS: " << NTHREADS << "\n"
                << "STRUCTURE:   SOA \n"
                << "MAX ITERATION: " << iterationMax << "\n\n"
                << "-------------------PARAMETERS----------------\n\n"
                << "tau: " << tau
                << ", correponding viscosity: " << cs2<T, L>() * (tau - 0.5) << "\n"
                << "beta:" << beta
                << ", correponding viscosity: " << cs2<T, L>() * (1/(2 * beta) - 0.5) << "\n";
    }
  }

  template<class T, LatticeType L>
    vector<T, CACHE_LINE> generate_initDensity() {
    vector<T, CACHE_LINE> initDensityR(s_g<T, L>(), initDensityValue);

    switch(initDensityType){
    case InitDensityType::homogeneous:{
      break;
    }
    case InitDensityType::peak:{
      const T densityPeakValue = 3.0 * initDensityValue;
      int centerX = static_cast<int>((lX_g<T, L>()-1)*0.4);
      int centerY = static_cast<int>((lY_g<T, L>()-1)*0.3);
      int centerZ = static_cast<int>((lZ_g<T, L>()-1)*0.2);
      int idx = idx_gF<T, L>(centerX, centerY, centerZ);
      initDensityR[idx] = densityPeakValue;
      break;
    }
    default:{
      BOOST_LOG_TRIVIAL(error) << "Wrong type of density initialization.";
    }
    }
    return initDensityR;
  }

  template<class T, LatticeType L>
    vector<MathVector<T, dimD<T, L>()>, CACHE_LINE> generate_initVelocity() {
    vector<MathVector<T, dimD<T, L>()>, CACHE_LINE> initVelocityXR(s_g<T, L>(),
                                                                   MathVector<T, dimD<T, L>()>{{initVelocityXValue}});


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

  template<class T, LatticeType L>
    vector<T, CACHE_LINE> generate_initAlpha() {
    vector<T, CACHE_LINE> initAlphaR(s_g<T, L>(), 2.);

    return initAlphaR;
  }

  template<class T, LatticeType L>
    vector<T, CACHE_LINE> generate_initDistributionStart(const vector<T, CACHE_LINE>& density,
                                                         const vector<MathVector<T, dimD<T, L>()>, CACHE_LINE>& velocity) {
    vector<T, CACHE_LINE> initDistributionR(dimQ<T, L>()*s_g<T, L>(), 0.0);

    for(int iX = 0; iX < lX_g<T, L>(); ++iX) {
      for(int iY = 0; iY < lY_g<T, L>(); ++iY) {
        for(int iZ = 0; iZ < lZ_g<T, L>(); ++iZ) {
          int idx = idx_gF<T, L>(iX, iY, iZ);
          T velocity2 = velocity[idx].norm2();
          for(int iQ = 0; iQ < dimQ<T, L>(); ++iQ) {
            //T density = density[idx];
            initDistributionR[idxPop_gF<T, L>(iX, iY, iZ, iQ)] = computeEquilibrium<T, L>(iQ, density[idx], velocity[idx], velocity2);
          }
        }
      }
    }
    return initDistributionR;
  }

  template<class T, LatticeType L>
    vector<T, CACHE_LINE> generate_initDistributionRestart() {
    vector<T, CACHE_LINE> initDistributionR(dimQ<T, L>()*s_g<T, L>(), 0.0);
    std::ostringstream number;
    number << startIteration;

    std::string inputFilename = "../../output/outputBackup/" + std::string(prefix) + "-" + number.str() + ".vtr";
    vector<T, CACHE_LINE> distributionVTK = readVTK<T, L>(inputFilename, "Distribution");

    for(int iX = 0; iX < lX_g<T, L>(); ++iX) {
      for(int iY = 0; iY < lY_g<T, L>(); ++iY) {
        for(int iZ = 0; iZ < lZ_g<T, L>(); ++iZ) {
          for(int iQ = 0; iQ < dimQ<T, L>(); ++iQ) {
            initDistributionR[idxPop_gF<T, L>(iX, iY, iZ, iQ)] = distributionVTK[dimQ<T, L>()*(lZ_g<T, L>()*(lY_g<T, L>()*iX + iY) + iZ)+ iQ];
          }
        }
      }
    }

    return initDistributionR;
  }

  template<class T, LatticeType L>
    vector<T, CACHE_LINE> generate_initDistribution(const vector<T, CACHE_LINE>& density,
                                                    const vector<MathVector<T, dimD<T, L>()>, CACHE_LINE>& velocity) {
    if(!startIteration) {
      return generate_initDistributionStart<T, L>(density, velocity);
    }
    else {
      return generate_initDistributionRestart<T, L>();
    }
  }

  template<class T, LatticeType L>
    std::array<std::shared_ptr<Force<T, L>>, numberForces> convert_forcesArray() {
    std::array<std::shared_ptr<Force<T, L>>, numberForces> forcesVectorR;
    for(int k = 0; k < numberForces; ++k) {
      BOOST_LOG_TRIVIAL(info) << "Initializing Force";
      forcesVectorR[k] = Create<T, L>(forceTypeArray[k],
                                      MathVector<T, dimD<T, L>()>({forceAmplitudeXArray[k],
                                            forceAmplitudeYArray[k]}),
                                      MathVector<T, dimD<T, L>()>({forceWaveLengthXArray[k],
                                            forceWaveLengthYArray[k]}));
    }
    return forcesVectorR;
  }

  template<class T, LatticeType L>
    std::vector<std::shared_ptr<BoundaryCondition<T, L>> > convert_boundaryConditionsVector() {
    std::vector<std::shared_ptr<BoundaryCondition<T, L>> > boundaryConditionsVectorR;
    for(int k = 0; k < numberBCs; ++k) {
      BOOST_LOG_TRIVIAL(info) << "Initializing Boundary condition";
      boundaryConditionsVectorR.push_back(Create<T, L>(boundaryTypeArray[k],
                                                       boundaryPositionArray[k],
                                                       boundaryStartArray[k],
                                                       boundaryEndArray[k],
                                                       boundaryPressureArray[k],
                                                       MathVector<T, dimD<T, L>()>{{boundaryVelocityXArray[k],
                                                             boundaryVelocityYArray[k]}}));
    }
    return boundaryConditionsVectorR;
  }

  template<class T, LatticeType L>
    std::vector<std::shared_ptr<BoundaryCondition<T, L>> > localize_BoundaryConditions(const std::vector<std::shared_ptr<BoundaryCondition<T, L>> >& boundaryConditionsArray_g,
                                                                                       const int  mpi_rank) {
    std::vector<std::shared_ptr<BoundaryCondition<T, L>> > boundaryConditionsVectorR;

    return boundaryConditionsVectorR;
  }

  template<class T, LatticeType L>
    std::array<std::shared_ptr<Output<T, L> > , numberOutputs> convert_outputsArray() {
    std::array<std::shared_ptr<Output<T, L> > , numberOutputs> outputsArrayR;
    for(int k = 0; k < numberOutputs; ++k) {
      BOOST_LOG_TRIVIAL(info) << "Initializing Output";
      outputsArrayR[k] = Create<T, L>(prefix, outputTypeArray[k]);
    }
    return outputsArrayR;
  }

  template<class T, LatticeType L>
    class Init {
  public:
    int iterationStart;
    LocalField<T, L> localField;
    GlobalField<T, L> globalField;
    std::shared_ptr<Solver<T, L>> solver;
    std::shared_ptr<Forcing<T, L>> forcing;
    BoundaryConditions<T, L> boundaryConditions;
    Forces<T, L> forces;
    Outputs<T, L> outputs;

  Init(int iterationStart_in,
       LocalField<T, L> localField_in,
       GlobalField<T, L> globalField_in,
       std::shared_ptr<Solver<T, L>> solver_in,
       std::shared_ptr<Forcing<T, L>> forcing_in,
       BoundaryConditions<T, L> boundaryConditions_in,
       Forces<T, L> forces_in,
       Outputs<T, L> outputs_in)
    : iterationStart(iterationStart_in)
      , localField(localField_in)
      , globalField(globalField_in)
      , solver(solver_in)
      , forcing(forcing_in)
      , boundaryConditions(boundaryConditions_in)
      , forces(forces_in)
      , outputs(outputs_in)
    {}
  };

  template <class T, LatticeType L>
    Init<T, L> init_Simulation(const int mpi_rank) {
    constexpr int iterationStart = 1;

    LocalField<T, L> localField;

    vector<T, CACHE_LINE> initNextDensity = generate_initDensity<T, L>();
    vector<MathVector<T, dimD<T, L>()>, CACHE_LINE> initNextVelocity = generate_initVelocity<T, L>();
    vector<T, CACHE_LINE> initNextAlpha = generate_initAlpha<T, L>();


    vector<T, CACHE_LINE> initNextDistribution = generate_initDistribution<T, L>(initNextDensity,
                                                                                 initNextVelocity);

    vector<T, CACHE_LINE> initPreviousDensity = initNextDensity;
    vector<MathVector<T, dimD<T, L>()>, CACHE_LINE> initPreviousVelocity = initNextVelocity;

    GlobalField<T, L> globalField = GlobalField<T, L>(initNextDensity,
                                                      initNextVelocity,
                                                      initNextAlpha,
                                                      initNextDistribution,
                                                      initPreviousDensity,
                                                      initPreviousVelocity);

    std::shared_ptr<Solver<T, L>> solver = Create<T, L>(solverMethod);

    std::shared_ptr<Forcing<T, L>> forcing = Create<T, L>(forcingMethod);

    std::vector<std::shared_ptr<BoundaryCondition<T, L>>> boundaryConditionsVector = convert_boundaryConditionsVector<T, L>();
    BoundaryConditions<T, L> boundaryConditions = BoundaryConditions<T, L>(localize_BoundaryConditions<T, L>(boundaryConditionsVector, mpi_rank));

    Forces<T, L> forces = Forces<T, L>(convert_forcesArray<T, L>());

    Outputs<T, L> outputs = Outputs<T, L>(convert_outputsArray<T, L>());

    return Init<T, L>(iterationStart,
                      localField,
                      globalField,
                      solver,
                      forcing,
                      boundaryConditions,
                      forces,
                      outputs);
  }
}

#endif // INIT_H
