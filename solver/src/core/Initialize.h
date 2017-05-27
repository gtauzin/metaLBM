#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include "Input.h"
#include "Options.h"
#include "Lattice.h"
#include "MathVector.h"
#include "Helpers.h"

namespace lbm {

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

}

#endif // INITIALIZE_H
