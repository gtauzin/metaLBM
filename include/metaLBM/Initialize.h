#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <iostream>
#include <string>

#include "Commons.h"
#include "Options.h"
#include "Lattice.h"
#include "Domain.h"
#include "DynamicArray.cuh"
#include "MathVector.h"
#include "Field.h"
#include "Equilibrium.h"
#include "Reader.h"

namespace lbm {

  template<class T>
  Field<T, 1, DomainType::Global, Architecture::CPU, true> initGlobalDensity() {
    INSTRUMENT_ON("initGlobalDensity<T>",2)

      Field<T, 1, DomainType::Global, Architecture::CPU, true> densityFieldR("density", initDensityValue);

    switch(initDensityT){
    case InitDensityType::Homogeneous: {
      break;
    }
    case InitDensityType::Peak: {
      const T densityPeakValue = 3.0 * initDensityValue;
      MathVector<unsigned int, 3> center;

      center[d::X] = static_cast<unsigned int>((gD::length()[d::X]-1)* (T) 0.4);
      center[d::Y] = static_cast<unsigned int>((gD::length()[d::Y]-1)* (T) 0.3);
      center[d::Z] = static_cast<unsigned int>((gD::length()[d::Z]-1)* (T) 0.2);
      densityFieldR.setGlobalValue(gD::getIndex(center), densityPeakValue);
      break;
    }
    default: {
      std::cout << "Wrong type of density initialization.";
    }
    }
    return densityFieldR;
  }

  template<class T>
  Field<T, L::dimD, DomainType::Global, Architecture::CPU, true> initGlobalVelocity() {
    INSTRUMENT_ON("initGlobalVelocity<T>",2)

    MathVector<T, L::dimD> initVelocityVectorProjected{{ (T) 0 }};
    initVelocityVectorProjected = Project<T, T, L::dimD>::Do(initVelocityVector);

    Field<T, L::dimD, DomainType::Global, Architecture::CPU, true> velocityFieldR("velocity",
                                           initVelocityVectorProjected);

    switch(initVelocityT){
    case InitVelocityType::Homogeneous: {
      break;
    }

    default:{
      std::cout << "Wrong type of velocity initialization.";
    }
    }
    return velocityFieldR;
  }

  template<class T>
  Field<T, 1, DomainType::Global, Architecture::CPU, true> initGlobalAlpha() {
    INSTRUMENT_ON("initGlobalAlpha<T>",2)

    Field<T, 1, DomainType::Global, Architecture::CPU, true> alphaFieldR("alpha", (T) 2);
    return alphaFieldR;
  }

  template<class T>
  DynamicArray<T, Architecture::CPU> initGlobalDistributionStart(const Field<T, 1, DomainType::Global, Architecture::CPU, true>& densityField,
                                                                 const Field<T, L::dimD, DomainType::Global, Architecture::CPU, true>& velocityField) {
    INSTRUMENT_ON("initGlobalDistributionStart<T>",2)

    DynamicArray<T, Architecture::CPU> distributionR(L::dimQ*gD::volume());

    Equilibrium_ equilibrium;
    MathVector<unsigned int, 3> iP;
    for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++) {
      for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++) {
        for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++) {
          iP =  MathVector<unsigned int, 3>({iX, iY, iZ});

          equilibrium.setVariables(densityField.getGlobalValue(iP),
                                   velocityField.getGlobalVector(iP));


          for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
            distributionR[gQD::getIndex(MathVector<unsigned int, 3>({iX, iY, iZ}), iQ)]
              = equilibrium.calculate(iQ);
          }
        }
      }
    }
    return distributionR;
  }

  template<class T>
  DynamicArray<T, Architecture::CPU> initGlobalDistributionRestart() {
    INSTRUMENT_ON("initGlobalDistributionRestart<T>",2)

      Reader<T, L::dimQ, inputOutput, inputOutputType, inputOutputDataFormat> reader(prefix);
    return reader.readArray("distribution", startIteration);
  }


  template<class T>
    DynamicArray<T, Architecture::CPU> initGlobalDistribution(const Field<T, 1, DomainType::Global, Architecture::CPU, true>& densityField,
                                                              const Field<T, L::dimD, DomainType::Global, Architecture::CPU, true>& velocityField) {
    INSTRUMENT_ON("initGlobalDistribution<T>",2)

    if(startIteration == 0) {
      return initGlobalDistributionStart<T>(densityField, velocityField);
    }
    else {
      return initGlobalDistributionRestart<T>();
    }
  }




  template<class T>
  Field<T, 1, DomainType::Local, Architecture::CPU, true> initLocalDensity() {
    INSTRUMENT_ON("initLocalDensity<T>",2)

      Field<T, 1, DomainType::Local,
            Architecture::CPU, true> densityFieldR("density", initDensityValue);

    switch(initDensityT){
    case InitDensityType::Homogeneous: {
      break;
    }
    case InitDensityType::Peak: {
      const T densityPeakValue = 3.0 * initDensityValue;
      MathVector<unsigned int, 3> center;

      center[d::X] = static_cast<unsigned int>((gD::length()[d::X]-1)* (T) 0.4);
      center[d::Y] = static_cast<unsigned int>((gD::length()[d::Y]-1)* (T) 0.3);
      center[d::Z] = static_cast<unsigned int>((gD::length()[d::Z]-1)* (T) 0.2);
      densityFieldR.setLocalValue(gD::getIndex(center), densityPeakValue);
      break;
    }
    default: {
      std::cout << "Wrong type of density initialization.";
    }
    }
    return densityFieldR;
  }

  template<class T>
  Field<T, L::dimD, DomainType::Local, Architecture::CPU, true> initLocalVelocity() {
    INSTRUMENT_ON("initLocalVelocity<T>",2)

    MathVector<T, L::dimD> initVelocityVectorProjected{{ (T) 0 }};
    initVelocityVectorProjected = Project<T, T, L::dimD>::Do(initVelocityVector);

    Field<T, L::dimD, DomainType::Local,
          Architecture::CPU, true> velocityFieldR("velocity",
                                                  initVelocityVectorProjected);

    switch(initVelocityT){
    case InitVelocityType::Homogeneous: {
      break;
    }

    default:{
      std::cout << "Wrong type of velocity initialization.";
    }
    }
    return velocityFieldR;
  }

  template<class T>
  Field<T, 1, DomainType::Local, Architecture::CPU, true> initLocalAlpha() {
    INSTRUMENT_ON("initLocalAlpha<T>",2)

    Field<T, 1, DomainType::Local, Architecture::CPU, true> alphaFieldR("alpha", (T) 2);
    return alphaFieldR;
  }

  template<class T>
  DynamicArray<T, Architecture::CPU> initLocalDistributionStart(const Field<T, 1, DomainType::Local, Architecture::CPU, true>& densityField,
                                                                const Field<T, L::dimD, DomainType::Local, Architecture::CPU, true>& velocityField) {
    INSTRUMENT_ON("initLocalDistributionStart<T>",2)

    DynamicArray<T, Architecture::CPU> distributionR(L::dimQ * lD::volume());

    Equilibrium_ equilibrium;
    MathVector<unsigned int, 3> iP;
    for(unsigned int iZ = lD::start()[d::Z]; iZ < lD::end()[d::Z]; iZ++) {
      for(unsigned int iY = lD::start()[d::Y]; iY < lD::end()[d::Y]; iY++) {
        for(unsigned int iX = lD::start()[d::X]; iX < lD::end()[d::X]; iX++) {
          iP =  MathVector<unsigned int, 3>({iX, iY, iZ});

          equilibrium.setVariables(densityField.getLocalHostValue(iP),
                                   velocityField.getLocalHostVector(iP));

          for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
            distributionR[lD::getIndex(iP, iQ)] = equilibrium.calculate(iQ);
          }
        }
      }
    }
    return distributionR;
  }

  template<class T>
  DynamicArray<T, Architecture::CPU> initLocalDistributionRestart() {
    INSTRUMENT_ON("initLocalDistributionRestart<T>",2)

      Reader<T, L::dimQ, inputOutput, inputOutputType, inputOutputDataFormat> reader(prefix);
    return reader.readArray("distribution", startIteration);
  }


  template<class T>
    DynamicArray<T, Architecture::CPU> initLocalDistribution(const Field<T, 1, DomainType::Local, Architecture::CPU, true>& densityField,
                                                             const Field<T, L::dimD, DomainType::Local, Architecture::CPU, true>& velocityField) {
    INSTRUMENT_ON("initLocalDistribution<T>",2)

    if(startIteration == 0) {
      return initLocalDistributionStart<T>(densityField, velocityField);
    }
    else {
      return initLocalDistributionRestart<T>();
    }
  }


}

#endif // INITIALIZE_H
