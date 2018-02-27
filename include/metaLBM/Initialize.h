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
#include "Force.h"
#include "Equilibrium.h"
#include "Reader.h"

namespace lbm {

  template<class T, Architecture architecture>
  Field<T, 1, DomainType::GlobalSpace, architecture, true> initGlobalDensity() {
    INSTRUMENT_ON("initGlobalDensity<T>",2)

      Field<T, 1, DomainType::GlobalSpace, architecture, true> densityFieldR("density", initDensityValue);

    switch(initDensityT){
    case InitDensityType::Homogeneous: {
      break;
    }
    case InitDensityType::Peak: {
      const T densityPeakValue = 3.0 * initDensityValue;
      MathVector<unsigned int, 3> center;

      center[d::X] = static_cast<unsigned int>((gSD::length()[d::X]-1)* (T) 0.4);
      center[d::Y] = static_cast<unsigned int>((gSD::length()[d::Y]-1)* (T) 0.3);
      center[d::Z] = static_cast<unsigned int>((gSD::length()[d::Z]-1)* (T) 0.2);
      densityFieldR.setGlobalValue(gSD::getIndex(center), densityPeakValue);
      break;
    }
    default: {
      std::cout << "Wrong type of density initialization.";
    }
    }
    return densityFieldR;
  }

  template<class T, Architecture architecture>
  Field<T, L::dimD, DomainType::GlobalSpace, architecture, true> initGlobalVelocity() {
    INSTRUMENT_ON("initGlobalVelocity<T>",2)

    MathVector<T, L::dimD> initVelocityVectorProjected{{ (T) 0 }};
    initVelocityVectorProjected = Project<T, T, L::dimD>::Do(initVelocityVector);

    Field<T, L::dimD, DomainType::GlobalSpace, architecture, true> velocityFieldR("velocity",
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

  template<class T, Architecture architecture>
  Field<T, 1, DomainType::GlobalSpace, architecture, true> initGlobalAlpha() {
    INSTRUMENT_ON("initGlobalAlpha<T>",2)

    Field<T, 1, DomainType::GlobalSpace, architecture, true> alphaFieldR("alpha", (T) 2);
    return alphaFieldR;
  }

  template<class T, Architecture architecture>
  DynamicArray<T, Architecture::CPU> initGlobalDistributionStart(const Field<T, 1, DomainType::GlobalSpace, architecture, true>& densityField,
                                                                 const Field<T, L::dimD, DomainType::GlobalSpace, architecture, true>& velocityField) {
    INSTRUMENT_ON("initGlobalDistributionStart<T>",2)

    DynamicArray<T, Architecture::CPU> distributionR(L::dimQ*gSD::volume());

    Equilibrium_ equilibrium;
    MathVector<unsigned int, 3> iP;
    for(unsigned int iZ = gSD::start()[d::Z]; iZ < gSD::end()[d::Z]; iZ++) {
      for(unsigned int iY = gSD::start()[d::Y]; iY < gSD::end()[d::Y]; iY++) {
        for(unsigned int iX = gSD::start()[d::X]; iX < gSD::end()[d::X]; iX++) {
          iP =  MathVector<unsigned int, 3>({iX, iY, iZ});

          equilibrium.setVariables(densityField.getGlobalValue(iP),
                                   velocityField.getGlobalVector(iP));


          for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
            distributionR[gSQD::getIndex(MathVector<unsigned int, 3>({iX, iY, iZ}), iQ)]
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


  template<class T, Architecture architecture>
    DynamicArray<T, Architecture::CPU> initGlobalDistribution(const Field<T, 1, DomainType::GlobalSpace, architecture, true>& densityField,
                                                              const Field<T, L::dimD, DomainType::GlobalSpace, architecture, true>& velocityField) {
    INSTRUMENT_ON("initGlobalDistribution<T>",2)

    if(startIteration == 0) {
      return initGlobalDistributionStart<T>(densityField, velocityField);
    }
    else {
      return initGlobalDistributionRestart<T>();
    }
  }




  template<class T, Architecture architecture>
  Field<T, 1, DomainType::LocalSpace, architecture, true> initLocalDensity() {
    INSTRUMENT_ON("initLocalDensity<T>",2)

      Field<T, 1, DomainType::LocalSpace,
            architecture, true> densityFieldR("density", initDensityValue);

    switch(initDensityT){
    case InitDensityType::Homogeneous: {
      break;
    }
    case InitDensityType::Peak: {
      const T densityPeakValue = 3.0 * initDensityValue;
      MathVector<unsigned int, 3> center;

      center[d::X] = static_cast<unsigned int>((gSD::length()[d::X]-1)* (T) 0.4);
      center[d::Y] = static_cast<unsigned int>((gSD::length()[d::Y]-1)* (T) 0.3);
      center[d::Z] = static_cast<unsigned int>((gSD::length()[d::Z]-1)* (T) 0.2);
      densityFieldR.setLocalValue(gSD::getIndex(center), densityPeakValue);
      break;
    }
    default: {
      std::cout << "Wrong type of density initialization.";
    }
    }
    return densityFieldR;
  }

  template<class T, Architecture architecture>
  Field<T, L::dimD, DomainType::LocalSpace, architecture, true> initLocalVelocity() {
    INSTRUMENT_ON("initLocalVelocity<T>",2)

    MathVector<T, L::dimD> initVelocityVectorProjected{{ (T) 0 }};
    initVelocityVectorProjected = Project<T, T, L::dimD>::Do(initVelocityVector);

    Field<T, L::dimD, DomainType::LocalSpace,
          architecture, true> velocityFieldR("velocity",
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


  template<class T, Architecture architecture>
  Field<T, L::dimD, DomainType::LocalSpace,
        architecture, true> initLocalForce(const MathVector<int, 3>& rankMPI) {
    INSTRUMENT_ON("initLocalForce<T>",2)
    Field<T, L::dimD, DomainType::LocalSpace,
            architecture, true> forceFieldR("force", forceAmplitude);
    Force<T, forceT> force;
    MathVector<unsigned int, 3> iP;
    for(unsigned int iZ = lSD::start()[d::Z]; iZ < lSD::end()[d::Z]; iZ++) {
      for(unsigned int iY = lSD::start()[d::Y]; iY < lSD::end()[d::Y]; iY++) {
        for(unsigned int iX = lSD::start()[d::X]; iX < lSD::end()[d::X]; iX++) {
          iP =  MathVector<unsigned int, 3>({iX, iY, iZ});
          force.setForce(iP, gSD::offset(rankMPI));
          forceFieldR.setLocalVector(lSD::getIndex(iP), force.getForce());
        }
      }
    }

    if(forceT == ForceType::ConstantShell) {


    }

    else {
      std::cout << "Wrong type of force initialization.";
    }
    return forceFieldR;
  }



  template<class T, Architecture architecture>
  Field<T, 1, DomainType::LocalSpace, architecture, true> initLocalAlpha() {
    INSTRUMENT_ON("initLocalAlpha<T>",2)

    Field<T, 1, DomainType::LocalSpace, architecture, true> alphaFieldR("alpha", (T) 2);
    return alphaFieldR;
  }

  template<class T, Architecture architecture>
  DynamicArray<T, Architecture::CPU> initLocalDistributionStart(const Field<T, 1, DomainType::LocalSpace, architecture, true>& densityField,
                                                                const Field<T, L::dimD, DomainType::LocalSpace, architecture, true>& velocityField) {
    INSTRUMENT_ON("initLocalDistributionStart<T>",2)

    DynamicArray<T, Architecture::CPU> distributionR(L::dimQ * lSD::volume());

    Equilibrium_ equilibrium;
    MathVector<unsigned int, 3> iP;
    for(unsigned int iZ = lSD::start()[d::Z]; iZ < lSD::end()[d::Z]; iZ++) {
      for(unsigned int iY = lSD::start()[d::Y]; iY < lSD::end()[d::Y]; iY++) {
        for(unsigned int iX = lSD::start()[d::X]; iX < lSD::end()[d::X]; iX++) {
          iP =  MathVector<unsigned int, 3>({iX, iY, iZ});

          equilibrium.setVariables(densityField.getLocalHostValue(iP),
                                   velocityField.getLocalHostVector(iP));

          for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
            distributionR[lSD::getIndex(iP, iQ)] = equilibrium.calculate(iQ);
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


  template<class T, Architecture architecture>
    DynamicArray<T, Architecture::CPU> initLocalDistribution(const Field<T, 1, DomainType::LocalSpace, architecture, true>& densityField,
                                                             const Field<T, L::dimD, DomainType::LocalSpace, architecture, true>& velocityField) {
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
