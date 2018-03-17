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
  Field<T, 1, DomainType::LocalSpace,
        architecture, true> initLocalDensity(const unsigned int numberElements,
                                             const MathVector<int, 3>& rankMPI) {
    INSTRUMENT_ON("initLocalDensity<T>",2)

      Field<T, 1, DomainType::LocalSpace,
      architecture, true> densityFieldR("density", numberElements, initDensityValue);

    switch(initDensityT){
    case InitDensityType::Homogeneous: {
      break;
    }
    case InitDensityType::Peak: {
      if(rankMPI[d::X] == 0) {
      const T densityPeakValue = 3.0 * initDensityValue;
      MathVector<unsigned int, 3> center;

      center[d::X] = static_cast<unsigned int>((lSD::sLength()[d::X]-1)* (T) 0.4);
      center[d::Y] = static_cast<unsigned int>((lSD::sLength()[d::Y]-1)* (T) 0.3);
      center[d::Z] = static_cast<unsigned int>((lSD::sLength()[d::Z]-1)* (T) 0.2);

      densityFieldR.setLocalValue(center, densityPeakValue);
      }
      break;
    }
    default: {
      std::cout << "Wrong type of density initialization.";
    }
    }
    return densityFieldR;
  }

  template<class T, Architecture architecture>
  Field<T, L::dimD, DomainType::LocalSpace,
        architecture, true> initLocalVelocity(const unsigned int numberElements) {
    INSTRUMENT_ON("initLocalVelocity<T>",2)

    MathVector<T, L::dimD> initVelocityVectorProjected{{ (T) 0 }};
    initVelocityVectorProjected = Project<T, T, L::dimD>::Do(initVelocityVector);

    Field<T, L::dimD, DomainType::LocalSpace,
      architecture, true> velocityFieldR("velocity", numberElements,
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
        architecture, true> initLocalForce(const unsigned int numberElements,
                                           const MathVector<int, 3>& rankMPI) {
    INSTRUMENT_ON("initLocalForce<T>",2)

      Field<T, L::dimD, DomainType::LocalSpace,
        architecture, true> forceFieldR("force", numberElements, 0);


      if(forceT == ForceType::ConstantShell) {
        Force<T, forceT> force(forceAmplitude, forceWaveLength, forcekMin, forcekMax);

        force.setLocalForceArray(forceFieldR.getMultiData(), gFD::offset(rankMPI)[d::X]);

      }

      else {
        Force<T, forceT> force(forceAmplitude, forceWaveLength, forcekMin, forcekMax);
        MathVector<unsigned int, 3> iP;
        for(unsigned int iZ = lSD::sStart()[d::Z]; iZ < lSD::sEnd()[d::Z]; iZ++) {
          for(unsigned int iY = lSD::sStart()[d::Y]; iY < lSD::sEnd()[d::Y]; iY++) {
            for(unsigned int iX = lSD::sStart()[d::X]; iX < lSD::sEnd()[d::X]; iX++) {
              iP =  MathVector<unsigned int, 3>({iX, iY, iZ});
              force.setForce(iP, gSD::sOffset(rankMPI));
              forceFieldR.setLocalVector(iP, force.getForce());
            }
          }
        }
      }
   return forceFieldR;
  }

  template<class T, Architecture architecture>
  Field<T, 1, DomainType::LocalSpace,
        architecture, true> initLocalAlpha(const unsigned int numberElements) {
    INSTRUMENT_ON("initLocalAlpha<T>",2)

    Field<T, 1, DomainType::LocalSpace,
          architecture, true> alphaFieldR("alpha", numberElements, (T) 2);
    return alphaFieldR;
  }

  template<class T, Architecture architecture>
  MultiDynamicArray<T, Architecture::CPU, L::dimQ> initLocalDistributionStart(const Field<T, 1, DomainType::LocalSpace, architecture, true>& densityField,
                                                                              const Field<T, L::dimD, DomainType::LocalSpace, architecture, true>& velocityField) {
    INSTRUMENT_ON("initLocalDistributionStart<T>",2)

      Field<T, L::dimQ, DomainType::LocalSpace,
            architecture, true> distributionFieldR("distribution",
                                                   densityField.numberElements);

    Equilibrium_ equilibrium;
    MathVector<unsigned int, 3> iP;
    for(unsigned int iZ = lSD::sStart()[d::Z]; iZ < lSD::sEnd()[d::Z]; iZ++) {
      for(unsigned int iY = lSD::sStart()[d::Y]; iY < lSD::sEnd()[d::Y]; iY++) {
        for(unsigned int iX = lSD::sStart()[d::X]; iX < lSD::sEnd()[d::X]; iX++) {
          iP =  MathVector<unsigned int, 3>({iX, iY, iZ});

          equilibrium.setVariables(densityField.getLocalValue(iP),
                                   velocityField.getLocalVector(iP));

          for(unsigned int iQ = 0; iQ < L::dimQ; ++iQ) {
            distributionFieldR.setLocalValue(iP, equilibrium.calculate(iQ), iQ);
          }
        }
      }
    }

    return distributionFieldR.getLocalArray();
  }

  template<class T>
  MultiDynamicArray<T, Architecture::CPU, L::dimQ> initLocalDistributionRestart(const unsigned int numberElements) {
    INSTRUMENT_ON("initLocalDistributionRestart<T>",2)

      Reader<T, L::dimQ, inputOutput, inputOutputType, inputOutputDataFormat> reader(prefix);
    return reader.readArray("distribution", startIteration);
  }


  template<class T, Architecture architecture>
  MultiDynamicArray<T, Architecture::CPU, L::dimQ> initLocalDistribution(const Field<T, 1, DomainType::LocalSpace, architecture, true>& densityField,
                                                             const Field<T, L::dimD, DomainType::LocalSpace, architecture, true>& velocityField) {
    INSTRUMENT_ON("initLocalDistribution<T>",2)

    if(startIteration == 0) {
      return initLocalDistributionStart<T>(densityField, velocityField);
    }
    else {
      return initLocalDistributionRestart<T>(densityField.numberElements);
    }
  }


}

#endif // INITIALIZE_H
