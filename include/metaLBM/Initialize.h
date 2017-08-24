#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <iostream>
#include <string>

#include "Commons.h"
#include "Options.h"
#include "Lattice.h"
#include "Domain.h"
#include "MathVector.h"
#include "Field.h"
#include "Distribution.h"
#include "Equilibrium.h"
#include "Reader.h"

namespace lbm {

  template<class T>
  LocalizedField<T, 1> initGlobalDensity() {
    SCOREP_INSTRUMENT_ON("initGlobalDensity<T>")

    Field<T, 1, true> densityFieldR("density", initDensityValue);

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
    return densityFieldR.getGlobalField();
  }

  template<class T>
  LocalizedField<T, L::dimD> initGlobalVelocity() {
    SCOREP_INSTRUMENT_ON("initGlobalVelocity<T>")

    MathVector<T, L::dimD> initVelocityVectorProjected{{ (T) 0 }};
    initVelocityVectorProjected = Project<T, L::dimD>::Do(initVelocityVector);

    Field<T, L::dimD, true> velocityFieldR("velocity",
                                           initVelocityVectorProjected);


    switch(initVelocityT){
    case InitVelocityType::Homogeneous: {
      break;
    }

    default:{
      std::cout << "Wrong type of velocity initialization.";
    }
    }
    return velocityFieldR.getGlobalField();
  }

  template<class T>
  LocalizedField<T, 1> initGlobalAlpha() {
    SCOREP_INSTRUMENT_ON("initGlobalAlpha<T>")

    Field<T, 1, true> alphaFieldR("alpha", (T) 2);
    return alphaFieldR.getGlobalField();
  }

  template<class T>
  LocalizedField<T, L::dimQ> initGlobalDistributionStart(const Field<T, 1, true>& densityField,
                                                         const Field<T, L::dimD, true>& velocityField) {
    SCOREP_INSTRUMENT_ON("initGlobalDistributionStart<T>")

    LocalizedField<T, L::dimQ> distributionR("distribution", gD::volume());

    Equilibrium_ equilibrium;
    MathVector<unsigned int, 3> iP;
    for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++) {
      for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++) {
        for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++) {
          iP = {iX, iY, iZ};

          equilibrium.setVariables(densityField.getGlobalValue(iP),
                                   velocityField.getGlobalVector(iP));


          UnrolledFor<0, L::dimQ>::Do([&] (unsigned int iQ) {
              distributionR[gQD::getIndex({iX, iY, iZ}, iQ)]
                = equilibrium.calculate(iQ);
            });
        }
      }
    }
    return distributionR;
  }

  template<class T>
  LocalizedField<T, L::dimQ> initGlobalDistributionRestart() {
    SCOREP_INSTRUMENT_ON("initGlobalDistributionRestart<T>")

    Reader<T, L::dimQ, ReaderType::VTR> reader(prefix);
    return reader.readLocalizedField("distribution", startIteration);
  }


  template<class T>
  LocalizedField<T, L::dimQ> initGlobalDistribution(const Field<T, 1, true>& densityField,
                                                    const Field<T, L::dimD, true>& velocityField) {
    SCOREP_INSTRUMENT_ON("initGlobalDistribution<T>")

    if(startIteration == 0) {
      return initGlobalDistributionStart<T>(densityField, velocityField);
    }
    else {
      return initGlobalDistributionRestart<T>();
    }
  }

}

#endif // INITIALIZE_H
