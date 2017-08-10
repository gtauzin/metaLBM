#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <iostream>
#include <string>

#include "Input.h"
#include "Options.h"
#include "Lattice.h"
#include "Domain.h"
#include "MathVector.h"
#include "Field.h"
#include "Equilibrium.h"
#include "Reader.h"
#include "Force.h"
#include "Boundary.h"
#include "Writer.h"

namespace lbm {

  template<class T>
  LocalizedField<T, 1> initGlobalDensity() {
    LocalizedField<T, 1> densityR("density", gD::volume(), initDensityValue);

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
      densityR[gD::getIndex(center)] = densityPeakValue;
      break;
    }
    default: {
      std::cout << "Wrong type of density initialization.";
    }
    }
    return densityR;
  }

  template<class T>
  LocalizedField<T, L::dimD> initGlobalVelocity() {
    MathVector<T, L::dimD> initVelocityValueProjected{{ (T) 0 }};

    initVelocityValueProjected = Project<T, L::dimD>::Do(initVelocityValue);

    LocalizedField<T, L::dimD> velocityR("velocity", gD::volume(),
                                         initVelocityValueProjected);


    switch(initVelocityT){
    case InitVelocityType::Homogeneous: {
      break;
    }

    default:{
      std::cout << "Wrong type of velocity initialization.";
    }
    }
    return velocityR;
  }

  template<class T>
  LocalizedField<T, 1> initGlobalAlpha() {
    LocalizedField<T, 1> alphaR("alpha", gD::volume(), (T) 2);
    return alphaR;
  }

  template<class T>
  LocalizedField<T, L::dimQ> initGlobalDistributionStart(const LocalizedField<T, 1>& globalDensity,
                                                         const LocalizedField<T, L::dimD>& globalVelocity) {
    LocalizedField<T, L::dimQ> distributionR("distribution", gD::volume());

    Equilibrium_ equilibrium;

    for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++) {
      for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++) {
        for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++) {
          unsigned int indexGlobal = gD::getIndex({iX, iY, iZ});

          equilibrium.setVariables(globalDensity.getField(indexGlobal),
                                   globalVelocity.getVector(indexGlobal));

          UnrolledFor<0, L::dimQ>::Do([&] (unsigned int iQ) {
              distributionR[gD::getIndex({iX, iY, iZ}, iQ)]
                = equilibrium.compute(iQ);
          });
        }
      }
    }
    return distributionR;
  }

  template<class T>
  LocalizedField<T, L::dimQ> initGlobalDistributionRestart() {
    Reader<T, L::dimQ, ReaderType::VTR> reader(prefix);
    return reader.readField("distribution", startIteration);
  }

  //initialized un Localized field vide avec un nom et le lire ensuite...

  template<class T>
  LocalizedField<T, L::dimQ> initGlobalDistribution(const LocalizedField<T, 1>& globalDensity,
                                                    const LocalizedField<T, L::dimD>& globalVelocity) {
    if(!startIteration) {
      return initGlobalDistributionStart<T>(globalDensity, globalVelocity);
    }
    else {
      return initGlobalDistributionRestart<T>();
    }
  }

}

#endif // INITIALIZE_H
