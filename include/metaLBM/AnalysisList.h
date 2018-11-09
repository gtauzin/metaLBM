#ifndef ANALYSISLIST_H
#define ANALYSISLIST_H

#include <fstream>
#include <ostream>
#include <string>

#include "Analysis.h"
#include "Commons.h"
#include "Communication.h"
#include "Computation.h"
#include "FieldList.h"
#include "Helpers.h"
#include "Options.h"
#include "Writer.h"

namespace lbm {

/**
 * AnalysisList containing all analyses.
 *
 * @tparam T datatype.
 * @tparam Architecture on which the code is executed.
 */

template <class T, Architecture architecture>
class ScalarAnalysisList {
 public:
  TotalEnergy<T> totalEnergy;
  TotalEnstrophy<T> totalEnstrophy;
  Communication_& communication;
  ScalarAnalysisWriter_ scalarAnalysisWriter;
  Computation<Architecture::CPU, L::dimD> computationLocal;

  ScalarAnalysisList(FieldList<T, architecture>& fieldList_in,
                     Communication_& communication_in,
                     const unsigned int scalarAnalysisStep_in,
                     const unsigned int startIteration_in)
    : totalEnergy(fieldList_in.density.getData(FFTWInit::numberElements),
                  fieldList_in.velocity.getData(FFTWInit::numberElements))
    , totalEnstrophy(fieldList_in.vorticity.getData(FFTWInit::numberElements))
    , communication(communication_in)
    , scalarAnalysisWriter(prefix, "observables", startIteration_in, scalarAnalysisStep_in)
    , computationLocal(lSD::sStart(), lSD::sEnd())
  {
    if (MPIInit::rank[d::X] == 0) {
      writeAnalysesHeader();
    }
  }

  inline bool getIsAnalyzed(const unsigned int iteration) {
    return scalarAnalysisWriter.getIsAnalyzed(iteration);
  }

  inline void writeAnalyses(const unsigned int iteration) {
    resetAnalyses();

    computationLocal.Do([&] LBM_HOST(const Position& iP) {
      totalEnergy(iP);
      totalEnstrophy(iP);
    });
    computationLocal.synchronize();

    normalizeAnalyses();
    reduceAnalyses();

    if (MPIInit::rank[d::X] == 0) {
      T scalarList[] = {totalEnergy.scalar, totalEnstrophy.scalar};
      scalarAnalysisWriter.openFile(iteration);
      scalarAnalysisWriter.writeAnalysis<T, 2>(iteration, scalarList);
      scalarAnalysisWriter.closeFile();
    }
  }

 private:
  inline void resetAnalyses() {
    totalEnergy.reset();
    totalEnstrophy.reset();
  }

  inline void reduceAnalyses() {
    communication.reduce(&(totalEnergy.scalar), 1);
    communication.reduce(&(totalEnstrophy.scalar), 1);
  }

  inline void normalizeAnalyses() {
    totalEnergy.normalize();
    totalEnstrophy.normalize();
  }

  inline void writeAnalysesHeader() {
    std::string header = "iteration";
    header += " " + std::string(totalEnergy.analysisName);
    header += " " + std::string(totalEnstrophy.analysisName);
    scalarAnalysisWriter.writeHeader(header);
  }
};

template <class T, Architecture architecture>
class SpectralAnalysisList {
 private:
  WaveNumber iK;
  unsigned int kNorm;

 public:
  PowerSpectra<T, gFD::maxWaveNumber()> energySpectra;
  PowerSpectra<T, gFD::maxWaveNumber()> forcingSpectra;
  Communication_& communication;
  SpectralAnalysisWriter_ spectralAnalysisWriter;
  Position offset;
  Computation<Architecture::CPU, L::dimD> computationFourier;

  SpectralAnalysisList(FieldList<T, architecture>& fieldList_in,
                       Communication_& communication_in,
                       const unsigned int spectralAnalysisStep_in,
                       const unsigned int startIteration_in)
    : energySpectra(fieldList_in.velocity.getData(FFTWInit::numberElements),
                    globalLengthPtrdiff_t, "energy_spectra")
    , forcingSpectra(fieldList_in.force.getData(FFTWInit::numberElements),
                     globalLengthPtrdiff_t, "forcing_spectra")
    , communication(communication_in)
    , spectralAnalysisWriter(prefix, "spectra", startIteration_in, spectralAnalysisStep_in)
    , offset(gFD::offset(MPIInit::rank))
    , computationFourier(lFD::start(), lFD::end())
  {
    if (MPIInit::rank[d::X] == 0) {
      writeAnalysesHeader();
    }
  }

  inline bool getIsAnalyzed(const unsigned int iteration) {
    return spectralAnalysisWriter.getIsAnalyzed(iteration);
  }

  inline void writeAnalyses(const unsigned int iteration) {
    resetAnalyses();

    forwardTransformAnalyses();
    computationFourier.Do([&] LBM_HOST(const Position& iFP) {
      auto index = lFD::getIndex(iFP);

      iK[d::X] = iFP[d::X] + offset[d::X] <= gSD::sLength()[d::X] / 2
                     ? iFP[d::X] + offset[d::X]
                     : iFP[d::X] + offset[d::X] - gSD::sLength()[d::X];
      iK[d::Y] = iFP[d::Y] <= gSD::sLength()[d::Y] / 2
                     ? iFP[d::Y]
                     : iFP[d::Y] - gSD::sLength()[d::Y];
      iK[d::Z] = iFP[d::Z] <= gSD::sLength()[d::Z] / 2
                     ? iFP[d::Z]
                     : iFP[d::Z] - gSD::sLength()[d::Z];
      kNorm = iK.norm();

      energySpectra(iFP, index, iK, kNorm);
      forcingSpectra(iFP, index, iK, kNorm);

    });
    computationFourier.synchronize();

    backwardTransformAnalyses();

    normalizeAnalyses();
    reduceAnalyses();

    if(MPIInit::rank[d::X] == 0) {
      T* spectraList[2] = {energySpectra.spectra, forcingSpectra.spectra};
      spectralAnalysisWriter.openFile(iteration);
      spectralAnalysisWriter.writeAnalysis<T, 2, gFD::maxWaveNumber()>(
          iteration, spectraList);
      spectralAnalysisWriter.closeFile();
    }
  }

 private:
  inline void resetAnalyses() {
    energySpectra.reset();
    forcingSpectra.reset();
  }

  inline void forwardTransformAnalyses() {
    energySpectra.executeForward();
    forcingSpectra.executeForward();
  }

  inline void backwardTransformAnalyses() {
    energySpectra.executeBackward();
    forcingSpectra.executeBackward();
  }

  inline void reduceAnalyses() {
    communication.reduce(energySpectra.spectra, gFD::maxWaveNumber());
    communication.reduce(forcingSpectra.spectra, gFD::maxWaveNumber());
  }

  inline void normalizeAnalyses() { energySpectra.normalize(); }

  inline void writeAnalysesHeader() {
    std::string header = "iteration wavenumber";
    header += " " + std::string(energySpectra.analysisName);
    header += " " + std::string(forcingSpectra.analysisName);
    spectralAnalysisWriter.writeHeader(header);
  }
};

 template <class T>
 class PerformanceAnalysisList {
 private:
    double initialMass;
    double currentMass;
    double differenceMass;

    double computationTime;
    double communicationTime;
    double writeFieldTime;
    double writeAnalysisTime;
    double totalTime;
    double mLUPS;
    ScalarAnalysisWriter_ scalarAnalysisWriter;

 public:
   PerformanceAnalysisList(const unsigned int performanceAnalysisStep_in,
                           const unsigned int startIteration_in)
     : initialMass(0.0)
     , currentMass(0.0)
     , differenceMass(0.0)
     , computationTime(0.0)
     , communicationTime(0.0)
     , writeFieldTime(0.0)
     , writeAnalysisTime(0.0)
     , totalTime(0.0)
     , mLUPS(0.0)
     , scalarAnalysisWriter(prefix, "performances", startIteration_in,
                             performanceAnalysisStep_in)
   {
     if (MPIInit::rank[d::X] == 0) {
       writeAnalysesHeader();
     }
   }

    inline void setInitialMass(const double mass) {
      initialMass = mass;
    }

    inline void updateMass(const double mass) {
      currentMass = mass;
      differenceMass = fabs(initialMass - currentMass) / initialMass;
    }

    inline void updateComputationTime(const double computationDuration) {
      computationTime += computationDuration;
      totalTime += computationDuration;
    }

    inline void updateCommunicationTime(const double communicationDuration) {
      communicationTime += communicationDuration;
      totalTime += communicationDuration;
    }

    inline void updateWriteFieldTime(const double writeFieldDuration) {
      writeFieldTime += writeFieldDuration;
      totalTime += writeFieldDuration;
    }

    inline void updateWriteAnalysisTime(const double writeAnalysisDuration) {
      writeAnalysisTime += writeAnalysisDuration;
      totalTime += writeAnalysisDuration;
    }

    inline void updateMLUPS(const double numberIteration) {
      mLUPS = gSD::sVolume() * numberIteration / totalTime /(1e6);
    }

    inline bool getIsAnalyzed(const unsigned int iteration) {
      return scalarAnalysisWriter.getIsAnalyzed(iteration);
    }

    inline void writeAnalysesHeader() {
      std::string header = "iteration computation_time communication_time write_field_time analysis_time total_time MLUPS mass_difference";
      scalarAnalysisWriter.writeHeader(header);
    }

    inline void writeAnalyses(const unsigned int iteration) {
      if (MPIInit::rank[d::X] == 0) {
        double scalarList[] = {computationTime, communicationTime, writeFieldTime,
                               writeAnalysisTime, totalTime, mLUPS, differenceMass};
        scalarAnalysisWriter.openFile(iteration);
        scalarAnalysisWriter.writeAnalysis<T, 7>(iteration, scalarList);
        scalarAnalysisWriter.closeFile();
      }
    }

    inline double getComputationTime() {
      return computationTime;
    }

     inline double getCommunicationTime() {
      return communicationTime;
    }

     inline double getWriteAnalysisTime() {
      return writeAnalysisTime;
    }

     inline double getWriteFieldTime() {
      return writeFieldTime;
    }

    inline double getTotalTime() {
      return totalTime;
    }

    inline double getMLUPS() {
      return mLUPS;
    }

    inline double getInitialMass() {
      return initialMass;
    }

    inline double getFinalMass() {
      return currentMass;
    }

    inline double getDifferenceMass() {
      return differenceMass;
    }

 };

typedef PerformanceAnalysisList<dataT> PerformanceAnalysisList_;

}  // namespace lbm

#endif  // ANALYSISLIST_H
