#ifndef ANALYSISLIST_H
#define ANALYSISLIST_H

#include <string>
#include <fstream>
#include <ostream>

#include "Commons.h"
#include "Options.h"
#include "Helpers.h"
#include "Writer.h"
#include "Analysis.h"
#include "FieldList.h"
#include "Communication.h"
#include "Computation.h"

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
                       Communication_& communication_in)
      : totalEnergy(fieldList_in.density.getLocalData(),
                    fieldList_in.velocity.getLocalData())
      , totalEnstrophy(fieldList_in.vorticity.getLocalData())
      , communication(communication_in)
      , scalarAnalysisWriter(prefix)
      , computationLocal(lSD::sStart(), lSD::sEnd())
    {
      if(communication.rankMPI[d::X] == 0) {
        writeAnalysesHeader();
      }
    }

    inline bool getIsAnalyzed(const unsigned int iteration) {
      return scalarAnalysisWriter.getIsAnalyzed(iteration);
    }

    inline void writeAnalyses(const unsigned int iteration) {
      resetAnalyses();
      computationLocal.Do
        ([&] HOST (const Position& iP) {
          totalEnergy(iP);
          totalEnstrophy(iP);
        });

      normalizeAnalyses();
      reduceAnalyses();

      if(communication.rankMPI[d::X] == 0) {
        T scalarList[] = {totalEnergy.scalar, totalEnstrophy.scalar};
        scalarAnalysisWriter.openFile(iteration);
        scalarAnalysisWriter.writeAnalysis<2>(iteration, scalarList);
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
    EnergySpectra<T, gFD::maxWaveNumber()> energySpectra;
    Communication_& communication;
    SpectralAnalysisWriter_ spectralAnalysisWriter;
    Position offset;
    Computation<Architecture::CPU, L::dimD> computationFourier;


    SpectralAnalysisList(FieldList<T, architecture>& fieldList_in,
                         Communication_& communication_in)
      : energySpectra(fieldList_in.velocity.getLocalData())
      , communication(communication_in)
      , spectralAnalysisWriter(prefix)
      , offset(gFD::offset(communication.rankMPI))
      , computationFourier(lFD::start(), lFD::end())
    {
      if(communication.rankMPI[d::X] == 0) {
        writeAnalysesHeader();
      }
    }

    inline bool getIsAnalyzed(const unsigned int iteration) {
      return spectralAnalysisWriter.getIsAnalyzed(iteration);
    }

    inline void writeAnalyses(const unsigned int iteration) {
      resetAnalyses();
      computationFourier.Do
        ([&] HOST (const Position& iFP) {
          auto index = lFD::getIndex(iFP);

          iK[d::X] = iFP[d::X]+offset[d::X] <= gSD::sLength()[d::X]/2 ?
            iFP[d::X]+offset[d::X] : iFP[d::X]+offset[d::X]-gSD::sLength()[d::X];
          iK[d::Y] = iFP[d::Y] <= gSD::sLength()[d::Y]/2 ?
            iFP[d::Y] : iFP[d::Y]-gSD::sLength()[d::Y];
          iK[d::Z] = iFP[d::Z] <= gSD::sLength()[d::Z]/2 ?
            iFP[d::Z] : iFP[d::Z]-gSD::sLength()[d::Z];
          kNorm = iK.norm();

          energySpectra(iFP, index, iK, kNorm);
        });

      normalizeAnalyses();
      reduceAnalyses();

      if(communication.rankMPI[d::X] == 0) {
        T * spectraList[1] = {energySpectra.spectra};
        spectralAnalysisWriter.openFile(iteration);
        spectralAnalysisWriter.writeAnalysis<1, gFD::maxWaveNumber()>(iteration, spectraList);
        spectralAnalysisWriter.closeFile();
      }
    }

  private:
    inline void resetAnalyses() {
      energySpectra.reset();
    }

    inline void reduceAnalyses() {
      communication.reduce(energySpectra.spectra, gFD::maxWaveNumber());
    }

    inline void normalizeAnalyses() {
      energySpectra.normalize();
    }

    inline void writeAnalysesHeader() {
      std::string header = "iteration wavenumber";
      header += " " + std::string(energySpectra.analysisName);
      spectralAnalysisWriter.writeHeader(header);
    }

  };


} // namespace lbm

#endif // ANALYSISLIST_H
