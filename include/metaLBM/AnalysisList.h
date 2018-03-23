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
    ScalarAnalysisWriter_& scalarAnalysisWriter;
    Communication_& communication;
    Computation<Architecture::CPU, L::dimD> computationLocal;

    ScalarAnalysisList(FieldList<T, architecture>& fieldList_in,
                       ScalarAnalysisWriter_& scalarAnalysisWriter_in,
                       Communication_& communication_in)
      : totalEnergy(fieldList_in.density.getMultiData(),
                    fieldList_in.velocity.getMultiData())
      , totalEnstrophy(fieldList_in.vorticity.getMultiData())
      , scalarAnalysisWriter(scalarAnalysisWriter_in)
      , communication(communication_in)
      , computationLocal(lSD::sStart(), lSD::sEnd())
    {
      writeAnalysesHeader();
    }

    HOST
    void operator()(const Position& iP) {
      totalEnergy(iP);
      totalEnstrophy(iP);
    }

    inline void writeAnalyses(const unsigned int iteration) {
      resetAnalyses();
      computationLocal.Do(*this);
      normalizeAnalyses();
      reduceAnalyses();

      T scalarList[] = {totalEnergy.scalar, totalEnstrophy.scalar};

      scalarAnalysisWriter.writeAnalysis<2>(iteration, scalarList);
    }

    inline void resetAnalyses() {
      totalEnergy.reset();
      totalEnstrophy.reset();
    }

    inline void reduceAnalyses() {
      communication.reduce(totalEnergy.scalar);
      communication.reduce(totalEnstrophy.scalar);
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
    SpectralAnalysisWriter_& spectralAnalysisWriter;
    Communication_& communication;
    Position offset;
    Computation<Architecture::CPU, L::dimD> computationFourier;


    SpectralAnalysisList(FieldList<T, architecture>& fieldList_in,
                         SpectralAnalysisWriter_& spectralAnalysisWriter_in,
                         Communication_& communication_in)
      : energySpectra(fieldList_in.velocity.getMultiData())
      , spectralAnalysisWriter(spectralAnalysisWriter_in)
      , communication(communication_in)
      , offset(gFD::offset(communication.getRankMPI()))
      , computationFourier(lFD::start(), lFD::end())
    {
      writeAnalysesHeader();
    }

    HOST
    void operator()(const Position& iFP) {
      auto index = lFD::getIndex(iFP);

      iK[d::X] = iFP[d::X]+offset[d::X] <= gSD::sLength()[d::X]/2 ?
        iFP[d::X]+offset[d::X] : iFP[d::X]+offset[d::X]-gSD::sLength()[d::X];
      iK[d::Y] = iFP[d::Y] <= gSD::sLength()[d::Y]/2 ?
        iFP[d::Y] : iFP[d::Y]-gSD::sLength()[d::Y];
      iK[d::Z] = iFP[d::Z] <= gSD::sLength()[d::Z]/2 ?
        iFP[d::Z] : iFP[d::Z]-gSD::sLength()[d::Z];
      kNorm = iK.norm();

      energySpectra(iFP, index, iK, kNorm);
    }

    //template<unsigned int MaxWaveNumber>
    inline void writeAnalyses(const unsigned int iteration) {
      resetAnalyses();
      computationFourier.Do(*this);
      normalizeAnalyses();
      reduceAnalyses();

      T * spectraList[] = {energySpectra.spectra};

      spectralAnalysisWriter.writeAnalysis<1, gFD::maxWaveNumber()>(iteration, spectraList);
    }

    inline void resetAnalyses() {
      energySpectra.reset();
    }

    inline void reduceAnalyses() {
      communication.reduce(energySpectra.spectra);
    }

    inline void normalizeAnalyses() {
      energySpectra.normalize();
    }

    inline void writeAnalysesHeader() {
      std::string header = "iteration";
      header += " " + std::string(energySpectra.analysisName);
      spectralAnalysisWriter.writeHeader(header);
    }

  };


} // namespace lbm

#endif // ANALYSISLIST_H
