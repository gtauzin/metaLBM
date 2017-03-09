#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <memory>
#include <array>

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
#include "lattice.h"

namespace lbm {

  template <class T, LatticeType L>
    class Output {
  protected:
    const std::string outputFolder;
    const std::string filePrefix;
    const std::string fileExtension;

    std::string getFileName(const int fileNumber){
      std::ostringstream number;

      number << fileNumber;
      std::string fileName = "../../output/" + outputFolder
        + filePrefix + number.str() + fileExtension;

      return fileName;
    }

  public:
  Output(const std::string& outputFolder_in,
         const std::string& filePrefix_in,
         const std::string& fileExtension_in)
    : outputFolder(outputFolder_in)
      , filePrefix(filePrefix_in)
      , fileExtension(fileExtension_in)
    {}

    virtual void write(Field<T, L>& field, const int fileNumber) = 0;
  };

  template <class T, LatticeType L>
    class Output_VTR : public Output<T, L> {

  protected:
    void writeHeader(std::ofstream& fileVTR) {
      fileVTR << "<?xml version=\"1.0\"?>\n";
      fileVTR << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
      fileVTR << "<RectilinearGrid WholeExtent=\"0 " << lX_g<T, L>()-1 << " 0 " << lY_g<T, L>()-1 << " 0 " << lZ_g<T, L>()-1 << " 0\">\n";
      fileVTR << "\t<Piece Extent=\"0 "
              << lX_g<T, L>()-1 << " 0 "
              << lY_g<T, L>()-1 << " 0 "
              << lY_g<T, L>()-1 << " 0\">\n";

      // start writing point data
      fileVTR << "\t\t<PointData>\n";
    }

    void writeScalarField(std::ofstream& fileVTR,
                          const std::string& scalarFieldName,
                          const vector<T, CACHE_LINE>& scalarField,
                          std::false_type) {
    }

    void writeScalarField(std::ofstream& fileVTR,
                          const std::string& scalarFieldName,
                          const vector<T, CACHE_LINE>& scalarField,
                          std::true_type) {
      fileVTR << "\t\t\t<DataArray type=\"Float32\" Name=\"" << scalarFieldName << "\" format=\"ascii\">\n";
      for(int iZ = 0; iZ < lZ_g<T, L>(); iZ++) {
        for(int iY = 0; iY < lY_g<T, L>(); iY++) {
          for(int iX = 0; iX < lX_g<T, L>(); iX++) {
            int idx = idx_gF<T, L>(iX, iY, iZ);
            fileVTR << "\t\t\t\t" << scalarField[idx] << "\n";
          }
        }
      }
      fileVTR << "\t\t\t</DataArray>\n";

    }

    void writeVectorField(std::ofstream& fileVTR,
                          const std::string& vectorFieldName,
                          const vector<T, CACHE_LINE>& scalarFieldX,
                          const vector<T, CACHE_LINE>& scalarFieldY,
                          const vector<T, CACHE_LINE>& scalarFieldZ,
                          std::false_type) {
    }

    void writeVectorField(std::ofstream& fileVTR,
                          const std::string& vectorFieldName,
                          const vector<MathVector<T, dimD<T, L>()>, CACHE_LINE>& vectorField,
                          std::true_type) {
      fileVTR << "\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"" << vectorFieldName << "\" format=\"ascii\">\n";
      for(int iZ = 0; iZ < lZ_g<T, L>(); iZ++) {
        for(int iY = 0; iY < lY_g<T, L>(); iY++) {
          for(int iX = 0; iX < lX_g<T, L>(); iX++) {
            int idx = idx_gF<T, L>(iX, iY, iZ);
            fileVTR << "\t\t\t\t" << vectorField[idx][d::X]
                    << " " << vectorField[idx][d::Y] << " "
                    << " " << vectorField[idx][d::Z] << "\n";
          }
        }
      }
      fileVTR << "\t\t\t</DataArray>\n";
    }

    void writeDistribution(std::ofstream& fileVTR,
                           const vector<T, CACHE_LINE>& distribution) {
      fileVTR << "\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"" << dimQ<T, L>()
              << "\" Name=\"Distribution" << "\" format=\"ascii\">" << std::endl;
      for(int iZ = 0; iZ < lZ_g<T, L>(); iZ++) {
        for(int iY = 0; iY < lY_g<T, L>(); iY++) {
          for(int iX = 0; iX < lX_g<T, L>(); iX++) {
            fileVTR << "\t\t\t\t" << distribution[idxPop_gF<T, L>(iX, iY, iZ, 0)] << " ";

            UnrolledFor<1, dimQ<T, L>()-1>::Do([&] (int iQ) {
                fileVTR << distribution[idxPop_gF<T, L>(iX, iY, iZ, iQ)] << " ";
              });

            fileVTR << distribution[idxPop_gF<T, L>(iX, iY, iZ, dimQ<T, L>()-1)] << std::endl;
          }
        }
      }
      fileVTR << "\t\t\t</DataArray>" << std::endl;
    }

    void writeFooter(std::ofstream& fileVTR) {
      fileVTR << "\t\t</PointData>\n";

      fileVTR << "\t\t<CellData />\n";

      fileVTR << "\t\t<Coordinates>\n";
      fileVTR << "\t\t\t<DataArray type=\"Float32\" Name=\"X\" format=\"ascii\">\n";
      for(int iX = 0; iX < lX_g<T, L>(); iX++){
        fileVTR << "\t\t\t\t" << iX+1 << "\n";
      }
      fileVTR << "\t\t\t</DataArray>\n";
      fileVTR << "\t\t\t<DataArray type=\"Float32\" Name=\"Y\" format=\"ascii\">\n";
      for(int iY = 0; iY < lY_g<T,L>(); iY++){
        fileVTR << "\t\t\t\t" << iY+1 << "\n";
      }
      fileVTR << "\t\t\t</DataArray>\n";
      fileVTR << "\t\t\t<DataArray type=\"Float32\" Name=\"Y\" format=\"ascii\">\n";
      for(int iZ = 0; iZ < lZ_g<T,L>(); iZ++){
        fileVTR << "\t\t\t\t" << iZ+1 << "\n";
      }
      fileVTR << "\t\t\t</DataArray>\n";
      fileVTR << "\t\t</Coordinates>\n";

      fileVTR << "\t</Piece>\n";
      fileVTR << "</RectilinearGrid>\n";
      fileVTR << "</VTKFile>\n";
    }

  public:
  Output_VTR(const std::string& filesPrefix)
    : Output<T, L>("outputVTR/", filesPrefix+"-", ".vtr")
      {}

  Output_VTR(const std::string& outputFolder_in,
             const std::string& filesPrefix)
    : Output<T, L>(outputFolder_in, filesPrefix, ".vtr")
      {}

    void write(Field<T, L>& field, const int fileNumber) {
      std::ofstream fileVTR;
      std::string fileName = this->getFileName(fileNumber);
      fileVTR.open(fileName);

      if(fileVTR){
        writeHeader(fileVTR);
        writeScalarField(fileVTR, "Density", field.nextDensity,
                         std::integral_constant<bool, writeNextDensity>());
        writeScalarField(fileVTR, "PreviousDensity", field.previousDensity,
                         std::integral_constant<bool, writePreviousDensity>());
        writeVectorField(fileVTR, "Velocity", field.nextVelocity,
                         std::integral_constant<bool, writeNextVelocity>());
        writeVectorField(fileVTR, "PreviousVelocity", field.previousVelocity,
                         std::integral_constant<bool, writePreviousVelocity>());
        writeScalarField(fileVTR, "Alpha", field.nextAlpha,
                         std::integral_constant<bool, writeNextAlpha>());
        writeFooter(fileVTR);
        fileVTR.close();
      }
      else{
        BOOST_LOG_TRIVIAL(error) << "Output VTR file: could not open " << fileName;
      }
    }

  };

  template <class T, LatticeType L>
    class Output_BackupVTR : public Output_VTR<T, L> {
  public:
  Output_BackupVTR(const std::string& filesPrefix)
    : Output_VTR<T, L>("outputBackup/", filesPrefix+"-")
      {}

    void write(Field<T, L>& field, const int fileNumber) {
      int iteration = fileNumber*writeStep;

      if(iteration%backupStep != 0) {
        return;
      }

      std::ofstream fileVTR;
      std::string fileName = this->getFileName(iteration);
      fileVTR.open(fileName);

      if(fileVTR) {
        this->writeHeader(fileVTR);
        this->writeDistribution(fileVTR, field.nextDistribution);
        this->writeFooter(fileVTR);
        fileVTR.close();
      }
      else {
        BOOST_LOG_TRIVIAL(error) << "Output VTR file: could not open " << fileName;
      }
    }
  };

  template<class T, LatticeType L>
    std::shared_ptr<Output<T, L>> Create(const std::string& filesPrefix,
                                         const OutputType outputType){
    switch(outputType) {
    case OutputType::vtr: {
      return std::shared_ptr<Output<T, L>>(new Output_VTR<T, L>(filesPrefix));
    }
    case OutputType::backup: {
      return std::shared_ptr<Output<T, L>>(new Output_BackupVTR<T, L>(filesPrefix));
    }
    default:{
      BOOST_LOG_TRIVIAL(error) << "Unknown type of output.";
      return nullptr;
    }
    }
  }


  template <class T, LatticeType L>
    class Outputs{
  private:
    std::array<std::shared_ptr<Output<T, L>>, numberOutputs> outputsArray;

  public:
  Outputs(const std::array<std::shared_ptr<Output<T, L>>, numberOutputs>& outputsArray_in)
    : outputsArray(outputsArray_in)
    {}


    inline void write(Field<T, L>& field, const int fileNumber){
      for(auto output : outputsArray){
        output->write(field, fileNumber);
      }
    }

  };

}

#endif // OUTPUT_H
