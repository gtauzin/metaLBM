#ifndef WRITER_H
#define WRITER_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "Commons.h"
#include "Options.h"
#include "StaticArray.h"
#include "MathVector.h"
#include "Domain.h"
#include "Field.h"
#include "Distribution.h"

namespace lbm {

  template <class T, WriterType writerType>
  class Writer {};

  template <class T>
  class Writer<T, WriterType::Generic> {
  protected:
    const std::string writeFolder;
    const std::string writerFolder;
    const std::string fileExtension;
    const std::string filePrefix;
    const std::string fileFormat;

    const MathVector<int, 3> rankMPI;

    const bool isSerial;
    bool isWritten;

    std::ofstream file;

    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in,
           const std::string& fileFormat_in,
           const MathVector<int, 3> rankMPI_in,
           const bool isSerial_in)
      : writeFolder("../output/")
      , writerFolder(writerFolder_in)
      , fileExtension(fileExtension_in)
      , filePrefix(filePrefix_in)
      , fileFormat(fileFormat_in)
      , rankMPI(rankMPI_in)
      , isSerial(isSerial_in)
      , isWritten(false)
    {}

    std::string getFileName(const unsigned int iteration) {
      std::ostringstream number;

      number << iteration;
      std::string fileName =  writeFolder + writerFolder
        + filePrefix + "-" + number.str() + fileExtension;

      return fileName;
    }

    inline bool getIsSerial() {
      return isSerial;
    }

    inline bool getIsWritten(const unsigned int iteration) {
      return (iteration % writeStep) == 0;
    }

  };

  template <class T>
  class Writer<T, WriterType::VTR>
    : public Writer<T, WriterType::Generic> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, WriterType::Generic>("outputVTR/", filePrefix_in,
                                       ".vtr", "ascii", rankMPI_in, true)
    {}

    using Writer<T, WriterType::Generic>::getIsSerial;
    using Writer<T, WriterType::Generic>::getIsWritten;

    inline void openFile(const unsigned int iteration) {
      std::string fileName = getFileName(iteration);
      file.open(fileName, std::ofstream::out | std::ofstream::trunc);

      if(file) {
        writeHeader();
      }
      else {
        std::cout << "Write VTR: could not open file " << fileName << std::endl;
      }
    }

    inline void closeFile() {
      writeFooter();
      file.close();
    }

    void writeField(const std::string& fieldName,
                    const unsigned int numberComponents,
                    const T * RESTRICT field,
                    const gD& globalDomain) {

      file << "\t\t\t<DataArray type=\"Float32\" "
           << "NumberOfComponents=\"" << numberComponents << "\" "
           << "Name=\"" << fieldName << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++) {
        for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++) {
          for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++) {
            MathVector<unsigned int, 3> iP = {iX, iY, iZ};

            file << "\t\t\t\t";
            for(unsigned int iC = 0; iC < numberComponents; ++iC) {
              file << field[globalDomain.getIndex(iP, iC)] << " ";
            }
             file << std::endl;
          }
        }
      }
      file << "\t\t\t</DataArray>\n";
    }


  private:
    using Writer<T, WriterType::Generic>::file;
    using Writer<T, WriterType::Generic>::fileFormat;
    using Writer<T, WriterType::Generic>::rankMPI;


    std::string getFileName(const unsigned int iteration) {
      if(rankMPI[d::X] == 0) {
        return Writer<T, WriterType::Generic>::getFileName(iteration);
      }

      else {
        return "/dev/null";
      }
    }

    void writeHeader() {
      file << "<?xml version=\"1.0\"?>\n";
      file << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" "
           << "byte_order=\"LittleEndian\">\n";
      file << "<RectilinearGrid WholeExtent=\""
           << gD::start()[d::X] << " " << gD::end()[d::X]-1 << " "
           << gD::start()[d::Y] << " " << gD::end()[d::Y]-1 << " "
           << gD::start()[d::Z] << " " << gD::end()[d::Z]-1 << " \">\n";
      file << "\t<Piece Extent=\""
           << gD::start()[d::X] << " " << gD::end()[d::X]-1 << " "
           << gD::start()[d::Y] << " " << gD::end()[d::Y]-1 << " "
           << gD::start()[d::Z] << " " << gD::end()[d::Z]-1 << " \">\n";
      file << "\t\t<PointData>\n";
    }

    void writeFooter() {
      file << "\t\t</PointData>\n" << "\t\t<CellData />\n" << "\t\t<Coordinates>\n";
      file << "\t\t\t<DataArray type=\"Float32\" Name=\"X\" format=\"ascii\">\n";
      for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++){
        file << "\t\t\t\t" << iX+1 << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t\t<DataArray type=\"Float32\" Name=\"Y\" format=\"ascii\">\n";
      for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++){
        file << "\t\t\t\t" << iY+1 << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t\t<DataArray type=\"Float32\" Name=\"Z\" format=\"ascii\">\n";
      for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++){
        file << "\t\t\t\t" << iZ+1 << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t</Coordinates>\n" << "\t</Piece>\n" << "</RectilinearGrid>\n";
      file << "</VTKFile>\n";
    }

  };

  typedef Writer<dataT, writerT> Writer_;

}

#endif // WRITER_H
