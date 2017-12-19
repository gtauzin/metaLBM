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

  template <class T, WriterType writerType, WriterFileFormat writerFileFormat>
  class Writer {};

  template <class T>
    class Writer<T, WriterType::Generic, WriterFileFormat::Generic> {
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
    class Writer<T, WriterType::VTR, WriterFileFormat::Generic>
    : public Writer<T, WriterType::Generic, WriterFileFormat::Generic> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in,
           const std::string fileFormat_in)
      : Writer<T, WriterType::Generic, WriterFileFormat::Generic>("outputVTR/", filePrefix_in,
                                                                  ".vtr", fileFormat_in,
                                                                  rankMPI_in, true)
    {}

    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getIsSerial;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getIsWritten;


    inline void closeFile() {
      writeFooter();
      file.close();
    }

    template<unsigned int NumberComponents>
      void writeField(const Field<T, NumberComponents, Architecture::CPU, true> field) {
      INSTRUMENT_ON("Writer<T, WriterType::VTR, WriterFileFormat::ascii>::writeField",2)

      file << "\t\t\t<DataArray type=\"Float32\" "
           << "NumberOfComponents=\"" << NumberComponents << "\" "
           << "Name=\"" << field.fieldName << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++) {
        for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++) {
          for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++) {
            MathVector<unsigned int, 3> iP = {iX, iY, iZ};
            file << "\t\t\t\t";
            for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
              file << field.getGlobalValue(iP, iC) << " ";
            }
             file << std::endl;
          }
        }
      }
      file << "\t\t\t</DataArray>\n";
    }

    template<unsigned int NumberComponents>
    void writeField(const Field<T, NumberComponents, Architecture::CPU, false> field) {
    }


  protected:
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::file;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::fileFormat;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::rankMPI;


    std::string getFileName(const unsigned int iteration) {
      if(rankMPI[d::X] == 0) {
        return Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getFileName(iteration);
      }

      else {
        return "/dev/null";
      }
    }

    void writeHeader() {
      INSTRUMENT_ON("Writer<T, WriterType::VTR, WriterFileFormat::Generic>::writeHeader",2)

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
      INSTRUMENT_ON("Writer<T, WriterType::VTR, WriterFileFormat::Generic>::writeFooter",2)

      file << "\t\t</PointData>\n" << "\t\t<CellData />\n" << "\t\t<Coordinates>\n";
      file << "\t\t\t<DataArray type=\"Float32\" "
           << "Name=\"" << "X" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++){
        file << "\t\t\t\t" << iX+1 << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t\t<DataArray type=\"Float32\" "
           << "Name=\"" << "Y" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++){
        file << "\t\t\t\t" << iY+1 << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t\t<DataArray type=\"Float32\" "
           << "Name=\"" << "Z" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++){
        file << "\t\t\t\t" << iZ+1 << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t</Coordinates>\n" << "\t</Piece>\n" << "</RectilinearGrid>\n";
      file << "</VTKFile>\n";
    }

  };

  template <class T>
    class Writer<T, WriterType::VTR, WriterFileFormat::ascii>
    : public Writer<T, WriterType::VTR, WriterFileFormat::Generic> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, WriterType::VTR, WriterFileFormat::Generic>(filePrefix_in,
                                                              rankMPI_in,
                                                              "ascii")
    {}

    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::getIsSerial;
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::getIsWritten;

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

    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::closeFile;

    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::writeField;


  private:
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::file;
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::fileFormat;
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::rankMPI;

    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::getFileName;
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::writeHeader;
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::writeFooter;
  };


    template <class T>
    class Writer<T, WriterType::VTR, WriterFileFormat::binary>
    : public Writer<T, WriterType::VTR, WriterFileFormat::Generic> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, WriterType::VTR, WriterFileFormat::Generic>(filePrefix_in,
                                                              rankMPI_in,
                                                              "binary")
    {}

    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::getIsSerial;
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::getIsWritten;

    inline void openFile(const unsigned int iteration) {
      std::string fileName = getFileName(iteration);
      file.open(fileName, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
      //file.write();

      if(file) {
        writeHeader();
      }
      else {
        std::cout << "Write VTR: could not open file " << fileName << std::endl;
      }
    }

    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::closeFile;

    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::writeField;


  private:
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::file;
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::fileFormat;
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::rankMPI;

    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::getFileName;
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::writeHeader;
    using Writer<T, WriterType::VTR, WriterFileFormat::Generic>::writeFooter;
  };



  template <class T>
    class Writer<T, WriterType::HDF5, WriterFileFormat::Generic>
    : public Writer<T, WriterType::Generic, WriterFileFormat::Generic> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in,
           const std::string fileFormat_in)
      : Writer<T, WriterType::Generic, WriterFileFormat::Generic>("outputHDF5/",
                                                                  filePrefix_in,
                                                                  ".h5", fileFormat_in,
                                                                  rankMPI_in, false)
    {}

    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getIsSerial;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getIsWritten;


    inline void closeFile() {
      writeFooter();
      file.close();
    }

    template<unsigned int NumberComponents>
      void writeField(const Field<T, NumberComponents, Architecture::CPU, true> field) {

      for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++) {
        for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++) {
          for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++) {
            MathVector<unsigned int, 3> iP = {iX, iY, iZ};
            file << "\t\t\t\t";
            for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
              file << field.getGlobalValue(iP, iC) << " ";
            }
             file << std::endl;
          }
        }
      }

    }

    template<unsigned int NumberComponents>
    void writeField(const Field<T, NumberComponents, Architecture::CPU, false> field) {
    }


  private:
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::file;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::fileFormat;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::rankMPI;


    std::string getFileName(const unsigned int iteration) {
      if(rankMPI[d::X] == 0) {
        return Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getFileName(iteration);
      }

      else {
        return "/dev/null";
      }
    }

    void writeHeader() {
    }

    void writeFooter() {
    }

  };


  typedef Writer<dataT, writerT, writerFileFormatT> Writer_;

}

#endif // WRITER_H
