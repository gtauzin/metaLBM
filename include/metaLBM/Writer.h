#ifndef WRITER_H
#define WRITER_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#ifdef HDF5_FOUND
  #include <hdf5.h>
#endif

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
      if(rankMPI[d::X] == 0) {
        std::ostringstream number;

        number << iteration;
        std::string fileName =  writeFolder + writerFolder
          + filePrefix + "-" + number.str() + fileExtension;

        return fileName;
      }

      else {
        return "/dev/null";
      }
    }

    inline bool getIsSerial() {
      return isSerial;
    }

    inline bool getIsWritten(const unsigned int iteration) {
      return (iteration % writeStep) == 0;
    }

  };


  template <class T>
  class Writer<T, WriterType::Generic, WriterFileFormat::ascii>
    : public Writer<T, WriterType::Generic, WriterFileFormat::Generic> {
  public:
    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in,
           const MathVector<int, 3> rankMPI_in,
           const bool isSerial_in)
      : Writer<T, WriterType::Generic, WriterFileFormat::Generic>(writerFolder_in,
                                                                  filePrefix_in,
                                                                  fileExtension_in,
                                                                  "ascii",
                                                                  rankMPI_in,
                                                                  isSerial_in)
    {}

    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getIsSerial;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getIsWritten;

  protected:
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::file;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::fileFormat;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::rankMPI;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getFileName;

    inline void open(const unsigned int iteration) {
      std::string fileName = getFileName(iteration);
      file.open(fileName, std::ofstream::out | std::ofstream::trunc);

      if(!file) {
        std::cout << "Could not open file " << fileName << std::endl;
      }

    }

    template<class U>
    inline void write(const U data) {
      file << data;
    }

  };


  template <class T>
  class Writer<T, WriterType::Generic, WriterFileFormat::binary>
    : public Writer<T, WriterType::Generic, WriterFileFormat::Generic> {
  public:
    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in,
           const MathVector<int, 3> rankMPI_in,
           const bool isSerial_in)
      : Writer<T, WriterType::Generic, WriterFileFormat::Generic>(writerFolder_in,
                                                                  filePrefix_in,
                                                                  fileExtension_in,
                                                                  "binary",
                                                                  rankMPI_in,
                                                                  isSerial_in)
    {}

    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getIsSerial;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getIsWritten;


  protected:
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::file;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::fileFormat;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::rankMPI;
    using Writer<T, WriterType::Generic, WriterFileFormat::Generic>::getFileName;

    inline void open(const unsigned int iteration) {
      std::string fileName = getFileName(iteration);
      file.open(fileName, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
    }

    template<class U>
    inline void write(const U data) {
      //file << data;
    }

  };


  template <class T, WriterFileFormat writerFileFormat>
    class Writer<T, WriterType::VTR, writerFileFormat>
    : public Writer<T, WriterType::Generic, writerFileFormat> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, WriterType::Generic, writerFileFormat>("outputVTR/", filePrefix_in,
                                                         ".vtr", rankMPI_in, true)
    {}

    using Writer<T, WriterType::Generic, writerFileFormat>::getIsSerial;
    using Writer<T, WriterType::Generic, writerFileFormat>::getIsWritten;

    inline void openFile(const unsigned int iteration) {
      open(iteration);
      writeHeader();
    }


    inline void closeFile() {
      writeFooter();
      file.close();
    }

    template<unsigned int NumberComponents>
      void writeField(const Field<T, NumberComponents, Architecture::CPU, true> field) {
      INSTRUMENT_ON("Writer<T, WriterType::VTR, writerFileFromat>::writeField<NumberComponents>",3)

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
              write(field.getGlobalValue(iP, iC));
              file << " ";
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


  private:
    using Writer<T, WriterType::Generic, writerFileFormat>::file;
    using Writer<T, WriterType::Generic, writerFileFormat>::fileFormat;
    using Writer<T, WriterType::Generic, writerFileFormat>::rankMPI;

    using Writer<T, WriterType::Generic, writerFileFormat>::open;
    using Writer<T, WriterType::Generic, writerFileFormat>::write;


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
      file << "\t\t</PointData>\n" // << "\t\t<CellData />\n"
           << "\t\t<Coordinates>\n";
      file << "\t\t\t<DataArray type=\"Float32\" "
           << "Name=\"" << "X" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++){
        file << "\t\t\t\t";
        write(iX+1);
        file << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t\t<DataArray type=\"Float32\" "
           << "Name=\"" << "Y" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++){
        file << "\t\t\t\t";
        write(iY+1);
        file << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t\t<DataArray type=\"Float32\" "
           << "Name=\"" << "Z" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++){
        file << "\t\t\t\t";
        write(iZ+1);
        file << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t</Coordinates>\n" << "\t</Piece>\n" << "</RectilinearGrid>\n";
      file << "</VTKFile>\n";
    }

  };


  #ifdef HDF5_FOUND
  template <class T, WriterFileFormat writerFileFormat>
    class Writer<T, WriterType::HDF5, writerFileFormat>
    : public Writer<T, WriterType::Generic, writerFileFormat> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, WriterType::Generic, writerFileFormat>("outputHDF5/", filePrefix_in,
                                                         ".h5", rankMPI_in, true)
    {}

    using Writer<T, WriterType::Generic, writerFileFormat>::getIsSerial;
    using Writer<T, WriterType::Generic, writerFileFormat>::getIsWritten;

    inline void openFile(const unsigned int iteration) {
      open(iteration);
    }


    inline void closeFile() {
      statusHDF5 = H5Fclose(fileHDF5);
    }

    template<unsigned int NumberComponents>
    void writeField(const Field<T, NumberComponents, Architecture::CPU, true> field) {
      dataSpaceHDF5 = H5Screate_simple(L::dimD, NumberComponents * length_g(), NULL);

      propertyListHDF5 = H5Pcreate(H5P_DATASET_CREATE);

      dataSetHDF5 = H5Dcreate(fileHDF5, field.fieldName, H5T_STD_I32BE,
                              dataSpaceHDF5, H5P_DEFAULT, propertyListHDF5, H5P_DEFAULT);

      statusHDF5 = H5Dwrite(dataSetHDF5, H5T_STD_I32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            field.globalData());

      statusHDF5 = H5Dclose(dataSetHDF5);
      statusHDF5 = H5Sclose(dataSpaceHDF5);
      statusHDF5 = H5Pclose(propertyListHDF5);

    }

    template<unsigned int NumberComponents>
    void writeField(const Field<T, NumberComponents, Architecture::CPU, false> field) {
    }


  private:
    using Writer<T, WriterType::Generic, writerFileFormat>::fileFormat;
    using Writer<T, WriterType::Generic, writerFileFormat>::rankMPI;
    using Writer<T, WriterType::Generic, writerFileFormat>::getFileName;

    hid_t fileHDF5;
    hid_t datasetHDF5;
    herr_t statusHDF5;
    hid_t dataSetHDF5;
    hid_t dataSpaceHDF5;
    hid_t propertyListHDF5;

    inline void open(const unsigned int iteration) {
      std::string fileName = getFileName(iteration);
      fileHDF5 = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      if(!fileHDF5) {
        std::cout << "Could not open file " << fileName << std::endl;
      }

    }

  };

#endif // HDF5_FOUND

  typedef Writer<dataT, writerT, writerFileFormatT> Writer_;

}

#endif // WRITER_H
