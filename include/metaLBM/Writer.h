#ifndef WRITER_H
#define WRITER_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#ifdef USE_HDF5
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

  template <class T, InputOutput inputOutput, InputOutputType inputOutputType,
            InputOutputDataFormat inputOutputDataFormat>
  class Writer {};

  template <class T>
  class Writer<T, InputOutput::Generic, InputOutputType::Generic,
               InputOutputDataFormat::Generic> {
  protected:
    const std::string writeFolder;
    const std::string writerFolder;
    const std::string fileExtension;
    const std::string filePrefix;
    const std::string fileFormat;

    const MathVector<int, 3> rankMPI;

    bool isWritten;

    std::ofstream file;

    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in,
           const std::string& fileFormat_in,
           const MathVector<int, 3> rankMPI_in)
      : writeFolder("../output/")
      , writerFolder(writerFolder_in)
      , fileExtension(fileExtension_in)
      , filePrefix(filePrefix_in)
      , fileFormat(fileFormat_in)
      , rankMPI(rankMPI_in)
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

  public:
    inline bool getIsWritten(const unsigned int iteration) {
      return (iteration % writeStep) == 0;
    }

    inline bool getIsBackedUp(const unsigned int iteration) {
      return (iteration % backUpStep) == 0;
    }

  };


  template <class T>
  class Writer<T, InputOutput::Generic, InputOutputType::Generic,
               InputOutputDataFormat::ascii>
    : public Writer<T, InputOutput::Generic, InputOutputType::Generic,
                    InputOutputDataFormat::Generic> {
  public:
    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, InputOutput::Generic, InputOutputType::Generic,
               InputOutputDataFormat::Generic>(writerFolder_in,
                                          filePrefix_in,
                                          fileExtension_in,
                                          "ascii",
                                          rankMPI_in)
    {}

    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getIsWritten;

  protected:
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::file;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::fileFormat;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::rankMPI;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getFileName;

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
  class Writer<T, InputOutput::Generic, InputOutputType::Generic,
               InputOutputDataFormat::binary>
    : public Writer<T, InputOutput::Generic, InputOutputType::Generic,
                    InputOutputDataFormat::Generic> {
  public:
    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, InputOutput::Generic, InputOutputType::Generic,
               InputOutputDataFormat::Generic>(writerFolder_in,
                                          filePrefix_in,
                                          fileExtension_in,
                                          "binary",
                                          rankMPI_in)
    {}

    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getIsWritten;


  protected:
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::file;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::fileFormat;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::rankMPI;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getFileName;

    inline void open(const unsigned int iteration) {
      std::string fileName = getFileName(iteration);
      file.open(fileName, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
    }

    template<class U>
    inline void write(U data) {
      file.write(reinterpret_cast<char*>(&data), sizeof(data));
    }

  };


  template <class T, InputOutputDataFormat inputOutputDataFormat>
    class Writer<T, InputOutput::VTR, InputOutputType::Serial, inputOutputDataFormat>
    : public Writer<T, InputOutput::Generic, InputOutputType::Generic, inputOutputDataFormat> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, InputOutput::Generic, InputOutputType::Generic,
               inputOutputDataFormat>("outputVTR/", filePrefix_in,
                                      ".vtr", rankMPI_in)
    {}

    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::getIsWritten;

    inline bool getIsSerial() {
      return true;
    }

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
      INSTRUMENT_ON("Writer<T, InputOutput::VTR, writerFileFromat>::writeField<NumberComponents>",3)

      file << "<DataArray type=\"Float32\" "
           << "NumberOfComponents=\"" << NumberComponents << "\" "
           << "Name=\"" << field.fieldName << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++) {
        for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++) {
          for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++) {
            MathVector<unsigned int, 3> iP = {iX, iY, iZ};
            for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
              write(field.getGlobalValue(iP, iC));
              file << " ";
            }
            //file << std::endl;
          }
        }
      }
      file << "</DataArray>\n";
    }

    template<unsigned int NumberComponents>
    void writeField(const Field<T, NumberComponents, Architecture::CPU, false> field) {
    }


  private:
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::file;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::fileFormat;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::rankMPI;

    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::open;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::write;


    void writeHeader() {
      INSTRUMENT_ON("Writer<T, InputOutput::VTR, InputOutputDataFormat::Generic>::writeHeader",2)

      file << "<?xml version=\"1.0\"?>\n";
      file << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" "
           << "byte_order=\"LittleEndian\">\n";
      file << "<RectilinearGrid WholeExtent=\""
           << gD::start()[d::X] << " " << gD::end()[d::X]-1 << " "
           << gD::start()[d::Y] << " " << gD::end()[d::Y]-1 << " "
           << gD::start()[d::Z] << " " << gD::end()[d::Z]-1 << " \">\n";
      file << "<Piece Extent=\""
           << gD::start()[d::X] << " " << gD::end()[d::X]-1 << " "
           << gD::start()[d::Y] << " " << gD::end()[d::Y]-1 << " "
           << gD::start()[d::Z] << " " << gD::end()[d::Z]-1 << " \">\n";
      file << "<PointData>\n";
    }

    void writeFooter() {
      INSTRUMENT_ON("Writer<T, InputOutput::VTR, InputOutputDataFormat::Generic>::writeFooter",2)

      file << "</PointData>\n" << "<CellData />\n" << "<Coordinates>\n";

      file << "<DataArray type=\"Float32\" "
           << "Name=\"" << "X" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iX = gD::start()[d::X]; iX < gD::end()[d::X]; iX++){
        write(iX+1);
        file << "\n";
      }
      file << "</DataArray>\n";
      file << "<DataArray type=\"Float32\" "
           << "Name=\"" << "Y" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iY = gD::start()[d::Y]; iY < gD::end()[d::Y]; iY++){
        write(iY+1);
        file << "\n";
      }
      file << "</DataArray>\n";
      file << "<DataArray type=\"Float32\" "
           << "Name=\"" << "Z" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iZ = gD::start()[d::Z]; iZ < gD::end()[d::Z]; iZ++){
        write(iZ+1);
        file << "\n";
      }
      file << "</DataArray>\n";
      file << "</Coordinates>\n" << "</Piece>\n" << "</RectilinearGrid>\n";
      file << "</VTKFile>\n";
    }

  };


  #ifdef USE_HDF5
  template <class T, InputOutputDataFormat inputOutputDataFormat>
  class Writer<T, InputOutput::HDF5, InputOutputType::Serial, inputOutputDataFormat>
    : public Writer<T, InputOutput::Generic, InputOutputType::Serial, inputOutputDataFormat> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, InputOutput::Generic,  InputOutputType::Serial,
               inputOutputDataFormat>("outputHDF5/", filePrefix_in,
                                      ".h5", rankMPI_in)
    {}

    inline bool getIsSerial() {
      return true;
    }

    using Writer<T, InputOutput::Generic, InputOutputType::Serial,
                 inputOutputDataFormat>::getIsWritten;

    inline void openFile(const unsigned int iteration) {
      open(iteration);
    }


    inline void closeFile() {
      statusHDF5 = H5Gclose(groupHDF5);
      statusHDF5 = H5Fclose(fileHDF5);
    }

    template<unsigned int NumberComponents>
    void writeField(Field<T, NumberComponents, Architecture::CPU, true> field) {
      // TODO: Need to create a data compound
      /* hid_t hdf5_pop_type = H5Tcreate (H5T_COMPOUND, sizeof (trivel)); */
      /* H5Tinsert (hdf5_pop_type, "x", HOFFSET (trivel, x), H5_TYPE_IN_USE); */
      /* H5Tinsert (hdf5_pop_type, "y", HOFFSET (trivel, y), H5_TYPE_IN_USE); */
      /* H5Tinsert (hdf5_pop_type, "z", HOFFSET (trivel, z), H5_TYPE_IN_USE); */

      propertyListHDF5 = H5Pcreate(H5P_DATASET_CREATE);
      dataSpaceHDF5 = H5Screate_simple(L::dimD, field.globalArray().size(), NULL);
      dataSetHDF5 = H5Dcreate2(groupHDF5, field.fieldName.c_str(),
                               H5T_NATIVE_DOUBLE, dataSpaceHDF5, H5P_DEFAULT,
                               propertyListHDF5, H5P_DEFAULT);

      statusHDF5 = H5Dwrite(dataSetHDF5, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            field.globalData());

      statusHDF5 = H5Dclose(dataSetHDF5);
      statusHDF5 = H5Sclose(dataSpaceHDF5);
      statusHDF5 = H5Pclose(propertyListHDF5);

    }

    template<unsigned int NumberComponents>
    void writeField(Field<T, NumberComponents, Architecture::CPU, false> field) {
    }


  private:
    using Writer<T, InputOutput::Generic, InputOutputType::Serial,
                 inputOutputDataFormat>::fileFormat;
    using Writer<T, InputOutput::Generic, InputOutputType::Serial,
                 inputOutputDataFormat>::rankMPI;
    using Writer<T, InputOutput::Generic, InputOutputType::Serial,
                 inputOutputDataFormat>::getFileName;

    hid_t fileHDF5;
    hid_t datasetHDF5;
    herr_t statusHDF5;
    hid_t groupHDF5;
    hid_t dataSetHDF5;
    hid_t dataSpaceHDF5;
    hid_t propertyListHDF5;

    inline void open(const unsigned int iteration) {
      std::string fileName = getFileName(iteration);
      fileHDF5 = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      groupHDF5 = H5Gcreate2(fileHDF5, "/field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      if(!fileHDF5) {
        std::cout << "Could not open file " << fileName << std::endl;
      }

    }

  };

#endif // USE_HDF5

  typedef Writer<dataT, inputOutput, inputOutputType, inputOutputDataFormat> Writer_;

}

#endif // WRITER_H
