#ifndef WRITE_H
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
      std::ostringstream number;

      number << iteration;
      std::string fileName =  writeFolder + writerFolder
        + filePrefix + "-" + number.str() + fileExtension;

      return fileName;
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
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getIsBackedUp;

  protected:
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::file;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::fileFormat;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::rankMPI;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getFileName;

    inline void open(const std::string& fileName) {
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
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getIsBackedUp;

  protected:
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::file;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::fileFormat;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::rankMPI;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getFileName;

    inline void open(const std::string& fileName) {
      file.open(fileName, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
    }

    template<class U>
    inline void write(U data) {
      file.write(reinterpret_cast<char*>(&data), sizeof(data));
    }

  };


    template <class T, InputOutputDataFormat inputOutputDataFormat>
  class Writer<T, InputOutput::DAT, InputOutputType::Serial, inputOutputDataFormat>
    : public Writer<T, InputOutput::Generic, InputOutputType::Generic, inputOutputDataFormat> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, InputOutput::Generic, InputOutputType::Generic,
               inputOutputDataFormat>("outputDAT/", filePrefix_in,
                                      ".dat", rankMPI_in)
    {}

    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::getIsWritten;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::getIsBackedUp;

    inline void openFile(const unsigned int iteration) {
      std::string fileName = "/dev/null";

      if(rankMPI[d::X] == 0) {
        fileName = getFileName(iteration);
      }

      open(fileName);
      writeHeader();

    }

    inline void closeFile() {
      file.close();
    }

    template<unsigned int NumberComponents>
    void writeAnalysis(const Field<T, NumberComponents, DomainType::GlobalSpace,
                       Architecture::CPU, true>& field) {
      INSTRUMENT_ON("Writer<T, InputOutput::VTR, writerFileFromat>::writeField<NumberComponents>",3)

      for(unsigned int iZ = gSD::start()[d::Z]; iZ < gSD::end()[d::Z]; iZ++) {
        for(unsigned int iY = gSD::start()[d::Y]; iY < gSD::end()[d::Y]; iY++) {
          for(unsigned int iX = gSD::start()[d::X]; iX < gSD::end()[d::X]; iX++) {
            MathVector<unsigned int, 3> iP = {iX, iY, iZ};
            for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
              write(field.getGlobalValue(iP, iC));
              file << " ";
            }
            //file << std::endl;
          }
        }
      }
    }

    template<unsigned int NumberComponents, DomainType initDomainType>
    void writeField(const Field<T, NumberComponents, initDomainType,
                    Architecture::CPU, false>& field) {
    }


  private:
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::file;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::fileFormat;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::rankMPI;

    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::getFileName;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::open;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::write;


    void writeHeader() {
      INSTRUMENT_ON("Writer<T, InputOutput::VTR, InputOutputDataFormat::Generic>::writeHeader",2)

        file << "<?xml version=\"1.0\"?>\n";
      file << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" "
           << "byte_order=\"LittleEndian\">\n";
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
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::getIsBackedUp;

    inline void openFile(const unsigned int iteration) {
      std::string fileName = "/dev/null";

      if(rankMPI[d::X] == 0) {
        fileName = getFileName(iteration);
      }

      open(fileName);
      writeHeader();

    }

    inline void closeFile() {
      writeFooter();
      file.close();
    }

    template<unsigned int NumberComponents, Architecture architecture>
    void writeField(const Field<T, NumberComponents, DomainType::GlobalSpace,
                    architecture, true>& field) {
      INSTRUMENT_ON("Writer<T, InputOutput::VTR, writerFileFromat>::writeField<NumberComponents>",3)

        file << "<DataArray type=\"Float32\" "
             << "NumberOfComponents=\"" << NumberComponents << "\" "
             << "Name=\"" << field.fieldName << "\" "
             << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iZ = gSD::start()[d::Z]; iZ < gSD::end()[d::Z]; iZ++) {
        for(unsigned int iY = gSD::start()[d::Y]; iY < gSD::end()[d::Y]; iY++) {
          for(unsigned int iX = gSD::start()[d::X]; iX < gSD::end()[d::X]; iX++) {
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

    template<unsigned int NumberComponents, DomainType initDomainType, Architecture architecture>
    void writeField(const Field<T, NumberComponents, initDomainType,
                    architecture, false>& field) {
    }

  private:
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::file;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::fileFormat;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::rankMPI;

    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::getFileName;
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
           << gSD::start()[d::X] << " " << gSD::end()[d::X]-1 << " "
           << gSD::start()[d::Y] << " " << gSD::end()[d::Y]-1 << " "
           << gSD::start()[d::Z] << " " << gSD::end()[d::Z]-1 << " \">\n";
      file << "<Piece Extent=\""
           << gSD::start()[d::X] << " " << gSD::end()[d::X]-1 << " "
           << gSD::start()[d::Y] << " " << gSD::end()[d::Y]-1 << " "
           << gSD::start()[d::Z] << " " << gSD::end()[d::Z]-1 << " \">\n";
      file << "<PointData>\n";
    }

    void writeFooter() {
      INSTRUMENT_ON("Writer<T, InputOutput::VTR, InputOutputDataFormat::Generic>::writeFooter",2)

        file << "</PointData>\n" << "<CellData />\n" << "<Coordinates>\n";

      file << "<DataArray type=\"Float32\" "
           << "Name=\"" << "X" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iX = gSD::start()[d::X]; iX < gSD::end()[d::X]; iX++){
        write(iX+1);
        file << "\n";
      }
      file << "</DataArray>\n";
      file << "<DataArray type=\"Float32\" "
           << "Name=\"" << "Y" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iY = gSD::start()[d::Y]; iY < gSD::end()[d::Y]; iY++){
        write(iY+1);
        file << "\n";
      }
      file << "</DataArray>\n";
      file << "<DataArray type=\"Float32\" "
           << "Name=\"" << "Z" << "\" "
           << "format=\"" + fileFormat + "\">\n";
      for(unsigned int iZ = gSD::start()[d::Z]; iZ < gSD::end()[d::Z]; iZ++){
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
    : public Writer<T, InputOutput::Generic, InputOutputType::Generic,
                    inputOutputDataFormat> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, InputOutput::Generic,  InputOutputType::Generic,
               inputOutputDataFormat>("outputHDF5/", filePrefix_in,
                                      ".h5", rankMPI_in)
      , writerXDMF(filePrefix_in, rankMPI_in)
    {}

    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::getIsWritten;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::getIsBackedUp;

    inline void openFile(const unsigned int iteration) {
      if(rankMPI[d::X] == 0) {
        std::string fileName = getFileName(iteration);
        open(fileName);

        writerXDMF.openFile(iteration);
      }
    }


    inline void closeFile() {
      if(rankMPI[d::X] == 0) {
        // statusHDF5 = H5Gclose(groupHDF5);
        statusHDF5 = H5Fclose(fileHDF5);

        writerXDMF.closeFile();
      }
    }

    template<unsigned int NumberComponents, Architecture architecture>
    void writeField(Field<T, NumberComponents, DomainType::GlobalSpace,
                    architecture, true>& field) {
      INSTRUMENT_ON("Writer<T, InputOutput::HDF5, writerFileFromat>::writeField<NumberComponents>",3)

      if(rankMPI[d::X] == 0) {

        propertyListHDF5 = H5Pcreate(H5P_DATASET_XFER);
        dataSpaceHDF5 = H5Screate_simple(L::dimD,
                                         Project<hsize_t,
                                         unsigned int, L::dimD>::Do(gSD::length()).data(),
                                         NULL);

        for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
          dataSetHDF5 = H5Dcreate2(fileHDF5, (field.fieldName+std::to_string(iC)).c_str(),
                                   H5T_NATIVE_DOUBLE, dataSpaceHDF5, H5P_DEFAULT,
                                   H5P_DEFAULT, H5P_DEFAULT);

          statusHDF5 = H5Dwrite(dataSetHDF5, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                propertyListHDF5,
                                field.getGlobalData()+ iC*gSD::volume());
          statusHDF5 = H5Dclose(dataSetHDF5);
        }

        statusHDF5 = H5Sclose(dataSpaceHDF5);
        statusHDF5 = H5Pclose(propertyListHDF5);

        writerXDMF.writeField(field);
      }
    }

    template<unsigned int NumberComponents, DomainType initDomainType,
             Architecture architecture>
    void writeField(Field<T, NumberComponents, initDomainType,
                    architecture, false>& field) {
    }

  protected:
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::fileFormat;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::rankMPI;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::getFileName;

    hid_t fileHDF5;
    hid_t datasetHDF5;
    herr_t statusHDF5;
    hid_t groupHDF5;
    hid_t dataSetHDF5;
    hid_t dataSpaceHDF5;
    hid_t propertyListHDF5;
    Writer<T, InputOutput::XDMF, InputOutputType::Serial, inputOutputDataFormat> writerXDMF;

  private:
    inline void open(const std::string& fileName) {
      propertyListHDF5 = H5Pcreate(H5P_FILE_ACCESS);
      fileHDF5 = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                           propertyListHDF5);
      H5Pclose(propertyListHDF5);

      // groupHDF5 = H5Gcreate2(fileHDF5, "/field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      if(!fileHDF5) {
        std::cout << "Could not open file " << fileName << std::endl;
      }
    }

  };


  template <class T, InputOutputDataFormat inputOutputDataFormat>
  class Writer<T, InputOutput::HDF5, InputOutputType::Parallel, inputOutputDataFormat>
    : public Writer<T, InputOutput::HDF5, InputOutputType::Serial, inputOutputDataFormat> {
  public:
    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::Writer;

    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::getIsWritten;
    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::getIsBackedUp;

    inline void openFile(const unsigned int iteration) {
      std::string fileName = getFileName(iteration);
      open(fileName);

      if(rankMPI[d::X] == 0) {
        writerXDMF.openFile(iteration);
      }
    }

    inline void closeFile() {
      // statusHDF5 = H5Gclose(groupHDF5);
      statusHDF5 = H5Fclose(fileHDF5);

      if(rankMPI[d::X] == 0) {
        writerXDMF.closeFile();
      }
    }

    template<unsigned int NumberComponents, Architecture architecture>
    void writeField(Field<T, NumberComponents, DomainType::LocalSpace,
                    architecture, true>& field) {
      INSTRUMENT_ON("Writer<T, InputOutput::HDF5, writerFileFromat>::writeField<NumberComponents>",3)

      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        fileSpaceHDF5 = H5Screate_simple(L::dimD,
                                         Project<hsize_t,
                                         unsigned int, L::dimD>::Do(gSD::length()).data(),
                                         NULL);

        dataSetHDF5 = H5Dcreate2(fileHDF5, (field.fieldName+std::to_string(iC)).c_str(),
                                 H5T_NATIVE_DOUBLE, fileSpaceHDF5, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);

        statusHDF5 = H5Sclose(fileSpaceHDF5);

        dataSpaceHDF5 = H5Screate_simple(L::dimD,
                                         Project<hsize_t,
                                         unsigned int, L::dimD>::Do(lSD::length()).data(),
                                         NULL);

        fileSpaceHDF5 = H5Dget_space(dataSetHDF5);

        H5Sselect_hyperslab(fileSpaceHDF5, H5S_SELECT_SET,
                            Project<hsize_t,
                            unsigned int, L::dimD>::Do(gSD::offset(rankMPI)).data(), NULL,
                            Project<hsize_t,
                            unsigned int, L::dimD>::Do(lSD::length()).data(), NULL);

        propertyListHDF5 = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(propertyListHDF5, H5FD_MPIO_COLLECTIVE);

        statusHDF5 = H5Dwrite(dataSetHDF5, H5T_NATIVE_DOUBLE,
                              dataSpaceHDF5, fileSpaceHDF5, propertyListHDF5,
                              field.getLocalData()+ iC*lSD::volume());

        statusHDF5 = H5Dclose(dataSetHDF5);
        statusHDF5 = H5Sclose(dataSpaceHDF5);
        statusHDF5 = H5Sclose(fileSpaceHDF5);
        statusHDF5 = H5Pclose(propertyListHDF5);
      }

      if(rankMPI[d::X] == 0) {
        writerXDMF.writeField(field);
      }

    }

    template<unsigned int NumberComponents, Architecture architecture>
    void writeField(Field<T, NumberComponents, DomainType::LocalSpace,
                    architecture, false>& field) {
    }

    inline void open(const std::string& fileName) {
      propertyListHDF5 = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(propertyListHDF5, MPI_COMM_WORLD, MPI_INFO_NULL);

      fileHDF5 = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                           propertyListHDF5);
      H5Pclose(propertyListHDF5);

      // groupHDF5 = H5Gcreate2(fileHDF5, "/field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      if(!fileHDF5) {
        std::cout << "Could not open file " << fileName << std::endl;
      }
    }

  private:
    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::fileFormat;
    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::rankMPI;
    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::getFileName;

    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::fileHDF5;
    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::statusHDF5;
    //using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
    //inputOutputDataFormat>::groupHDF5;
    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::dataSetHDF5;
    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::dataSpaceHDF5;
    hid_t fileSpaceHDF5;
    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::propertyListHDF5;

    using Writer<T, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>::writerXDMF;

  };


  template <class T, InputOutputDataFormat inputOutputDataFormat>
  class Writer<T, InputOutput::XDMF, InputOutputType::Serial, inputOutputDataFormat>
    : public Writer<T, InputOutput::Generic,
                    InputOutputType::Generic, inputOutputDataFormat> {
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in)
      : Writer<T, InputOutput::Generic, InputOutputType::Generic,
               inputOutputDataFormat>("outputHDF5/", filePrefix_in,
                                      ".xmf", rankMPI_in)
      , fileName("/dev/null")
      , fileNameHDF5("/dev/null")
    {}

    inline void openFile(const unsigned int iteration) {
      fileName = getFileName(iteration);
      fileNameHDF5 = getFileNameHDF5(iteration);
      open(fileName);
      writeHeader();
    }

    inline void closeFile() {
      writeFooter();
      file.close();
    }

    template<unsigned int NumberComponents, DomainType initDomainType,
             Architecture architecture>
    void writeField(const Field<T, NumberComponents, initDomainType,
                    architecture, true>& field) {
      INSTRUMENT_ON("Writer<T, InputOutput::XDMF, writerFileFromat>::writeField<NumberComponents>",3)

        for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
          file << "<Attribute Name=\"" << field.fieldName+std::to_string(iC) << "\" "
               << "AttributeType=\"Scalar\" Center=\"Node\">\n";
          file << "<DataItem Dimensions=\""
               << gSD::length()[d::X];

          for(unsigned int iD = 1; iD < L::dimD; ++iD) {
            file << " " << gSD::length()[iD];
          }
          file << "\" ";
          file << "NumberType=\"Double\" Precision=\"8\" Format=\"HDF\">\n";
          file << fileNameHDF5 << ":/" << field.fieldName+std::to_string(iC) << "\n";
          file << "</DataItem>\n";
          file << "</Attribute>\n";
        }
    }

    template<unsigned int NumberComponents, DomainType initDomainType,
             Architecture architecture>
    void writeField(const Field<T, NumberComponents, initDomainType,
                    architecture, false>& field) {
    }

  private:
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::file;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::filePrefix;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::fileFormat;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::rankMPI;

    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::getFileName;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::open;
    using Writer<T, InputOutput::Generic, InputOutputType::Generic,
                 inputOutputDataFormat>::write;

    std::string fileName;
    std::string fileNameHDF5;

    std::string getFileNameHDF5(const unsigned int iteration) {
      std::ostringstream number;

      number << iteration;
      std::string fileNameR =  filePrefix + "-" + number.str() + ".h5";

      return fileNameR;
    }

    void writeHeader() {
      INSTRUMENT_ON("Writer<T, InputOutput::XDMF, InputOutputDataFormat::Generic>::writeHeader",2)

        file << "<?xml version=\"1.0\"?>\n";
      file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
      file << "<Xdmf>\n";
      file << "<Domain>\n";
      file << "<Grid Name=\"grid\" GridType=\"Uniform\">\n";
      file << "<Topology TopologyType=\"" << L::dimD << "DCoRectMesh\" Dimensions=\""
           << gSD::length()[d::X];

      for(unsigned int iD = 1; iD < L::dimD; ++iD) {
        file << " " << gSD::length()[iD];
      }
      file << "\" />\n";
      file << "<Geometry GeometryType=\"Origin_Dx";
      if(L::dimD >= 2) {
        file << "Dy";
      }
      if(L::dimD == 3) {
        file << "Dz";
      }
      file << "\">\n";
      file << "<DataItem Dimensions=\"" << L::dimD
           << "\" NumberType=\"Integer\" Format=\"XML\">"
           << "0";

      for(unsigned int iD = 1; iD < L::dimD; ++iD) {
        file << " 0";
      }
      file << "</DataItem>\n";

      file << "<DataItem Dimensions=\"" << L::dimD
           << "\" NumberType=\"Integer\" Format=\"XML\">"
           << "1";

      for(unsigned int iD = 1; iD < L::dimD; ++iD) {
        file << " 1";
      }
      file << "</DataItem>\n";
      file << "</Geometry>\n";
    }

    void writeFooter() {
      INSTRUMENT_ON("Writer<T, InputOutput::XDMF, InputOutputDataFormat::Generic>::writeFooter",2)
        file << "</Grid>\n";
      file << "</Domain>\n";
      file << "</Xdmf>\n";
    }

  };

#endif // USE_HDF5

  typedef Writer<dataT, inputOutput, inputOutputType, inputOutputDataFormat> Writer_;

}

#endif // WRITER_H
