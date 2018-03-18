#ifndef WRITER_H
#define WRITER_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <hdf5.h>

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
    const MathVector<int, 3> sizeMPI;

    bool isWritten;

    std::ofstream file;

    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in,
           const std::string& fileFormat_in,
           const MathVector<int, 3>& rankMPI_in,
           const MathVector<int, 3>& sizeMPI_in)
      : writeFolder("../output/")
      , writerFolder(writerFolder_in)
      , fileExtension(fileExtension_in)
      , filePrefix(filePrefix_in)
      , fileFormat(fileFormat_in)
      , rankMPI(rankMPI_in)
      , sizeMPI(sizeMPI_in)
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
  private:
    using Base = Writer<T, InputOutput::Generic, InputOutputType::Generic,
                        InputOutputDataFormat::Generic>;

  public:
    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in,
           const MathVector<int, 3>& rankMPI_in,
           const MathVector<int, 3>& sizeMPI_in)
      : Base(writerFolder_in, filePrefix_in, fileExtension_in, "ascii",
             rankMPI_in, sizeMPI_in)
    {}

    using Base::getIsWritten;
    using Base::getIsBackedUp;

  protected:
    using Base::file;
    using Base::fileFormat;
    using Base::rankMPI;
    using Base::getFileName;

    inline void open(const std::string& fileName) {
      file.open(fileName, std::ofstream::out | std::ofstream::trunc);
      file.precision(16);

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
  private:
    using Base = Writer<T, InputOutput::Generic, InputOutputType::Generic,
                        InputOutputDataFormat::Generic>;

  public:
    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in,
           const MathVector<int, 3> rankMPI_in,
           const MathVector<int, 3> sizeMPI_in)
      : Writer<T, InputOutput::Generic, InputOutputType::Generic,
               InputOutputDataFormat::Generic>(writerFolder_in,
                                               filePrefix_in,
                                               fileExtension_in,
                                               "binary",
                                               rankMPI_in, sizeMPI_in)
    {}

    using Base::getIsWritten;
    using Base::getIsBackedUp;

  protected:
    using Base::file;
    using Base::fileFormat;
    using Base::rankMPI;
    using Base::sizeMPI;
    using Base::getFileName;

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
    : public Writer<T, InputOutput::Generic, InputOutputType::Generic,
                    inputOutputDataFormat> {
  private:
  using Base = Writer<T, InputOutput::Generic, InputOutputType::Generic,
                      inputOutputDataFormat>;
  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3> rankMPI_in,
           const MathVector<int, 3> sizeMPI_in)
      : Writer<T, InputOutput::Generic, InputOutputType::Generic,
               inputOutputDataFormat>("outputDAT/", filePrefix_in,
                                      ".dat", rankMPI_in, sizeMPI_in)
    {}

    using Base::getIsWritten;
    using Base::getIsBackedUp;

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
      void writeField(const unsigned int iteration, const T * data) {
      INSTRUMENT_ON("Writer<T, InputOutput::VTR, writerFileFromat>::writeField<NumberComponents>",3)

        write(iteration);

        for(auto i = 0; i < NumberComponents; ++i) {
          write(data[i]);
          file << " ";
        }
        file << std::endl;
    }

    template<unsigned int NumberComponents, DomainType initDomainType>
    void writeField(const Field<T, NumberComponents, initDomainType,
                    Architecture::CPU, false>& field) {
    }


  private:
    using Base::file;
    using Base::fileFormat;
    using Base::rankMPI;

    using Base::getFileName;
    using Base::open;
    using Base::write;


    void writeHeader() {
      INSTRUMENT_ON("Writer<T, InputOutput::VTR, InputOutputDataFormat::Generic>::writeHeader",2)

        file << "<?xml version=\"1.0\"?>\n";
      file << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" "
           << "byte_order=\"LittleEndian\">\n";
    }

  };


  template <class T, InputOutputDataFormat inputOutputDataFormat>
  class Writer<T, InputOutput::HDF5, InputOutputType::Parallel, inputOutputDataFormat>
    : public Writer<T, InputOutput::Generic, InputOutputType::Generic,
                    inputOutputDataFormat> {
  private:
    using Base = Writer<T, InputOutput::Generic, InputOutputType::Generic,
                        inputOutputDataFormat>;

    using Base::fileFormat;
    using Base::rankMPI;
    using Base::sizeMPI;
    using Base::getFileName;

    hid_t fileHDF5;
    hid_t datasetHDF5;
    herr_t statusHDF5;
    hid_t groupHDF5;
    hid_t dataSetHDF5;
    hid_t dataSpaceHDF5;
    hid_t fileSpaceHDF5;
    hid_t propertyListHDF5;
    Writer<T, InputOutput::XDMF, InputOutputType::Serial,
           inputOutputDataFormat> writerXDMF;

  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3>& rankMPI_in,
           const MathVector<int, 3>& sizeMPI_in)
      : Base("outputHDF5/", filePrefix_in,
             ".h5", rankMPI_in, sizeMPI_in)
      , writerXDMF(filePrefix_in, rankMPI_in, sizeMPI_in)
    {}

    using Base::getIsWritten;
    using Base::getIsBackedUp;

    inline void openFile(const unsigned int iteration) {
      std::string fileName = getFileName(iteration);
      open(fileName);

      if(rankMPI[d::X] == 0) {
        writerXDMF.openFile(iteration);
      }
    }

    inline void closeFile() {
      statusHDF5 = H5Fclose(fileHDF5);

      if(rankMPI[d::X] == 0) {
        writerXDMF.closeFile();
      }
    }

    template<unsigned int NumberComponents, Architecture architecture>
      void writeField(Field<T, NumberComponents, DomainType::LocalSpace,
                    architecture, true>& field) {
      INSTRUMENT_ON("Writer<T, InputOutput::HDF5, writerFileFromat>::writeField<NumberComponents>",3)

      propertyListHDF5 = H5Pcreate(H5P_DATASET_XFER);
      for(auto iC = 0; iC < NumberComponents; ++iC) {
        fileSpaceHDF5 = H5Screate_simple(L::dimD,
                                         Project<hsize_t,
                                         unsigned int, L::dimD>::Do(gSD::pLength()).data(),
                                         NULL);

        dataSetHDF5 = H5Dcreate2(fileHDF5, (field.fieldName+std::to_string(iC)).c_str(),
                                 H5T_NATIVE_DOUBLE, fileSpaceHDF5, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);

        statusHDF5 = H5Sclose(fileSpaceHDF5);

        dataSpaceHDF5 = H5Screate_simple(L::dimD,
                                         Project<hsize_t,
                                         unsigned int, L::dimD>::Do(lSD::pLength()).data(),
                                         NULL);

        fileSpaceHDF5 = H5Dget_space(dataSetHDF5);

        H5Sselect_hyperslab(fileSpaceHDF5, H5S_SELECT_SET,
                            Project<hsize_t,
                            unsigned int, L::dimD>::Do(gSD::pOffset(rankMPI)).data(), NULL,
                            Project<hsize_t,
                            unsigned int, L::dimD>::Do(lSD::pLength()).data(), NULL);

        H5Pset_dxpl_mpio(propertyListHDF5, H5FD_MPIO_COLLECTIVE);

        statusHDF5 = H5Dwrite(dataSetHDF5, H5T_NATIVE_DOUBLE,
                              dataSpaceHDF5, fileSpaceHDF5, propertyListHDF5,
                              field.getMultiData()[iC]);

        statusHDF5 = H5Dclose(dataSetHDF5);
        statusHDF5 = H5Sclose(dataSpaceHDF5);
        statusHDF5 = H5Sclose(fileSpaceHDF5);
      }

      statusHDF5 = H5Pclose(propertyListHDF5);

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

      if(!fileHDF5) {
        std::cout << "Could not open file " << fileName << std::endl;
      }
    }

  };


  template <class T, InputOutputDataFormat inputOutputDataFormat>
  class Writer<T, InputOutput::XDMF, InputOutputType::Serial, inputOutputDataFormat>
    : public Writer<T, InputOutput::Generic,
                    InputOutputType::Generic, inputOutputDataFormat> {
  private:
    using Base = Writer<T, InputOutput::Generic,
                        InputOutputType::Generic, inputOutputDataFormat>;

  public:
    Writer(const std::string& filePrefix_in,
           const MathVector<int, 3>& rankMPI_in,
           const MathVector<int, 3>& sizeMPI_in)
      : Base("outputHDF5/", filePrefix_in,
             ".xmf", rankMPI_in, sizeMPI_in)
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

        for(auto iC = 0; iC < NumberComponents; ++iC) {
          file << "<Attribute Name=\"" << field.fieldName+std::to_string(iC) << "\" "
               << "AttributeType=\"Scalar\" Center=\"Node\">\n";
          file << "<DataItem Dimensions=\""
               << gSD::pLength()[d::X];

          for(auto iD = 1; iD < L::dimD; ++iD) {
            file << " " << gSD::sLength()[iD];
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
    using Base::file;
    using Base::filePrefix;
    using Base::fileFormat;
    using Base::rankMPI;

    using Base::getFileName;
    using Base::open;
    using Base::write;

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
           << gSD::pLength()[d::X];

      for(auto iD = 1; iD < L::dimD; ++iD) {
        file << " " << gSD::sLength()[iD];
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

      for(auto iD = 1; iD < L::dimD; ++iD) {
        file << " 0";
      }
      file << "</DataItem>\n";

      file << "<DataItem Dimensions=\"" << L::dimD
           << "\" NumberType=\"Integer\" Format=\"XML\">"
           << "1";

      for(auto iD = 1; iD < L::dimD; ++iD) {
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


  typedef Writer<dataT, inputOutput, inputOutputType, inputOutputDataFormat> Writer_;

}

#endif // WRITER_H
