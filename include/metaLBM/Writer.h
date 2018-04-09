#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <hdf5.h>
#include <sys/stat.h>

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "Domain.h"
#include "Field.h"
#include "Distribution.h"

namespace lbm {

  template <class T, InputOutput inputOutput, InputOutputFormat inputOutputFormat>
  class Writer {};

  template <class T>
  class Writer<T, InputOutput::Generic, InputOutputFormat::Generic> {
  protected:
    const std::string writeFolder;
    const std::string writerFolder;
    const std::string fileExtension;
    const std::string filePrefix;
    const std::string fileFormat;

    std::ofstream file;

    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in,
           const std::string& fileFormat_in)
      : writeFolder("../output/")
      , writerFolder(writerFolder_in)
      , fileExtension(fileExtension_in)
      , filePrefix(filePrefix_in)
      , fileFormat(fileFormat_in)
    {}

    inline std::string getFileName(const unsigned int iteration) {
      std::ostringstream number;

      number << iteration;
      return writeFolder + writerFolder + filePrefix + "-" + number.str() + fileExtension;
    }

    inline std::string getFileName() {
      return writeFolder + writerFolder + filePrefix + fileExtension;
    }


  public:
    inline bool getIsWritten(const unsigned int iteration) {
      return (iteration % writeStep) == 0;
    }

  };


  template <class T>
  class Writer<T, InputOutput::Generic, InputOutputFormat::ascii>
    : public Writer<T, InputOutput::Generic, InputOutputFormat::Generic> {
  private:
    using Base = Writer<T, InputOutput::Generic, InputOutputFormat::Generic>;

  public:
    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in)
      : Base(writerFolder_in, filePrefix_in, fileExtension_in, "ascii")
    {}

    using Base::getIsWritten;

  protected:
    inline void openAndAppend(const std::string& fileName) {
      Base::file.open(fileName, std::ofstream::out | std::ofstream::app);
      Base::file.precision(16);

      if(!Base::file) {
        std::cout << "Could not open file " << fileName << std::endl;
      }

    }

    inline void openAndTruncate(const std::string& fileName) {
      Base::file.open(fileName, std::ofstream::out | std::ofstream::trunc);
      Base::file.precision(16);

      if(!Base::file) {
        std::cout << "Could not open file " << fileName << std::endl;
      }

    }

    template<class U>
    inline void write(const U data) {
      Base::file << data;
    }

  };


  template <class T>
  class Writer<T, InputOutput::Generic, InputOutputFormat::binary>
    : public Writer<T, InputOutput::Generic, InputOutputFormat::Generic> {
  private:
    using Base = Writer<T, InputOutput::Generic, InputOutputFormat::Generic>;

  public:
    Writer(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in)
      : Base(writerFolder_in, filePrefix_in, fileExtension_in,
             "binary")
    {}

    using Base::getIsWritten;
    using Base::getIsBackedUp;

  protected:
    inline void openAndAppend(const std::string& fileName) {
      Base::file.open(fileName, std::ofstream::out | std::ofstream::app
                      | std::ofstream::binary);
    }

    inline void openAndTruncate(const std::string& fileName) {
      Base::file.open(fileName, std::ofstream::out | std::ofstream::trunc
                      | std::ofstream::binary);
    }

    template<class U>
    inline void write(U data) {
      Base::file.write(reinterpret_cast<char*>(&data), sizeof(data));
    }

  };

  template <class T, InputOutputFormat inputOutputFormat>
  class ScalarAnalysisWriter
    : public Writer<T, InputOutput::Generic, inputOutputFormat> {
  private:
    using Base = Writer<T, InputOutput::Generic, inputOutputFormat>;

  public:
    ScalarAnalysisWriter(const std::string& filePrefix_in)
      : Base(filePrefix_in+"/", "observables", ".dat")
    {}

    inline bool getIsAnalyzed(const unsigned int iteration) {
      return (iteration % scalarAnalysisStep) == 0;
    }


    inline void openFile(const unsigned int iteration) {
      Base::openAndAppend(Base::getFileName());
    }

    inline void closeFile() {
      Base::file.close();
    }

    template<unsigned int NumberScalarAnalyses>
    void writeAnalysis(const unsigned int iteration, T * data) {
      { INSTRUMENT_ON("Writer<T, InputOutputFormat>::writeAnalysis<NumberScalarAnalysis>",3) }

      Base::write(iteration);
      Base::file << " ";

      for(auto iS = 0; iS < NumberScalarAnalyses; ++iS) {
        Base::write(data[iS]);
        Base::file << " ";
      }
      Base::file << std::endl;
    }

    void writeHeader(const std::string& header) {
      { INSTRUMENT_ON("AnalysisWriter::writeHeader",2) }
      Base::openAndTruncate(Base::getFileName());

      Base::file << header << std::endl;

      closeFile();
    }

  };

  template <class T, InputOutputFormat inputOutputFormat>
  class SpectralAnalysisWriter
    : public Writer<T, InputOutput::Generic, inputOutputFormat> {
  private:
    using Base = Writer<T, InputOutput::Generic, inputOutputFormat>;

  public:
    SpectralAnalysisWriter(const std::string& filePrefix_in)
      : Base(filePrefix_in+"/", "spectra", ".dat")
    {}

    inline bool getIsAnalyzed(const unsigned int iteration) {
      return (iteration % spectralAnalysisStep) == 0;
    }

    inline void openFile(const unsigned int iteration) {
      Base::openAndAppend(Base::getFileName());
    }

    inline void closeFile() {
      Base::file.close();
    }

    template<unsigned int NumberSpectralAnalyses, unsigned int MaxWaveNumber>
    void writeAnalysis(const unsigned int iteration, T * data[NumberSpectralAnalyses]) {
      { INSTRUMENT_ON("Writer<T, InputOutput::DAT, writerFileFromat>::writeAnalysis<NumberComponents>",3) }

      for(auto kNorm = 0; kNorm < MaxWaveNumber; ++kNorm) {
        Base::write(iteration);
        Base::file << " ";

        Base::write(kNorm);
          Base::file << " ";

          for(auto iS = 0; iS < NumberSpectralAnalyses; ++iS) {
            Base::write(data[iS][kNorm]);
            Base::file << " ";
          }
          Base::file << std::endl;
      }
    }

    void writeHeader(const std::string& header) {
      { INSTRUMENT_ON("AnalysisWriter::writeHeader",2) }
      Base::openAndTruncate(Base::getFileName());

      Base::file << header << std::endl;

      closeFile();
    }

  };


  template <class T, InputOutput inputOutput>
  class FieldWriter {};

  template <class T>
  class FieldWriter<T, InputOutput::HDF5>
    : public Writer<T, InputOutput::Generic, InputOutputFormat::Generic> {
  private:
    using Base = Writer<T, InputOutput::Generic, InputOutputFormat::Generic>;

  protected:
    hid_t fileHDF5;
    hid_t datasetHDF5;
    herr_t statusHDF5;
    hid_t groupHDF5;
    hid_t dataSetHDF5;
    hid_t dataSpaceHDF5;
    hid_t fileSpaceHDF5;
    hid_t propertyListHDF5;
    FieldWriter<T, InputOutput::XDMF> writerXDMF;
    const MathVector<int, 3> rankMPI;

  public:
    FieldWriter(const std::string& filePrefix_in,
                const MathVector<int, 3>& rankMPI_in,
                const std::string& name_in = "field")
      : Base(filePrefix_in+"/", name_in, ".h5", "binary")
      , rankMPI(rankMPI_in)
      , writerXDMF(filePrefix_in, name_in)
    {
      if(rankMPI[d::X] == 0) {
        int dirError = mkdir(Base::writeFolder.c_str(), S_IRWXU | S_IRWXG
                             | S_IROTH | S_IXOTH);
        dirError = mkdir((Base::writeFolder+Base::writerFolder).c_str(), S_IRWXU | S_IRWXG
                         | S_IROTH | S_IXOTH);
        if (dirError == -1)
          {
            //std::cout << "Error creating directory! It probably already exists..."
            //          << std::endl;
          }
      }

    }

    using Base::getIsWritten;

    inline void openFile(const unsigned int iteration) {
      std::string fileName = Base::getFileName(iteration);
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
    void writeField(Field<T, NumberComponents, architecture, true>& field) {
      { INSTRUMENT_ON("Writer<T, InputOutput::HDF5, writerFileFromat>::writeField<NumberComponents>",3) }

      std::string fieldName = field.fieldName;
      propertyListHDF5 = H5Pcreate(H5P_DATASET_XFER);

      for(auto iC = 0; iC < NumberComponents; ++iC) {
        if(NumberComponents > 1) {
          fieldName = field.fieldName + dName[iC];
        }
        fileSpaceHDF5 = H5Screate_simple(L::dimD,
                                         Project<hsize_t,
                                         unsigned int, L::dimD>::Do(gSD::pLength()).data(),
                                         NULL);

        dataSetHDF5 = H5Dcreate2(fileHDF5, (fieldName).c_str(),
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
                            unsigned int, L::dimD>::Do(gSD::pOffset(rankMPI)).data(),
                            NULL,
                            Project<hsize_t,
                            unsigned int, L::dimD>::Do(lSD::pLength()).data(), NULL);

        H5Pset_dxpl_mpio(propertyListHDF5, H5FD_MPIO_COLLECTIVE);

        statusHDF5 = H5Dwrite(dataSetHDF5, H5T_NATIVE_DOUBLE,
                              dataSpaceHDF5, fileSpaceHDF5, propertyListHDF5,
                              field.getLocalData(iC));

        statusHDF5 = H5Dclose(dataSetHDF5);
        statusHDF5 = H5Sclose(dataSpaceHDF5);
        statusHDF5 = H5Sclose(fileSpaceHDF5);

        if(rankMPI[d::X] == 0) {
          writerXDMF.write(fieldName, NumberComponents);
        }

      }

      statusHDF5 = H5Pclose(propertyListHDF5);
    }

    template<unsigned int NumberComponents, Architecture architecture>
    void writeField(Field<T, NumberComponents, architecture, false>& field) {
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


  template <class T, InputOutput inputOutput>
  class DistributionWriter {};

  template <class T>
  class DistributionWriter<T, InputOutput::HDF5>
    : public FieldWriter<T, InputOutput::HDF5> {
  private:
    using Base = FieldWriter<T, InputOutput::HDF5>;

  public:
    DistributionWriter(const std::string& filePrefix_in,
                       const MathVector<int, 3>& rankMPI_in)
      : Base(filePrefix_in, rankMPI_in, "distribution")
    {}

    inline bool getIsBackedUp(const unsigned int iteration) {
      return (iteration % backUpStep) == 0;
    }

    template<Architecture architecture>
    void writeDistribution(Distribution<T, architecture>& distribution) {
      { INSTRUMENT_ON("Writer<T, InputOutput::HDF5, writerFileFromat>::writeField<NumberComponents>",3) }

      Base::propertyListHDF5 = H5Pcreate(H5P_DATASET_XFER);
      for(auto iC = 0; iC < L::dimQ; ++iC) {
        Base::fileSpaceHDF5 = H5Screate_simple(L::dimD,
                                               Project<hsize_t,
                                               unsigned int, L::dimD>::Do(gSD::pLength()).data(),
                                               NULL);

        Base::dataSetHDF5 = H5Dcreate2(Base::fileHDF5, (distribution.fieldName+std::to_string(iC)).c_str(),
                                       H5T_NATIVE_DOUBLE, Base::fileSpaceHDF5, H5P_DEFAULT,
                                       H5P_DEFAULT, H5P_DEFAULT);

        Base::statusHDF5 = H5Sclose(Base::fileSpaceHDF5);

        Base::dataSpaceHDF5 = H5Screate_simple(L::dimD,
                                               Project<hsize_t,
                                               unsigned int, L::dimD>::Do(lSD::pLength()).data(),
                                               NULL);

        Base::fileSpaceHDF5 = H5Dget_space(Base::dataSetHDF5);

        H5Sselect_hyperslab(Base::fileSpaceHDF5, H5S_SELECT_SET,
                            Project<hsize_t,
                            unsigned int, L::dimD>::Do(gSD::pOffset(Base::rankMPI)).data(),
                            NULL,
                            Project<hsize_t,
                            unsigned int, L::dimD>::Do(lSD::pLength()).data(), NULL);

        H5Pset_dxpl_mpio(Base::propertyListHDF5, H5FD_MPIO_COLLECTIVE);

        Base::statusHDF5 = H5Dwrite(Base::dataSetHDF5, H5T_NATIVE_DOUBLE,
                                    Base::dataSpaceHDF5, Base::fileSpaceHDF5,
                                    Base::propertyListHDF5,
                                    distribution.getLocalData(iC));

        Base::statusHDF5 = H5Dclose(Base::dataSetHDF5);
        Base::statusHDF5 = H5Sclose(Base::dataSpaceHDF5);
        Base::statusHDF5 = H5Sclose(Base::fileSpaceHDF5);

        if(Base::rankMPI[d::X] == 0) {
          Base::writerXDMF.write(distribution.fieldName+std::to_string(iC), L::dimQ);
        }

      }

      Base::statusHDF5 = H5Pclose(Base::propertyListHDF5);

    }

    using Base::openFile;
    using Base::closeFile;

  };


  template <class T>
  class FieldWriter<T, InputOutput::XDMF>
    : public Writer<T, InputOutput::Generic, InputOutputFormat::ascii> {
  private:
    using Base = Writer<T, InputOutput::Generic, InputOutputFormat::ascii>;

  public:
    FieldWriter(const std::string& filePrefix_in,
                const std::string& name_in)
      : Base(filePrefix_in+"/", name_in,
             ".xmf")
      , fileName("/dev/null")
      , fileNameHDF5("/dev/null")
    {}

    inline void openFile(const unsigned int iteration) {
      fileName = Base::getFileName(iteration);
      fileNameHDF5 = getFileNameHDF5(iteration);
      Base::openAndTruncate(fileName);
      writeHeader();
    }

    inline void closeFile() {
      writeFooter();
      Base::file.close();
    }

    void write(const std::string& name, unsigned int numberComponents) {
      { INSTRUMENT_ON("Writer<T, InputOutput::XDMF, writerFileFromat>::writeField<NumberComponents>",3) }

      Base::file << "<Attribute Name=\"" << name << "\" "
                 << "AttributeType=\"Scalar\" Center=\"Node\">\n";
      Base::file << "<DataItem Dimensions=\""
                 << gSD::pLength()[d::X];

      for(auto iD = 1; iD < L::dimD; ++iD) {
        Base::file << " " << gSD::sLength()[iD];
      }
      Base::file << "\" ";
      Base::file << "NumberType=\"Double\" Precision=\"8\" Format=\"HDF\">\n";
      Base::file << fileNameHDF5 << ":/" << name << "\n";
      Base::file << "</DataItem>\n";
      Base::file << "</Attribute>\n";
    }

  private:
    std::string fileName;
    std::string fileNameHDF5;

    std::string getFileNameHDF5(const unsigned int iteration) {
      std::ostringstream number;

      number << iteration;
      std::string fileNameR =  Base::filePrefix + "-" + number.str() + ".h5";

      return fileNameR;
    }

    void writeHeader() {
      INSTRUMENT_ON("Writer<T, InputOutput::XDMF, InputOutputFormat::Generic>::writeHeader",2)

        Base::file << "<?xml version=\"1.0\"?>\n";
      Base::file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
      Base::file << "<Xdmf>\n";
      Base::file << "<Domain>\n";
      Base::file << "<Grid Name=\"grid\" GridType=\"Uniform\">\n";
      Base::file << "<Topology TopologyType=\"" << L::dimD << "DCoRectMesh\" Dimensions=\""
                 << gSD::pLength()[d::X];

      for(auto iD = 1; iD < L::dimD; ++iD) {
        Base::file << " " << gSD::sLength()[iD];
      }
      Base::file << "\" />\n"
                 << "<Geometry GeometryType=\"Origin_Dx";
      if(L::dimD >= 2) {
        Base::file << "Dy";
      }
      if(L::dimD == 3) {
        Base::file << "Dz";
      }
      Base::file << "\">\n"
                 << "<DataItem Dimensions=\"" << L::dimD
                 << "\" NumberType=\"Integer\" Format=\"XML\">"
                 << "0";

      for(auto iD = 1; iD < L::dimD; ++iD) {
        Base::file << " 0";
      }
      Base::file << "</DataItem>\n";

      Base::file << "<DataItem Dimensions=\"" << L::dimD
                 << "\" NumberType=\"Integer\" Format=\"XML\">"
                 << "1";

      for(auto iD = 1; iD < L::dimD; ++iD) {
        Base::file << " 1";
      }
      Base::file << "</DataItem>\n"
                 << "</Geometry>\n";
    }

    void writeFooter() {
      { INSTRUMENT_ON("Writer<T, InputOutput::XDMF, InputOutputFormat::Generic>::writeFooter",2) }
      Base::file << "</Grid>\n"
                 << "</Domain>\n"
                 << "</Xdmf>\n";
    }

  };


  typedef FieldWriter<dataT, InputOutput::HDF5> FieldWriter_;
  typedef DistributionWriter<dataT, InputOutput::HDF5> DistributionWriter_;
  typedef ScalarAnalysisWriter<dataT, InputOutputFormat::ascii> ScalarAnalysisWriter_;
  typedef SpectralAnalysisWriter<dataT, InputOutputFormat::ascii> SpectralAnalysisWriter_;
}
