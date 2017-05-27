#ifndef WRITER_H
#define WRITER_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include "Parameters.h"
#include "Options.h"
#include "StaticArray.h"
#include "MathVector.h"
#include "Field.h"
#include "Distribution.h"

namespace lbm {

  template <class T, WriterType writerType, DomainType domainType>
  class Writer {};

  template <class T>
  class Writer<T, WriterType::Generic, domainType> {
  protected:
    static constexpr Domain<domainType> domain;
    static constexpr std::string writeFolder = "../../output/";
    static constexpr std::string fileExtension = ".error";
    const std::string filePrefix;

    Writer(const std::string& filePrefix_in)
      : filePrefix(filePrefix_in)
    {}

    std::string getFileName(const unsigned int iteration){
      std::ostringstream number;

      number << fileNumber;
      std::string fileName =  writeFolder + writerFolder
        + filePrefix + number.str() + fileExtension;

      return fileName;
    }

  private:
    inline void writeField(Field<T, FieldType::Generic>& field) {
      // ABORT
      return;
    }

  public:
    inline void write(Fields<T>& field, const unsigned int iteration) {
      std::ofstream file;
      std::string fileName = getFileName(fiteration);
      file.open(fileName);

      if(file){
        writeFields(field);
      }
      else{
        BOOST_LOG_TRIVIAL(error) << "Writer VTR file: could not open " << fileName;
      }
    }
  };

  template <class T, DomainType domainType>
  class Writer<T, WriterType::VTR, domainType>
    : public Writer<T, WriterType::Generic, domainType> {
  private:
    using Writer<T, WriterType::Generic>::domain;
    using Writer<T, WriterType::Generic>::writeFolder;
    static constexpr std::string writerFolder = "outputVTR/";

    using Writer<T, WriterType::Generic>::getFileName;

    void writeFields(std::ofstream file, Field<T>& field) {
      writeHeader(file);


      writeFooter(file);
    }

  protected:
    static constexpr std::string fileExtension = ".vtr";

    void writeHeader(std::ofstream& file) {
      file << "<?xml version=\"1.0\"?>\n";
      file << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
      file << "<RectilinearGrid WholeExtent=\""
           << domain::start()[d::X] << " " << domain::end()[d::X]-1 << " "
           << domain::start()[d::Y] << " " << domain::end()[d::Y]-1 << " "
           << domain::start()[d::Z] << " " << domain::end()[d::Z]-1 << " \">\n";
      file << "\t<Piece Extent=\"";
           << domain::start()[d::X] << " " << domain::end()[d::X]-1 << " "
           << domain::start()[d::Y] << " " << domain::end()[d::Y]-1 << " "
           << domain::start()[d::Z] << " " << domain::end()[d::Z]-1 << " \">\n";
      file << "\t\t<PointData>\n";
    }

    void writeField(std::ofstream& file,
                    const Field<T>& field) {
      file << "\t\t\t<DataArray type=\"Float32\" "
           << "NumberOfComponents=\"" << field.numberComponents << "\" "
           << "Name=\"" << field.fieldName << "\" "
           << "format=\"ascii\">\n";
      for(int iZ = 0; iZ < L::lZ_g; iZ++) {
        for(int iY = 0; iY < L::lY_g; iY++) {
          for(int iX = 0; iX < L::lX_g; iX++) {
            MathVector<int, 3> iP = {iX, iY, iZ};
            file << field[iP]
          }
        }
      }
      file << "\t\t\t</DataArray>\n";
    }

    void writeFooter(std::ofstream& file) {
      file << "\t\t</PointData>\n" << "\t\t<CellData />\n" << "\t\t<Coordinates>\n";
      file << "\t\t\t<DataArray type=\"Float32\" Name=\"X\" format=\"ascii\">\n";
      for(int iX = domain::start()[d::X]; iX < domain::end()[d::X]; iX++){
        file << "\t\t\t\t" << iX+1 << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t\t<DataArray type=\"Float32\" Name=\"Y\" format=\"ascii\">\n";
      for(int iY = domain::start()[d::Y]; iY < domain::end()[d::Y]; iY++){
        file << "\t\t\t\t" << iY+1 << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t\t<DataArray type=\"Float32\" Name=\"Z\" format=\"ascii\">\n";
      for(int iZ = domain::start()[d::Z]; iZ < domain::end()[d::Z]; iZ++){
        file << "\t\t\t\t" << iZ+1 << "\n";
      }
      file << "\t\t\t</DataArray>\n";
      file << "\t\t</Coordinates>\n" << "\t</Piece>\n" << "</RectilinearGrid>\n";
      file << "</VTKFile>\n";
    }

  public:
    Writer(const std::string& filesPrefix_in)
      : Writer<T, WriterType::Generic>(filesPrefix_in)
    {}
  };

  template <class T, DomainType domainType>
  class Writers {
  private:
    StaticArray<Writer<T, WriterType::Generic, domainType>, numberWriters> writersArray;

  public:
    Writers(const StaticArray<Writer<T, WriterType::Generic, domainType>,
            numberWriters>& writersArray_in)
      : writersArray(writersArray_in)
    {}

    inline void write(Field<T>& field, const int fileNumber){
      for(auto writer : writersArray){
        writer.write(field, fileNumber);
      }
    }

  };

}

#endif // WRITER_H
