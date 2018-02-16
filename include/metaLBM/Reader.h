#ifndef READER_H
#define READER_H

#include <string>

#include <rapidxml.hpp>
#include <rapidxml_utils.hpp>

#include "Commons.h"
#include "Options.h"
#include "MathVector.h"
#include "Domain.h"
#include "Field.h"

namespace lbm {

  template<class T, unsigned int NumberComponents, InputOutput inputOutput,
           InputOutputType inputOutputType, InputOutputDataFormat inputOutputDataFormat>
  class Reader {};

  template<class T, unsigned int NumberComponents>
  class Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
               InputOutputDataFormat::Generic> {
  protected:
    const std::string readFolder;
    const std::string readerFolder;
    const std::string fileExtension;
    const std::string filePrefix;

    Reader(const std::string& readerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in)
      : readFolder("../../output/")
      , readerFolder(readerFolder_in)
      , fileExtension(fileExtension_in)
      , filePrefix(filePrefix_in)
    {}

    std::string getFileName(const unsigned int iteration){
      std::ostringstream number;

      number << iteration;
      std::string fileName =  readFolder + readerFolder
        + filePrefix + '-' + number.str() + fileExtension;

      return fileName;
    }

  public:
    inline DynamicArray<T, Architecture::CPU> readArray(const std::string& fieldName,
                                                        const unsigned int iteration) {
      // ABORT
      return DynamicArray<T, Architecture::CPU>();
    }
  };


  template <class T, unsigned int NumberComponents>
  class Reader<T, NumberComponents, InputOutput::VTR, InputOutputType::Serial,
               InputOutputDataFormat::ascii>
    : public Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
                    InputOutputDataFormat::Generic> {
  private:
    using Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getFileName;
    typedef Domain<DomainType::Global, partitionningT,
                   MemoryLayout::Generic, NumberComponents> gNCD;

  public:
    Reader(const std::string& filePrefix_in)
      : Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
               InputOutputDataFormat::Generic>("outputVTR/", filePrefix_in, ".vtr")
    {}

    inline DynamicArray<T, Architecture::CPU> readArray(const std::string& fieldName,
                                                        const unsigned int iteration) {
      INSTRUMENT_ON("Reader<T, NumberComponents, InputOutput::VTR>::readArray",3)

      DynamicArray<T, Architecture::CPU> arrayR(gD::volume());

      std::string fileName = getFileName(iteration);
      rapidxml::file<> xmlFile(fileName.c_str());
      rapidxml::xml_document<> document;
      document.parse<0>(xmlFile.data());

      rapidxml::xml_node<>* node = document.first_node("VTKFile")->first_node("RectilinearGrid")->first_node("Piece")->first_node("PointData")->first_node("DataArray");
      std::string content = node->value();
      std::istringstream file(content);

      std::string line;
      std::getline(file, line);

      for(unsigned int index = 0; index < gD::volume(); ++index) {
        MathVector<unsigned int, 3> iP;

        unsigned int indexTemporary = index;
        iP[d::X] = indexTemporary%gD::length()[d::X];
        indexTemporary = (indexTemporary-iP[d::X])/gD::length()[d::X];
        iP[d::Y] = indexTemporary%gD::length()[d::Y];
        indexTemporary = (indexTemporary-iP[d::Y])/gD::length()[d::Y];
        iP[d::Z] = indexTemporary;

        std::getline(file, line);
        std::string::size_type offset, offsetTemporary;
        T value = std::stod(line, &offset);
        offsetTemporary = offset;
        arrayR[gNCD::getIndex(iP, 0)] = value;

        for(unsigned int iQ = 1; iQ < NumberComponents; ++iQ) {
          value = std::stod(line.substr(offset), &offset);
          offset += offsetTemporary;
          offsetTemporary = offset;
          arrayR[gNCD::getIndex(iP, iQ)] = value;
        }
      }

      return arrayR;
    }

  };


  template <class T, unsigned int NumberComponents,
            InputOutputDataFormat inputOutputDataFormat>
    class Reader<T, NumberComponents, InputOutput::HDF5, InputOutputType::Serial,
                 inputOutputDataFormat>
    : public Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
                    InputOutputDataFormat::Generic> {
  private:
    using Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getFileName;
    typedef Domain<DomainType::Global, partitionningT,
                   MemoryLayout::Generic, NumberComponents> gNCD;

  public:
    Reader(const std::string& filePrefix_in)
      : Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
               InputOutputDataFormat::Generic>("outputHDF5/", filePrefix_in, ".h5")
    {}

    inline DynamicArray<T, Architecture::CPU> readArray(const std::string& fieldName,
                                                        const unsigned int iteration) {
      INSTRUMENT_ON("Reader<T, NumberComponents, InputOutput::VTR>::readArray",3)

      DynamicArray<T, Architecture::CPU> arrayR(gD::volume());

      return arrayR;
    }

  };


  template <class T, unsigned int NumberComponents,
            InputOutputDataFormat inputOutputDataFormat>
    class Reader<T, NumberComponents, InputOutput::HDF5, InputOutputType::Parallel,
                 inputOutputDataFormat>
    : public Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
                    InputOutputDataFormat::Generic> {
  private:
    using Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getFileName;
    typedef Domain<DomainType::Global, partitionningT,
                   MemoryLayout::Generic, NumberComponents> gNCD;

  public:
    Reader(const std::string& filePrefix_in)
      : Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
               InputOutputDataFormat::Generic>("outputHDF5/", filePrefix_in, ".h5")
    {}

    inline DynamicArray<T, Architecture::CPU> readArray(const std::string& fieldName,
                                                        const unsigned int iteration) {
      INSTRUMENT_ON("Reader<T, NumberComponents, InputOutput::VTR>::readArray",3)

      DynamicArray<T, Architecture::CPU> arrayR(gD::volume());

      return arrayR;
    }

  };


}


#endif // READER_H
