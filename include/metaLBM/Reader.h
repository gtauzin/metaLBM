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


  template <class T, unsigned int NumberComponents,
            InputOutputDataFormat inputOutputDataFormat>
    class Reader<T, NumberComponents, InputOutput::HDF5, InputOutputType::Parallel,
                 inputOutputDataFormat>
    : public Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
                    InputOutputDataFormat::Generic> {
  private:
    using Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
                 InputOutputDataFormat::Generic>::getFileName;

  public:
    Reader(const std::string& filePrefix_in)
      : Reader<T, NumberComponents, InputOutput::Generic, InputOutputType::Generic,
               InputOutputDataFormat::Generic>("outputHDF5/", filePrefix_in, ".h5")
    {}

    inline MultiDynamicArray<T, Architecture::CPU,
                             NumberComponents> readArray(const std::string& fieldName,
                                                         const unsigned int iteration) {
      { INSTRUMENT_ON("Reader<T, NumberComponents, InputOutput::VTR>::readArray",3) }

      MultiDynamicArray<T, Architecture::CPU, NumberComponents> arrayR(gSD::sVolume());



      return arrayR;
    }

  };


}


#endif // READER_H
