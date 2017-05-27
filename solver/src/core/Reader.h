#ifndef READ_H
#define READ_H

#include <string>

#include <rapidxml.hpp>
#include <rapidxml_utils.hpp>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/sources/severity_feature.hpp>
namespace logging = boost::log;

#include <Field.h>

namespace lbm {

  template <class T, unsigned int NumberComponents>
    Field<T, NumberComponents, true> readVTK(const std::string& inputFilename,
                                             const std::string& dataArrayName) {

    Field<T, NumberCOmponents, true> fieldR(s_g() * L::dimQ, dataArrayName);

    rapidxml::file<> xmlFile(inputFilename.c_str());
    rapidxml::xml_document<> doc;
    doc.parse<0>(xmlFile.data());

    rapidxml::xml_node<>* distribution_node = doc.first_node("VTKFile")->first_node("RectilinearGrid")
      ->first_node("Piece")->first_node("PointData")->first_node("DataArray");
    std::string content = distribution_node->value();
    std::istringstream f(content);

    std::string line;
    std::getline(f, line);

    for(int tuple = 0; tuple < s_g(); ++tuple) {
      std::getline(f,line);
      std::string::size_type offset, tmp_offset;
      T value = std::stod(line, &offset);
      tmp_offset = offset;
      output[tuple * L::dimQ] = value;

      for(int iQ = 1; iQ < L::dimQ; ++iQ) {
        value = std::stod(line.substr(offset), &offset);
        offset += tmp_offset;
        tmp_offset = offset;
        output[tuple * L::dimQ + iQ] = value;
      }
    }

    return FieldR;
  }
}


#endif // READ_H
