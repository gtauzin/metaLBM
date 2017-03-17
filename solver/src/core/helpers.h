#ifndef HELPERS_H
#define HELPERS_H

#include <memory>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include<iostream>

#include <rapidxml.hpp>
#include <rapidxml_utils.hpp>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include <boost/align/aligned_allocator.hpp>
#include <vector>
template<class T, std::size_t Alignment = 1>
#ifdef ALIGN_MEMORY
  using vector = std::vector<T,
  boost::alignment::aligned_allocator<T, Alignment> >;
#else
using vector = std::vector<T>;
#endif

#include "commons.h"
#include "structure.h"
#include "lattice.h"

namespace boost {
  namespace log{

    template <class T>
      inline logging::formatting_ostream& operator<<(log::formatting_ostream& ostream,
                                                     std::vector<T>& vector) {
      ostream << "[ ";
      for(auto i : vector){
        ostream << " " << i << " ,";
      }
      ostream << "]";
      return ostream;
    }

    template <class T, unsigned int size>
      inline logging::formatting_ostream& operator<<(log::formatting_ostream& ostream,
                                                     std::array<T, size>& array) {
      ostream << "[ ";
      lbm::UnrolledFor<0, size>::Do([&] (int i) {
          ostream << " " << i << " ,";
        });
      ostream << "]";
      return ostream;
    }

  }
}

namespace lbm {

  template <class T, LatticeType L>
    vector<T, CACHE_LINE> readVTK(const std::string& inputFilename,
                                  const std::string& dataArrayName) {

    vector<T, CACHE_LINE> output(s_g() * P::dimQ);

    rapidxml::file<> xmlFile(inputFilename.c_str());
    rapidxml::xml_document<> doc;
    doc.parse<0>(xmlFile.data());

    rapidxml::xml_node<>* distribution_node = doc.first_node("VTKFile")->first_node("RectilinearGrid")->first_node("Piece")->first_node("PointData")->first_node("DataArray");
    std::string content = distribution_node->value();
    std::istringstream f(content);

    std::string line;
    std::getline(f, line);

    for(int tuple = 0; tuple < s_g(); ++tuple) {
      std::getline(f,line);
      std::string::size_type offset, tmp_offset;
      T value = std::stod(line, &offset);
      tmp_offset = offset;
      output[tuple * P::dimQ] = value;

      for(int iQ = 1; iQ < P::dimQ; ++iQ) {
        value = std::stod(line.substr(offset), &offset);
        offset += tmp_offset;
        tmp_offset = offset;
        output[tuple * P::dimQ + iQ] = value;
      }
    }


    return output;
  }


  template <class T>
    struct RootFinderFunctor {
    public:
      RootFinderFunctor(){};

#pragma omp declare simd
      virtual T evaluateFunction(T const& x) = 0;

      virtual T evaluateDerivative(T const& x) = 0;

    };

#pragma omp declare simd
  template <class T>
    inline bool NewtonRaphsonSolver(std::shared_ptr<RootFinderFunctor<T>> functor,
                                    const T tolerance, const int iterationMax,
                                    T& xR, const T xMin, const T xMax) {

    BOOST_LOG_TRIVIAL(debug) << "xR: " << xR
                             << ", xMin: " << xMin
                             << ", xMax: " << xMax
                             << ", tolerance: " << tolerance
                             << ", iterationMax: " << iterationMax;


    T error = 1 + tolerance;
    T xStep = 0.0;

    for(int iteration = 1; iteration <= iterationMax; ++iteration) {
      xR = xR - xStep;

      T functionEvaluation = functor->evaluateFunction(xR);
      T derivativeEvaluation = functor->evaluateDerivative(xR);
      xStep = functionEvaluation/derivativeEvaluation;

      error = fabs(xStep);

      BOOST_LOG_TRIVIAL(debug) << "iteration: " << iteration
                               << ", error: " << error
                               << ", xR: " << xR;

      if(error <= tolerance) {
        if(xR > xMin && xR < xMax) {
          return true;
        }

        else {
          return false;
        }
      }
    }

    return false;

  }

#pragma omp declare simd
  template <class T>
    inline bool Bisection_NewtonRaphsonSolver(std::shared_ptr<RootFinderFunctor<T>> functor,
                                              const T tolerance, const int iterationMax,
                                              T& xR, const T xMin, const T xMax) {

    BOOST_LOG_TRIVIAL(debug) << "xR: " << xR
                             << ", xMin: " << xMin
                             << ", xMax: " << xMax
                             << ", tolerance: " << tolerance
                             << ", iterationMax: " << iterationMax;

    T xLow, xHigh;
    T function_xLow = functor->evaluateFunction(xMin);
    T function_xHigh = functor->evaluateFunction(xMax);

    if ((function_xLow > 0.0 && function_xHigh > 0.0)
        || (function_xLow < 0.0 && function_xHigh < 0.0)) {
      BOOST_LOG_TRIVIAL(debug) << "Root must be in [xMin, xMax]";
      return false;
    }

    if (function_xLow == 0.0) {
      xR = xMin;
      return true;
    }

    if (function_xHigh == 0.0) {
      xR = xMax;
      return true;
    }

    if (function_xLow < 0.0) {
      xLow = xMin;
      xHigh = xMax;
    }

    else {
      xLow = xMax;
      xHigh = xMin;
    }

    xR = 0.5 * (xMin + xMax);
    //T xR = 2.0;
    T xStepPrevious = fabs(xMax-xMin);
    T xStep = xStepPrevious;
    T functionEvaluation = functor->evaluateFunction(xR);
    T derivativeEvaluation = functor->evaluateDerivative(xR);

    for(int iteration = 1; iteration <= iterationMax; ++iteration) {

      if ((((xR-xHigh)*derivativeEvaluation-functionEvaluation)
           *((xR-xLow)*derivativeEvaluation-functionEvaluation) > 0.0)
          || (fabs(2.0*functionEvaluation) > fabs(xStepPrevious*derivativeEvaluation))) {

        xStepPrevious = xStep;
        xStep = 0.5 * (xHigh-xLow);
        xR = xLow + xStep;

        if (xLow == xR) {
          return true;
        }
      }

      else {
        xStepPrevious = xStep;
        xStep = functionEvaluation/derivativeEvaluation;
        T xTemp = xR;
        xR -= xStep;

        if (xTemp == xR) {
          return true;
        }
      }

      BOOST_LOG_TRIVIAL(debug) << "iteration: " << iteration
                               << ", error: " << fabs(xStep)
                               << ", xR: " << xR;

      if(fabs(xStep) <= tolerance) {
        return true;
      }

      functionEvaluation = functor->evaluateFunction(xR);
      derivativeEvaluation = functor->evaluateDerivative(xR);

      if (functionEvaluation < 0.0) {
        xLow = xR;
      }
      else {
        xHigh = xR;
      }
    }

    BOOST_LOG_TRIVIAL(debug) << "Maximum number of iterations exceeded";
    return false;

  }

}


#endif // HELPERS_H
