#ifndef LOGGING_H
#define LOGGING_H

#include <iostream>
#include <ostream>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
namespace logging = boost::log;

#include "StaticArray.h"
#include "DynamicArray.h"

namespace lbm {

  void initLogging(const int mpi_rank) {
#ifdef ENABLE_LOG
    std::cout << "Logging enabled!\n";
    BOOST_LOG_TRIVIAL(debug) << "Logging for debug starts.";

    std::ostringstream rank;
    rank << mpi_rank;
    logging::add_file_log(keywords::file_name = "../log/log_rank-" + rank.str() + ".log",
                          keywords::rotation_size = 10 * 1024 * 1024,
                          keywords::time_based_rotation
                          = sinks::file::rotation_at_time_point(0, 0, 0),
                          keywords::format = "[%TimeStamp%]: %Message%"
                          );

    logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::debug);

#else
    logging::core::get()->set_logging_enabled(false);

#endif
  }
}

namespace boost {
  namespace log{

    template <class T>
      inline logging::formatting_ostream& operator<<(log::formatting_ostream& ostream,
                                                     DynamicArray<T>& dArray) {
      ostream << dArray;
      return ostream;
    }

    template <class T, unsigned int Size>
      inline logging::formatting_ostream& operator<<(log::formatting_ostream& ostream,
                                                     StaticArrat<T, Size>& array) {
      ostream << sArray;
      return ostream;
    }

  }
}

#endif // LOGGING_H
