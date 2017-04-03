#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <array>
#include <numeric>
#include <memory>

#include <boost/align/aligned_allocator.hpp>
#include <vector>
template<class T, std::size_t Alignment = 1>
#ifdef ALIGN_MEMORY
  using vector = std::vector<T,
  boost::alignment::aligned_allocator<T, Alignment> >;
#else
  using vector = std::vector<T>;
#endif

#include "input.h"
#include "commons.h"

namespace lbm {

  template <class T>
  class Distribution {
  public:
    vector<T, CACHE_LINE> fPop;

    Distribution(const vector<T, CACHE_LINE>& fPop_in)
      : fPop(fPop_in)
    {}

    T& operator[] (int i) {
      return fPop[i];
    }

    const T& operator[] (int i) const {
      return fPop[i];
    }
  };
}

#endif // DISTRIBUTION_H
