#pragma once

#include <mpi.h>
#include <cmath>
#include <memory>

#ifdef USE_FFTW
#include <fftw3-mpi.h>
#endif

#include "Commons.h"
#include "MathVector.h"

namespace lbm {

template <class T>
LBM_HOST LBM_DEVICE constexpr T& Min(T& a, T& b) {
  return a > b ? b : a;
}

template <class T>
LBM_HOST LBM_DEVICE constexpr T& arrayMin_impl(T* begin, T* end) {
  return begin + 1 == end ? *begin : Min(*begin, arrayMin_impl(begin + 1, end));
}

template <class T, std::size_t N>
LBM_HOST LBM_DEVICE constexpr T& arrayMin(T (&arr)[N]) {
  return arrayMin_impl(arr, arr + N);
}

template <class T>
LBM_HOST LBM_DEVICE constexpr T& Max(T& a, T& b) {
  return a > b ? a : b;
}

template <class T>
LBM_HOST LBM_DEVICE constexpr T& arrayMax_impl(T* begin, T* end) {
  return begin + 1 == end ? *begin : Max(*begin, arrayMin_impl(begin + 1, end));
}

template <class T, std::size_t N>
LBM_HOST LBM_DEVICE constexpr T& arrayMax(T (&arr)[N]) {
  return arrayMax_impl(arr, arr + N);
}

template <int Begin, int End, int Step = 1>
struct UnrolledFor {
  template <typename F>
  LBM_DEVICE LBM_HOST static void Do(F f) {
    f(Begin);
    UnrolledFor<Begin + Step, End, Step>::Do(f);
  }
};

template <int End>
struct UnrolledFor<End, End> {
  template <typename F>
  LBM_DEVICE LBM_HOST static void Do(F f) {}
};

template <class T>
LBM_DEVICE LBM_HOST inline T PowerBase(T arg, int power) {
  if (power == 1) {
    return arg;
  } else if (power == 0) {
    return (T)1;
  } else if (power == -1) {
    return (T)1.0 / arg;
  }
  return (T)pow(arg, power);
}

template <class T, int power>
class Power {
 public:
  LBM_DEVICE LBM_HOST static inline T Do(const T arg) {
    return (T)pow(arg, power);
  }
};

template <class T>
class Power<T, 0> {
 public:
  LBM_DEVICE LBM_HOST static inline T Do(const T arg) { return (T)1; }
};

template <class T>
class Power<T, 1> {
 public:
  LBM_DEVICE LBM_HOST static inline T Do(const T arg) { return (T)arg; }
};

template <class T>
class Power<T, 2> {
 public:
  LBM_DEVICE LBM_HOST static inline T Do(const T arg) { return (T)arg * arg; }
};

template <class T>
class Power<T, 3> {
 public:
  LBM_DEVICE LBM_HOST static inline T Do(const T arg) {
    return (T)arg * arg * arg;
  }
};

template <class T>
class Power<T, 4> {
 public:
  LBM_DEVICE LBM_HOST static inline T Do(const T arg) {
    return (T)arg * arg * arg * arg;
  }
};

template <class T>
class Power<T, -1> {
 public:
  LBM_DEVICE LBM_HOST static inline T Do(const T arg) { return (T)1.0 / arg; }
};

}  // namespace lbm
