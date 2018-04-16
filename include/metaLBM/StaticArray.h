#pragma once

#include "Commons.h"
#include "Helpers.h"

namespace lbm {

template <class U, unsigned int Size>
class StaticArray {
 public:
  U sArray[Size];

  LBM_DEVICE LBM_HOST StaticArray<U, Size>& operator=(
      const StaticArray<U, Size> other) {
    for (auto i = 0; i < Size; ++i) {
      sArray[i] = other[i];
    }

    return *this;
  }

  LBM_DEVICE LBM_HOST StaticArray<U, Size>& operator=(const U other[Size]) {
    for (auto i = 0; i < Size; ++i) {
      sArray[i] = other[i];
    }

    return *this;
  }

  LBM_DEVICE LBM_HOST U& operator[](int i) { return sArray[i]; }

  LBM_DEVICE LBM_HOST const U& operator[](int i) const { return sArray[i]; }

  LBM_DEVICE LBM_HOST U* data() { return sArray; }

  LBM_DEVICE LBM_HOST const U* data() const { return sArray; }

  LBM_DEVICE LBM_HOST static constexpr unsigned int size() { return Size; }
};

template <class U, unsigned int Size>
std::ostream& operator<<(std::ostream& os, const StaticArray<U, Size>& sArray) {
  os << "[";
  for (auto i = 0; i < Size - 1; ++i) {
    os << sArray[i] << " ";
  }
  os << sArray[Size - 1] << "]";
  return os;
}

template <class U, unsigned int Size>
LBM_DEVICE LBM_HOST bool operator==(StaticArray<U, Size> const& lhs,
                                    StaticArray<U, Size> const& rhs) {
  for (auto i = 0; i < Size; ++i) {
    if (!(lhs[i] == rhs[i])) {
      return false;
    }
  }
  return true;
}

}  // namespace lbm
