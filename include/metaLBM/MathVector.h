#ifndef MATHVECTOR_H
#define MATHVECTOR_H

#include <cmath>
#include <fstream>
#include <sstream>

#include "StaticArray.h"
#include "Helpers.h"

namespace lbm {

  template<class U, unsigned int NumberComponents>
  class MathVector {
  public:
    StaticArray<U, NumberComponents> sArray;

    #pragma omp declare simd
    DEVICE HOST
    U& operator[] (int i) {
      return sArray[i];
    }

    #pragma omp declare simd
    DEVICE HOST
    const U& operator[] (int i) const {
      return sArray[i];
    }

    DEVICE HOST
    U * RESTRICT data() {
      return sArray.data();
    }

    DEVICE HOST
    const U * RESTRICT data() const {
      return sArray.data();
    }

    #pragma omp declare simd
    DEVICE HOST
    MathVector<U, NumberComponents>& operator=(const MathVector<U, NumberComponents> other){
      sArray = other.sArray;
      return *this;
    }

    #pragma omp declare simd
    DEVICE HOST
    MathVector<U, NumberComponents>& operator=(const U other[NumberComponents]){
      sArray = other;
      return *this;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline U sum() {
      U sumR = 0;
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
          sumR += sArray[iC];
      }

      return sumR;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline U dot(const MathVector<U, NumberComponents>& other){
      U dotR = 0;
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        dotR += sArray[iC]*other[iC];
      }

      return dotR;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline U norm2() const {
      U norm2R = 0;
      for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        norm2R += sArray[iC]*sArray[iC];
      }

      return norm2R;
    }

  };

  template<class U, unsigned int NumberComponents>
  HOST
  std::ostream& operator<<(std::ostream& os,
                           const MathVector<U, NumberComponents>& mV) {
    os << mV.sArray;
    return os;
  }

  template<class U, unsigned int NumberComponents>
  HOST
  std::ofstream& operator<<(std::ofstream& file,
                           const MathVector<U, NumberComponents>& mV){
    file << "\t\t\t\t";

    UnrolledFor<0, NumberComponents-1>::Do([&] HOST (int i) {
        file << mV[i] << " ";
    });

    file << mV[NumberComponents-1];

    return file;
  }


  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  bool operator==(MathVector<U, NumberComponents> const &lhs,
                  MathVector<U, NumberComponents> const &rhs) {
    return lhs.sArray == rhs.sArray;
  }


  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents>& operator+=(MathVector<U, NumberComponents>& lhs,
                                              const MathVector<U, NumberComponents>& rhs)
  {
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        lhs[iC] += rhs[iC];
    }

    return lhs;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents> operator+(const MathVector<U, NumberComponents>& mV_a,
                                            const MathVector<U, NumberComponents>& mV_b)
  {
    MathVector<U, NumberComponents> mV_result;
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
      mV_result[iC] = mV_a[iC] + mV_b[iC];
    }
    return mV_result;
  }


  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents>& operator*=(MathVector<U, NumberComponents>& mV,
                                              const U factor) {
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
      mV[iC] *= factor;
    }

    return mV;
  }

  template<class U, class V, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents> operator*(const MathVector<U, NumberComponents>& mV,
                                            const V factor)
  {
    MathVector<U, NumberComponents> mV_result;
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        mV_result[iC] = mV[iC] * (U) factor;
    }

    return mV_result;
  }

  template<class U, class V, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents> operator*(const V factor,
                                            const MathVector<U, NumberComponents>& mV)
  {
    MathVector<U, NumberComponents> mV_result;
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        mV_result[iC] = mV[iC] * (U) factor;
    }

    return mV_result;
  }


  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents>& operator/=(MathVector<U, NumberComponents>& mV,
                                              const U factor)
  {
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
      mV[iC] /= factor;
    }

    return mV;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents> operator/(const MathVector<U, NumberComponents>& mV,
                                            const U factor)
  {
    MathVector<U, NumberComponents> mV_result;
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
      mV_result[iC] = mV[iC] / factor;
    }

    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents>& operator-=(MathVector<U, NumberComponents>& lhs,
                                              const MathVector<U, NumberComponents>& rhs)
  {
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
      lhs[iC] -= rhs[iC];
    }

    return lhs;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents> operator-(const MathVector<U, NumberComponents>& mV_a,
                                            const MathVector<U, NumberComponents>& mV_b)
  {
    MathVector<U, NumberComponents> mV_result;
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
      mV_result[iC] = mV_a[iC] - mV_b[iC];
    }

    return mV_result;
  }

  template<unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, NumberComponents> operator-(const MathVector<unsigned int, NumberComponents>& mV_a,
                                                       const MathVector<unsigned int, NumberComponents>& mV_b)
  {
    MathVector<unsigned int, NumberComponents> mV_result;
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
      mV_result[iC] = mV_a[iC] - mV_b[iC];
    }

    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, NumberComponents> operator-(const MathVector<unsigned int, NumberComponents>& mV_a,
                                                       const MathVector<U, NumberComponents>& mV_b)
  {
    MathVector<unsigned int, NumberComponents> mV_result;
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
      mV_result[iC] = mV_a[iC] - (unsigned int) mV_b[iC];
    }

    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, NumberComponents> operator-(const MathVector<U, NumberComponents>& mV_a,
                                                       const MathVector<unsigned int, NumberComponents>& mV_b)
  {
    MathVector<unsigned int, NumberComponents> mV_result;
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
      mV_result[iC] = (unsigned int) mV_a[iC] - mV_b[iC];
    }

    return mV_result;
  }

  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, 3> operator-(const MathVector<unsigned int, 3>& mV_a,
                                        const MathVector<unsigned int, 3>& mV_b)
  {
    return {mV_a[0] - mV_b[0], mV_a[1] - mV_b[1], mV_a[2] - mV_b[2]};
  }

  template<unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, 3> operator-(const MathVector<unsigned int, 3>& mV_a,
                                        const MathVector<unsigned int, NumberComponents>& mV_b)
  {
    MathVector<unsigned int, 3> mV_result = mV_a;
    for(unsigned int iC = 0; iC < NumberComponents; ++iC) {
        mV_result[iC] = mV_a[iC] - mV_b[iC];
    }

    return mV_result;
  }

  template<unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, 3> operator-(const MathVector<unsigned int, NumberComponents>& mV_a,
                                        const MathVector<unsigned int, 3>& mV_b)
  {
    MathVector<unsigned int, 3> mV_result = -1* mV_b;
    for(unsigned int iC = 0; iC < 3; ++iC) {
      mV_result[iC] = mV_a[iC] - mV_b[iC];
    }

    return mV_result;
  }


  template<class T, unsigned int Dimension>
  struct Project {
    #pragma omp declare simd
    DEVICE HOST
    static inline MathVector<T, Dimension> Do(const MathVector<T, 3>& mV) {

      MathVector<T, Dimension> mVProjected{{ (T) 0 }};

      for(unsigned int iD = 0; iD < Dimension; ++iD) {
        mVProjected[iD] = mV[iD];
      }

    return mVProjected;
    }
  };

  template<class T, unsigned int Dimension>
  struct ProjectAndLeave1 {
    #pragma omp declare simd
    DEVICE HOST
    static inline MathVector<T, 3> Do(const MathVector<T, 3>& mV) {

      MathVector<T, 3> mVProjected = { (T) 1, (T) 1, (T) 1};

      for(unsigned int iD = 0; iD < Dimension; ++iD) {
          mVProjected[iD] = mV[iD];
      }

    return mVProjected;
    }
  };


}

#endif // MATHVECTOR_H
