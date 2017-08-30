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
      UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
          sumR += sArray[i];
        });

      return sumR;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline U norm2() const {
      U norm2R = 0;
      UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
          norm2R += sArray[i] * sArray[i];
        });

      return norm2R;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline U dot(MathVector<U, NumberComponents> other){
      U dotR = 0;

      UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
          dotR += sArray[i]*other[i];
        });

      return dotR;
    }

    #pragma omp declare simd
    DEVICE HOST
    inline U magnitude(){
      return sqrt(this->norm2());
    }

  };

  template<class U, unsigned int NumberComponents>
  std::ostream& operator<<(std::ostream& os,
                           const MathVector<U, NumberComponents>& mV) {
    os << mV.sArray;
    return os;
  }

  template<class U, unsigned int NumberComponents>
  std::ofstream& operator<<(std::ofstream& file,
                           const MathVector<U, NumberComponents>& mV){
    file << "\t\t\t\t";

    UnrolledFor<0, NumberComponents-1>::Do([&] HOST DEVICE (int i) {
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
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        lhs[i] += rhs[i];
      });

    return lhs;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents> operator+(const MathVector<U, NumberComponents>& mV_a,
                                            const MathVector<U, NumberComponents>& mV_b)
  {
    MathVector<U, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV_result[i] = mV_a[i] + mV_b[i];
      });
    return mV_result;
  }


  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents>& operator*=(MathVector<U, NumberComponents>& mV,
                                              const U factor)
  {
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV[i] *= factor;
      });

    return mV;
  }

  template<class U, class V, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents> operator*(const MathVector<U, NumberComponents>& mV,
                                            const V factor)
  {
    MathVector<U, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV_result[i] = mV[i] * (U) factor;
      });

    return mV_result;
  }

  template<class U, class V, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents> operator*(const V factor,
                                            const MathVector<U, NumberComponents>& mV)
  {
    MathVector<U, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV_result[i] = mV[i] * (U) factor;
      });

    return mV_result;
  }


  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents>& operator/=(MathVector<U, NumberComponents>& mV,
                                              const U factor)
  {
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV[i] /= factor;
      });

    return mV;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents> operator/(const MathVector<U, NumberComponents>& mV,
                                            const U factor)
  {
    MathVector<U, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV_result[i] = mV[i] / factor;
      });

    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents>& operator-=(MathVector<U, NumberComponents>& lhs,
                                              const MathVector<U, NumberComponents>& rhs)
  {
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        lhs[i] -= rhs[i];
      });

    return lhs;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<U, NumberComponents> operator-(const MathVector<U, NumberComponents>& mV_a,
                                            const MathVector<U, NumberComponents>& mV_b)
  {
    MathVector<U, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV_result[i] = mV_a[i] - mV_b[i];
      });

    return mV_result;
  }

  template<unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, NumberComponents> operator-(const MathVector<unsigned int, NumberComponents>& mV_a,
                                                       const MathVector<unsigned int, NumberComponents>& mV_b)
  {
    MathVector<unsigned int, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV_result[i] = mV_a[i] - mV_b[i];
      });

    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, NumberComponents> operator-(const MathVector<unsigned int, NumberComponents>& mV_a,
                                                       const MathVector<U, NumberComponents>& mV_b)
  {
    MathVector<unsigned int, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV_result[i] = mV_a[i] - (unsigned int) mV_b[i];
      });

    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, NumberComponents> operator-(const MathVector<U, NumberComponents>& mV_a,
                                                       const MathVector<unsigned int, NumberComponents>& mV_b)
  {
    MathVector<unsigned int, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV_result[i] = (unsigned int) mV_a[i] - mV_b[i];
      });

    return mV_result;
  }

  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, 3> operator-(const MathVector<unsigned int, 3>& mV_a,
                                        const MathVector<unsigned int, 3>& mV_b)
  {
    MathVector<unsigned int, 3> mV_result = mV_a;
    UnrolledFor<0, 3>::Do([&] HOST DEVICE (int i) {
        mV_result[i] = mV_a[i] - mV_b[i];
      });

    return mV_result;
  }

  template<unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, 3> operator-(const MathVector<unsigned int, 3>& mV_a,
                                        const MathVector<unsigned int, NumberComponents>& mV_b)
  {
    MathVector<unsigned int, 3> mV_result = mV_a;
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV_result[i] = mV_a[i] - mV_b[i];
      });

    return mV_result;
  }

  template<unsigned int NumberComponents>
  #pragma omp declare simd
  DEVICE HOST
  MathVector<unsigned int, 3> operator-(const MathVector<unsigned int, NumberComponents>& mV_a,
                                        const MathVector<unsigned int, 3>& mV_b)
  {
    MathVector<unsigned int, 3> mV_result = -1* mV_b;
    UnrolledFor<0, NumberComponents>::Do([&] HOST DEVICE (int i) {
        mV_result[i] = mV_a[i] - mV_b[i];
      });

    return mV_result;
  }


  template<class T, unsigned int Dimension>
  struct Project {
    #pragma omp declare simd
    DEVICE HOST
    static inline MathVector<T, Dimension> Do(const MathVector<T, 3>& mV) {

      MathVector<T, Dimension> mVProjected{{ (T) 0 }};

    UnrolledFor<0, Dimension>::Do([&] HOST DEVICE (unsigned int iD) {
      mVProjected[iD] = mV[iD];
     });

    return mVProjected;
    }
  };

  template<class T, unsigned int Dimension>
  struct ProjectAndLeave1 {
    #pragma omp declare simd
    DEVICE HOST
    static inline MathVector<T, 3> Do(const MathVector<T, 3>& mV) {

      MathVector<T, 3> mVProjected = { (T) 1, (T) 1, (T) 1};

      UnrolledFor<0, Dimension>::Do([&] HOST DEVICE (unsigned int iD) {
          mVProjected[iD] = mV[iD];
      });

    return mVProjected;
    }
  };


}

#endif // MATHVECTOR_H
