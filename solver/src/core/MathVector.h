#ifndef MATHVECTOR_H
#define MATHVECTOR_H

#include <cmath>

#include "StaticArray.h"
#include "Helpers.h"

namespace lbm {

  template<class U, unsigned int NumberComponents>
  class MathVector {
  public:
    StaticArray<U, NumberComponents> sArray;

    U& operator[] (int i) {
      return sArray[i];
    }

    const U& operator[] (int i) const {
      return sArray[i];
    }

    MathVector<U, NumberComponents>& operator=(const MathVector<U, NumberComponents> other){
      sArray = other.sArray;
      return *this;
    }

    MathVector<U, NumberComponents>& operator=(const U other[NumberComponents]){
      sArray = other;
      return *this;
    }

    inline U sum() {
      U sumR = 0;
      UnrolledFor<0, NumberComponents>::Do([&] (int i) {
          sumR += sArray[i];
        });

      return sumR;
    }

    inline U norm2() const {
      U norm2R = 0;
      UnrolledFor<0, NumberComponents>::Do([&] (int i) {
          norm2R += sArray[i] * sArray[i];
        });

      return norm2R;
    }

    inline U dot(MathVector<U, NumberComponents> other){
      U dotR = 0;

      UnrolledFor<0, NumberComponents>::Do([&] (int i) {
          dotR += sArray[i]*other[i];
        });

      return dotR;
    }

    inline U magnitude(){
      return sqrt(this->norm2());
    }

  };

  template<class U, unsigned int NumberComponents>
  std::ostream& operator<<(std::ostream& os,
                           const MathVector<U, NumberComponents>& mV){
    os << mV.sArray;
    return os;
  }

  template<class U, unsigned int NumberComponents>
  std::ofstream& operator<<(std::ofstream& file,
                           const MathVector<U, NumberComponents>& mV){
    file << "\t\t\t\t";

    UnrolledFor<0, NumberComponents-1>::Do([&] (int i) {
        file << mV[i] << " ";
    });

    file << mV[NumberComponents-1] << "\n";

    return file;
  }


  template<class U, unsigned int NumberComponents>
  bool operator==(MathVector<U, NumberComponents> const &lhs,
                  MathVector<U, NumberComponents> const &rhs) {
    return lhs.sArray == rhs.sArray;
  }


  template<class U, unsigned int NumberComponents>
  MathVector<U, NumberComponents>& operator+=(MathVector<U, NumberComponents>& lhs,
                                              const MathVector<U, NumberComponents>& rhs)
  {
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        lhs[i] += rhs[i];
      });

    return lhs;
  }

  template<class U, unsigned int NumberComponents>
  MathVector<U, NumberComponents> operator+(const MathVector<U, NumberComponents>& mV_a,
                                            const MathVector<U, NumberComponents>& mV_b)
  {
    MathVector<U, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        mV_result[i] = mV_a[i] + mV_b[i];
      });
    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  MathVector<int, 3> operator+(const MathVector<int, 3>& mV_a,
                               const MathVector<U, NumberComponents>& mV_b)
  {
    MathVector<int, 3> mV_result{{0}};
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        mV_result[i] = mV_a[i] + (int) mV_b[i];
      });
    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  MathVector<int, 3> operator+(const MathVector<U, NumberComponents>& mV_a,
                               const MathVector<int, 3>& mV_b)
  {
    MathVector<int, 3> mV_result{{0}};
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        mV_result[i] = (int) mV_a[i] + mV_b[i];
      });
    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  MathVector<int, 3> operator-(const MathVector<int, 3>& mV_a,
                               const MathVector<U, NumberComponents>& mV_b)
  {
    MathVector<int, 3> mV_result{{0}};
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        mV_result[i] = mV_a[i] - (int) mV_b[i];
      });
    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  MathVector<int, 3> operator-(const MathVector<U, NumberComponents>& mV_a,
                               const MathVector<int, 3>& mV_b)
  {
    MathVector<int, 3> mV_result{{0}};
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        mV_result[i] = (int) mV_a[i] - mV_b[i];
      });
    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  MathVector<U, NumberComponents>& operator*=(MathVector<U, NumberComponents>& mV,
                                              const U factor)
  {
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        mV[i] *= factor;
      });

    return mV;
  }

  template<class U, unsigned int NumberComponents>
  MathVector<U, NumberComponents> operator*(const MathVector<U, NumberComponents>& mV,
                                            const U factor)
  {
    MathVector<U, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        mV_result[i] = mV[i] * factor;
      });

    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  MathVector<U, NumberComponents> operator*(const U factor,
                                            const MathVector<U, NumberComponents>& mV)
  {
    MathVector<U, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        mV_result[i] = factor * mV[i];
      });

    return mV_result;
  }


  template<class U, unsigned int NumberComponents>
  MathVector<U, NumberComponents>& operator/=(MathVector<U, NumberComponents>& mV,
                                              const U factor)
  {
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        mV[i] /= factor;
      });

    return mV;
  }

  template<class U, unsigned int NumberComponents>
  MathVector<U, NumberComponents> operator/(const MathVector<U, NumberComponents>& mV,
                                            const U factor)
  {
    MathVector<U, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        mV_result[i] = mV[i] / factor;
      });

    return mV_result;
  }

  template<class U, unsigned int NumberComponents>
  MathVector<U, NumberComponents>& operator-=(MathVector<U, NumberComponents>& lhs,
                                              const MathVector<U, NumberComponents>& rhs)
  {
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        lhs[i] -= rhs[i];
      });

    return lhs;
  }

  template<class U, unsigned int NumberComponents>
  MathVector<U, NumberComponents> operator-(const MathVector<U, NumberComponents>& mV_a,
                                            const MathVector<U, NumberComponents>& mV_b)
  {
    MathVector<U, NumberComponents> mV_result;
    UnrolledFor<0, NumberComponents>::Do([&] (int i) {
        mV_result[i] = mV_a[i] - mV_b[i];
      });

    return mV_result;
  }

}

#endif // MATHVECTOR_H
