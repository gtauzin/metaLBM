#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>

namespace lbm {

  template<int Begin, int End, int Step = 1>
    struct UnrolledFor {
      template<typename F>
      static void Do (F f) {
        f(Begin);
        UnrolledFor<Begin+Step, End, Step>::Do(f);
      }
    };

  template<int End>
    struct UnrolledFor<End, End> {
    template<typename F>
      static void Do (F f) {
    }
  };


  template<class U, unsigned int size>
    class MathVector {
  public:
    U mV[size];

    U sum() {
      U sumR = 0;
      UnrolledFor<0, size>::Do([&] (int i) {
          sumR += mV[i];
        });

      return sumR;
    }

    U norm2() const {
      U norm2R = 0;
      UnrolledFor<0, size>::Do([&] (int i) {
          norm2R += mV[i] * mV[i];
        });

      return norm2R;
    }

    U dot(MathVector<U, size> other){
      U dotR = 0;

      UnrolledFor<0, size>::Do([&] (int i) {
          dotR += mV[i]*other[i];
        });

      return dotR;
    }

    U magnitude(){
      return sqrt(this->norm2());
    }

    MathVector<U, size>& operator=(const MathVector<U, size>& other){
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] = other[i];
        });

      return *this;
    }

    MathVector<U, size>& operator=(const U other[size]){
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] = other[i];
        });

      return *this;
    }

    U& operator[] (int i) {
      return mV[i];
    }

    const U& operator[] (int i) const {
      return mV[i];
    }

  };


  template<class U, unsigned int size>
    std::ostream& operator<<(std::ostream& os, const MathVector<U, size>& mV){
    os << "[ ";
    UnrolledFor<0, size>::Do([&] (int i) {
        os << mV[i] << " ";
      });
    os << "]\n";
    return os;
  }

  template<class U, unsigned int size>
    MathVector<U, size>& operator+=(MathVector<U, size>& mV_a, const MathVector<U, size>& mV_b)
    {
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_a[i] += mV_b[i];
        });

      return mV_a;
    }

  template<class U, unsigned int size>
    MathVector<U, size> operator+(const MathVector<U, size>& mV_a, const MathVector<U, size>& mV_b)
    {
      MathVector<U, size> mV_c;
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_c[i] = mV_a[i] + mV_b[i];
        });
      return mV_c;
    }

    template<class U, unsigned int size>
    MathVector<int, 3> operator+(const MathVector<int, 3>& mV_lhs,
                                 const MathVector<U, size>& mV_rhs)
    {
      MathVector<int, 3> mV_result{{0}};
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_result[i] = mV_lhs[i] + (int) mV_rhs[i];
        });
      return mV_result;
    }

    template<class U, unsigned int size>
    MathVector<int, 3> operator+(const MathVector<U, size>& mV_lhs,
                                 const MathVector<int, 3>& mV_rhs)
    {
      MathVector<int, 3> mV_result{{0}};
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_result[i] = (int) mV_lhs[i] + mV_rhs[i];
        });
      return mV_result;
    }

        template<class U, unsigned int size>
    MathVector<int, 3> operator-(const MathVector<int, 3>& mV_lhs,
                                 const MathVector<U, size>& mV_rhs)
    {
      MathVector<int, 3> mV_result{{0}};
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_result[i] = mV_lhs[i] - (int) mV_rhs[i];
        });
      return mV_result;
    }

    template<class U, unsigned int size>
    MathVector<int, 3> operator-(const MathVector<U, size>& mV_lhs,
                                 const MathVector<int, 3>& mV_rhs)
    {
      MathVector<int, 3> mV_result{{0}};
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_result[i] = (int) mV_lhs[i] - mV_rhs[i];
        });
      return mV_result;
    }

  template<class U, unsigned int size>
    MathVector<U, size>& operator*=(MathVector<U, size>& mV, const U factor)
    {
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] *= factor;
        });

      return mV;
    }

  template<class U, unsigned int size>
    MathVector<U, size> operator*(const MathVector<U, size>& mV_in, const U factor)
    {
      MathVector<U, size> mV;
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] = mV_in[i] * factor;
        });

      return mV;
    }

  template<class U, unsigned int size>
    MathVector<U, size> operator*(const U factor, const MathVector<U, size>& mV_in)
    {
      MathVector<U, size> mV;
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] = factor * mV_in[i];
        });

      return mV;
    }


  template<class U, unsigned int size>
    MathVector<U, size>& operator/=(MathVector<U, size>& mV, const U factor)
    {
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] /= factor;
        });

      return mV;
    }

  template<class U, unsigned int size>
    MathVector<U, size> operator/(const MathVector<U, size>& mV_in, const U factor)
    {
      MathVector<U, size> mV;
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] = mV_in[i] / factor;
        });

      return mV;
    }

  template<class U, unsigned int size>
    MathVector<U, size>& operator-=(MathVector<U, size>& mV_a, const MathVector<U, size>& mV_b)
    {
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_a[i] -= mV_b[i];
        });

      return mV_a;
    }

  template<class U, unsigned int size>
    MathVector<U, size> operator-(const MathVector<U, size>& mV_a, const MathVector<U, size>& mV_b)
    {
      MathVector<U, size> mV_c;
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_c[i] = mV_a[i] - mV_b[i];
        });

      return mV_c;
    }

  template<class U, unsigned int size>
    bool operator==(MathVector<U, size> const &lhs,
                    MathVector<U, size> const &rhs) {
    for(unsigned int i = 0; i < size; ++i){
      if (!(lhs[i] == rhs[i])) {
        return false;
      }
    }
    return true;
  }

}
#endif // STRUCTURE_H
