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


  template<class T, unsigned int size>
    class MathVector {
  public:
    T mV[size];


    T sum() {
      T sumR = 0;
      UnrolledFor<0, size>::Do([&] (int i) {
          sumR += mV[i];
        });

      return sumR;
    }

    T norm2() const {
      T norm2R = 0;
      UnrolledFor<0, size>::Do([&] (int i) {
          norm2R += mV[i] * mV[i];
        });

      return norm2R;
    }

    T dot(MathVector<T, size> other){
      T dotR = 0;

      UnrolledFor<0, size>::Do([&] (int i) {
          dotR += mV[i]*other[i];
        });

      return dotR;
    }

    T magnitude(){
      return sqrt(this->norm2());
    }

    MathVector<T, size>& operator=(const MathVector<T, size>& other){
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] = other[i];
        });

      return *this;
    }

    MathVector<T, size>& operator=(const T other[size]){
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] = other[i];
        });

      return *this;
    }

    T& operator[] (int i) {
      return mV[i];
    }

    const T& operator[] (int i) const {
      return mV[i];
    }

  };


  template<class T, unsigned int size>
    std::ostream& operator<<(std::ostream& os, const MathVector<T, size>& mV){
    os << "[ ";
    UnrolledFor<0, size>::Do([&] (int i) {
        os << mV[i] << " ";
      });
    os << "]\n";
    return os;
  }

  template<class T, unsigned int size>
    MathVector<T, size>& operator+=(MathVector<T, size>& mV_a, const MathVector<T, size>& mV_b)
    {
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_a[i] += mV_b[i];
        });

      return mV_a;
    }

  template<class T, unsigned int size>
    MathVector<T, size> operator+(const MathVector<T, size>& mV_a, const MathVector<T, size>& mV_b)
    {
      MathVector<T, size> mV_c;
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_c[i] = mV_a[i] + mV_b[i];
        });
      return mV_c;
    }

  template<class T, unsigned int size>
    MathVector<T, size>& operator*=(MathVector<T, size>& mV, const T factor)
    {
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] *= factor;
        });

      return mV;
    }

  template<class T, unsigned int size>
    MathVector<T, size> operator*(const MathVector<T, size>& mV_in, const T factor)
    {
      MathVector<T, size> mV;
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] = mV_in[i] * factor;
        });

      return mV;
    }

  template<class T, unsigned int size>
    MathVector<T, size> operator*(const T factor, const MathVector<T, size>& mV_in)
    {
      MathVector<T, size> mV;
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] = factor * mV_in[i];
        });

      return mV;
    }


  template<class T, unsigned int size>
    MathVector<T, size>& operator/=(MathVector<T, size>& mV, const T factor)
    {
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] /= factor;
        });

      return mV;
    }

  template<class T, unsigned int size>
    MathVector<T, size> operator/(const MathVector<T, size>& mV_in, const T factor)
    {
      MathVector<T, size> mV;
      UnrolledFor<0, size>::Do([&] (int i) {
          mV[i] = mV_in[i] / factor;
        });

      return mV;
    }

  template<class T, unsigned int size>
    MathVector<T, size>& operator-=(MathVector<T, size>& mV_a, const MathVector<T, size>& mV_b)
    {
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_a[i] -= mV_b[i];
        });

      return mV_a;
    }

  template<class T, unsigned int size>
    MathVector<T, size> operator-(const MathVector<T, size>& mV_a, const MathVector<T, size>& mV_b)
    {
      MathVector<T, size> mV_c;
      UnrolledFor<0, size>::Do([&] (int i) {
          mV_c[i] = mV_a[i] - mV_b[i];
        });

      return mV_c;
    }

  template<class T, unsigned int size>
    bool operator==(MathVector<T, size> const &lhs,
                    MathVector<T, size> const &rhs) {
    for(unsigned int i = 0; i < size; ++i){
      if (lhs[i] != rhs[i]) {
        return false;
      }
    }
    return true;
  }

}
#endif // STRUCTURE_H
