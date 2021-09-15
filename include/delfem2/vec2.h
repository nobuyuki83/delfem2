/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_VEC2_H
#define DFM2_VEC2_H

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>

#include "delfem2/dfm2_inline.h"

#if defined(_MSC_VER)
#  pragma warning( push )
// nonstandard extension used : nameless struct/union
#  pragma warning( disable : 4201 )
#endif

// -----------------------------------------------------

namespace delfem2 {

template<typename T>
T Dot2(const T w[2], const T v[2]);

template<typename T>
DFM2_INLINE T Length2(const T v[2]);

template<typename T>
T Distance2(const T v1[2], const T v2[2]);

template<typename T>
T SquareLength2(const T v[2]);

template<typename T>
T SquareDistance2(const T v1[2], const T v2[2]);

template<typename T>
void Normalize2(T w[2]);

// ----------------------------------------------------

template<typename T>
T Area_Tri2(
    const T v1[2],
    const T v2[2],
    const T v3[2]);

DFM2_INLINE bool IsCrossLines(
    const double po_s0[2],
    const double po_e0[2],
    const double po_s1[2],
    const double po_e1[2]);

template<typename T>
DFM2_INLINE void GaussianDistribution2(
    T noise[2]);

// ----------------------

template<typename T>
class CVec2;

template<typename T, typename T1>
CVec2<T> operator*(T1, const CVec2<T> &);

template<typename T, typename T1>
CVec2<T> operator*(const CVec2<T> &, T1);

template<typename T>
std::ostream &operator<<(std::ostream &output, const CVec2<T> &v);

template<typename T>
std::istream &operator>>(std::istream &input, CVec2<T> &v);

template<typename T>
T operator^(const CVec2<T> &lhs, const CVec2<T> &rhs);

template<typename T>
CVec2<T> operator/(const CVec2<T> &vec, double d); //! divide by real number

/**
 * @brief 2 dimensional vector class
 */
template<typename T>
class CVec2 {
 public:
  CVec2() : p{0.0, 0.0} {}
  CVec2(const CVec2 &) = default;
  CVec2(T x, T y) {
    this->p[0] = x;
    this->p[1] = y;
  }
  CVec2(const std::array<T,2> &arr) {
    this->p[0] = arr[0];
    this->p[1] = arr[1];
  }
  explicit CVec2(T x) {
    this->p[0] = x;
    this->p[1] = x;
  }
  // above: constructor / destructor
  // -------------------------------
  // below: operator
  CVec2 operator-() const {
    return CVec2(-p[0], -p[1]);
  }
  inline CVec2 &operator+=(const CVec2 &rhs) {
    p[0] += rhs.p[0];
    p[1] += rhs.p[1];
    return *this;
  }
  inline CVec2 &operator-=(const CVec2 &rhs) {
    p[0] -= rhs.p[0];
    p[1] -= rhs.p[1];
    return *this;
  }
  inline CVec2 &operator*=(double scale) {
    p[0] *= scale;
    p[1] *= scale;
    return *this;
  }
  inline CVec2 &operator/=(double d) {
    if (fabs(d) < 1.0e-6) {
      assert(0);
      return *this;
    }
    p[0] /= d;
    p[1] /= d;
    return *this;
  }
  inline CVec2 operator+(const CVec2 &rhs) const {
    CVec2 v = *this;
    v += rhs;
    return v;
  }
  inline CVec2 operator-(const CVec2 &rhs) const {
    CVec2 v = *this;
    v -= rhs;
    return v;
  }
  inline T operator[](int i) const {
    assert(i < 2);
    return p[i];
  }
  inline T &operator[](int i) {
    assert(i < 2);
    return p[i];
  }
  // above: operator
  // ---------------
  // below: function

  //! @details (named after Eigen and STL)
  T *data() { return p; }

  //! @details (named after Eigen and STL)
  const T *data() const { return p; }

  //! @brief squared Euclidian norm (named similalty to Eigen)
  T squaredNorm() const {
    return p[0] * p[0] + p[1] * p[1];
  }

  //! @brief in place normalization with Euclidian norm (named similarly to Eigen)
  void normalize() {
    const double mag = norm();
    p[0] /= mag;
    p[1] /= mag;
  }

  //! @brief Euclidian norm (named similalty to Eigen)
  T norm() const{
    return std::sqrt(p[0]*p[0]+p[1]*p[1]);
  }
  
  //! @brief normalizeation (named similarly to Eigen)
  CVec2 normalized() const {
    CVec2 r(*this);
    r.normalize();
    return r;
  }

  //! @brief set zero vector (named similarly to Eigen)
  inline void setZero() {
    p[0] = T(0);
    p[1] = T(0);
  }

  //! @brief dot product (named similarly to Eigen)
  T dot(const CVec2<T>& rhs) const {
    return x*rhs.x+y*rhs.y;
  }

  template <typename S>
  CVec2<S> cast() const {
    return CVec2<S>(
        static_cast<S>(p[0]),
        static_cast<S>(p[1]));
  };

  CVec2 Mat3Vec2_AffineProjection(const T *A) const {
    CVec2<T> y0;
    y0.p[0] = A[0] * p[0] + A[1] * p[1] + A[2];
    y0.p[1] = A[3] * p[0] + A[4] * p[1] + A[5];
    const T w = A[6] * p[0] + A[7] * p[1] + A[8];
    y0.p[0] /= w;
    y0.p[1] /= w;
    return y0;
  }
  CVec2 Mat3Vec2_AffineDirection(const T *A) const {
    CVec2<T> y0;
    y0.p[0] = A[0] * p[0] + A[1] * p[1];
    y0.p[1] = A[3] * p[0] + A[4] * p[1];
    return y0;
  }

  CVec2 Rotate(T t) const {
    CVec2 r;
    r.p[0] = cos(t) * p[0] - sin(t) * p[1];
    r.p[1] = sin(t) * p[0] + cos(t) * p[1];
    return r;
  }

 public:
  union {
    T p[2];
    struct {
      T x, y;
    };
  };
};
using CVec2d = CVec2<double>;
using CVec2f = CVec2<float>;
using CVec2i = CVec2<int>;

template<typename T>
CVec2<T> rotate(const CVec2<T> &p0, double theta);

template<typename T>
CVec2<T> rotate90(const CVec2<T> &p0);

// --------------------------------------------

template<typename T>
CVec2<T> Mat2Vec(const double A[4], const CVec2<T> &v);

//! @brief Area of the Triangle
template<typename T>
double Area_Tri(
    const CVec2<T> &v1,
    const CVec2<T> &v2,
    const CVec2<T> &v3);

template<typename T>
double Cross(
    const CVec2<T> &v1,
    const CVec2<T> &v2);

template<typename T>
double SquareLength(
    const CVec2<T> &point);

template<typename T>
double Length(
    const CVec2<T> &point);

//! @brief Length between two points
template<typename T>
T Distance(
    const CVec2<T> &ipo0,
    const CVec2<T> &ipo1);

//! @brief Length between two points
template<typename T>
double SquareDistance(
    const CVec2<T> &ipo0,
    const CVec2<T> &ipo1);

//! @brief Hight of a triangle : between v1 and line of v2-v3
template<typename T>
double TriHeight(
    const CVec2<T> &v1,
    const CVec2<T> &v2,
    const CVec2<T> &v3);

//! @brief compute dot product
template<typename T>
double Dot(const CVec2<T> &ipo0, const CVec2<T> &ipo1);

//! @details get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
template<typename T>
double FindNearestPointParameter_Line_Point(
    const CVec2<T> &po_c,
    const CVec2<T> &po_s,
    const CVec2<T> &po_e);

//! @details  get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
template<typename T>
CVec2<T> GetNearest_LineSeg_Point(
    const CVec2<T> &po_c,
    const CVec2<T> &po_s,
    const CVec2<T> &po_e);

//! @details  get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
template<typename T>
double GetDist_LineSeg_Point(
    const CVec2<T> &po_c,
    const CVec2<T> &po_s,
    const CVec2<T> &po_e);

template<typename T>
bool IsCross_LineSeg_LineSeg(
    const CVec2<T> &po_s0,
    const CVec2<T> &po_e0,
    const CVec2<T> &po_s1,
    const CVec2<T> &po_e1);

template<typename T>
double GetDist_LineSeg_LineSeg(
    const CVec2<T> &po_s0,
    const CVec2<T> &po_e0,
    const CVec2<T> &po_s1,
    const CVec2<T> &po_e1);

/**
 * @brief square root of circumradius
 */
template<typename T>
double SquareCircumradius(
    const CVec2<T> &p0,
    const CVec2<T> &p1,
    const CVec2<T> &p2);

/**
 * @brief center of the circumcircle
 */
template<typename T>
bool CenterCircumcircle(
    const CVec2<T> &p0,
    const CVec2<T> &p1,
    const CVec2<T> &p2,
    CVec2<T> &center);

/**
 * @brief check if Delaunay condition satisfied
 * @return
 * 0 : p3 is inside circum circle on the p0,p1,p2
 * 1 :       on
 * 2 :       outsdie
 */
template<typename T>
int DetDelaunay(
    const CVec2<T> &p0,
    const CVec2<T> &p1,
    const CVec2<T> &p2,
    const CVec2<T> &p3);


// ---------------------------------------------------------------


//! @brief Area of the Triangle (3 indexes and vertex array)
template<typename T>
double Area_Tri(
    int iv1,
    int iv2,
    int iv3,
    const std::vector<CVec2<T> > &point);

template<typename T>
void MakeMassMatrixTri(
    double M[9],
    double rho,
    const unsigned int aIP[3],
    const std::vector<CVec2<T> > &aVec2);

} // namespace delfem2

#if defined(_MSC_VER)
#  pragma warning( pop )
#endif

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/vec2.cpp"
#endif

#endif // VEC_2


