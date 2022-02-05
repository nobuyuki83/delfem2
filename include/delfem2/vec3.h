/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @detail The order of dependency in delfem2 is "vec2 -> mat2 -> vec3 -> quaternion -> mat3 -> mat4",
 */

#ifndef DFM2_VEC3_H
#define DFM2_VEC3_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <array>

#include "delfem2/dfm2_inline.h"

#if defined(_MSC_VER)
#  pragma warning( push )
#  pragma warning( disable : 4201 )
#endif

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

template<typename REAL>
void AverageTwo3(
    REAL po[3],
    const REAL p0[3], const REAL p1[3]);

template<typename REAL>
void AverageFour3(
    REAL po[3],
    const REAL p0[3], const REAL p1[3], const REAL p2[3], const REAL p3[3]);

/**
 * @func add values for 3-array (vo += vi)
 * @tparam REAL float and double
 * @param vo (out)
 * @param vi (in)
 */
template<typename REAL>
DFM2_INLINE void Add3(
    REAL vo[3],
    const REAL vi[3]);

// above: no dependency
// -------------------------------------------------------------
// below: definition of CVec3

template<typename T>
class CVec3;

template<typename T0, typename T1>
CVec3<T0> operator*(T1 d, const CVec3<T0> &rhs);

template<typename T0, typename T1>
CVec3<T0> operator/(const CVec3<T0> &vec, T1 d);

template<typename T>
std::ostream &operator<<(std::ostream &output, const CVec3<T> &v);

template<typename T>
std::istream &operator>>(std::istream &input, CVec3<T> &v);

template<typename T>
std::ostream &operator<<(std::ostream &output, const std::vector<CVec3<T> > &aV);

template<typename T>
std::istream &operator>>(std::istream &input, std::vector<CVec3<T> > &aV);

/**
 * @class 3 dimentional vector class
 */
template<typename T>
class CVec3 {
 public:
  CVec3(T vx, T vy, T vz) : p{vx, vy, vz} {}
  CVec3() : p{0.0, 0.0, 0.0} {}
  CVec3(const CVec3 &) = default;

  template<typename S>
  CVec3(const std::array<S, 3> &&v) : p{
      static_cast<T>(v[0]),
      static_cast<T>(v[1]),
      static_cast<T>(v[2])} {}

  template<typename S>
  explicit CVec3(const S *v) : p{
      static_cast<T>(v[0]),
      static_cast<T>(v[1]),
      static_cast<T>(v[2])} {}

  template<typename S>
  explicit CVec3(const std::vector<S> &v)  : p{v[0], v[1], v[2]} {}

  virtual ~CVec3() = default;

  // above: constructor / destructor
  // -------------------------------
  // below: operator

  inline CVec3 &operator=(const CVec3 &b) {
    if (this != &b) {
      x = b.x;
      y = b.y;
      z = b.z;
    }
    return *this;
  }
  inline CVec3 operator-() const { return ((T) (-1)) * (*this); }
  inline CVec3 operator+() const { return *this; }
  inline CVec3 operator+() { return *this; }
  inline CVec3 operator-() { return CVec3(-p[0], -p[1], -p[2]); }
  inline CVec3 &operator+=(const CVec3 &b) {
    x += b.x;
    y += b.y;
    z += b.z;
    return *this;
  }
  inline CVec3 &operator-=(const CVec3 &b) {
    x -= b.x;
    y -= b.y;
    z -= b.z;
    return *this;
  }
  inline CVec3 operator+(const CVec3 &b) const {
    return CVec3(x + b.x, y + b.y, z + b.z);
  }

  inline CVec3 operator-(const CVec3 &b) const {
    return CVec3(x - b.x, y - b.y, z - b.z);
  }

  inline CVec3 operator*(T b) const {
    return CVec3(x * b, y * b, z * b);
  }

  inline CVec3 &operator*=(T d) {
    x *= d;
    y *= d;
    z *= d;
    return *this;
  }
  inline CVec3 &operator/=(T d) {
    if (fabs(d) < NEARLY_ZERO) { return *this; }
    x /= d;
    y /= d;
    z /= d;
    return *this;
  }
  template<typename INDEX>
  inline T operator[](INDEX i) const {
    assert(i < 3);
    return p[i];
  }
  template<typename INDEX>
  inline T &operator[](INDEX i) {
    assert(i < 3);
    return p[i];
  }
  template<typename INDEX>
  inline T operator()(INDEX i) const {
    assert(i < 3);
    return p[i];
  }
  template<typename INDEX>
  inline T &operator()(INDEX i) {
    assert(i < 3);
    return p[i];
  }

  // above: operator
  // ------

  //! @details named after Eigen library
  // [[nodiscard]]
  CVec3 normalized() const {
    CVec3 r = (*this);
    r.normalize();
    return r;
  }

  //! @brief in-place normalization of vector
  //! @details named after Eigen library
  void normalize();

  //! @details named after Eigen library
  inline T norm() const {
    return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
  }

  //! @details named after Eigen library
  inline T squaredNorm() const {
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
  }

  //! @details named after Eigen library
  void setZero();

  //! @details named after Eigen library
  T dot(const CVec3 &rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; }

  /**
   * cross product
   * @param[in] b
   * @return (this) crossproduct b
   * @details named after Eigen library
   */
  CVec3 cross(const CVec3 &b) const {
    return CVec3(y * b.z - z * b.y,
                 z * b.x - x * b.z,
                 x * b.y - y * b.x);
  }

  CVec3 mult(const CVec3 &b) const { return CVec3(x * b.x, y * b.y, z * b.z); }

  template<typename S>
  CVec3<S> cast() const {
    return CVec3<S>(static_cast<S>(p[0]),
                    static_cast<S>(p[1]),
                    static_cast<S>(p[2]));
  };

  T *data() { return p; }
  const T *data() const { return p; }

  [[nodiscard]] std::vector<double> stlvec() const {
    std::vector<double> d = {p[0], p[1], p[2]};
    return d;
  }

  template<typename S>
  void SetVector(S vx, S vy, S vz) {
    p[0] = vx;
    p[1] = vy;
    p[2] = vz;
  }

  template<typename S>
  void CopyTo(S *v) const {
    v[0] = p[0];
    v[1] = p[1];
    v[2] = p[2];
  }

  template<typename S>
  void push_back_to_vector(std::vector<S> &vec) const {
    vec.push_back(p[0]);
    vec.push_back(p[1]);
    vec.push_back(p[2]);
  }

  template<typename S>
  void CopyToScale(S *v, S s) const {
    v[0] = p[0] * s;
    v[1] = p[1] * s;
    v[2] = p[2] * s;
  }

  template<typename S>
  void AddToScale(S *v, S s) const {
    v[0] += p[0] * s;
    v[1] += p[1] * s;
    v[2] += p[2] * s;
  }

  void Print() const {
    std::cout << p[0] << " " << p[1] << " " << p[2] << std::endl;
  }
  [[nodiscard]] bool isNaN() const {
    double s = p[0] + p[1] + p[2];
    return !(s > s - 1.0);
  }
  static CVec3 Axis(unsigned int idim) {
    CVec3 r(0, 0, 0);
    if (idim < 3) { r[idim] = 1; }
    return r;
  }
  
  using Scalar = T;

 public:
  union {
    T p[3];
    struct {
      T x, y, z;
    };
  };
};
using CVec3d = CVec3<double>;
using CVec3f = CVec3<float>;

// --------------------------------------------------------------------------
// rule about naming, the method starts "Set" change it self (not const)

template<typename T>
CVec3<T> Mat3Vec(
    const T M[9],
    const CVec3<T> &v);

template<typename T>
CVec3<T> Mat4Vec(
    const T M[16],
    const CVec3<T> &v);

template<typename T>
DFM2_INLINE CVec3<T> QuatVec(
    const T quat[4],
    const CVec3<T> &v0);

template<typename REAL>
CVec3<REAL> QuatConjVec(
    const REAL quat[4],
    const CVec3<REAL> &v0);

template<typename T>
bool operator==(const CVec3<T> &lhs, const CVec3<T> &rhs);

template<typename T>
bool operator!=(const CVec3<T> &lhs, const CVec3<T> &rhs);

// ---------------------------------------------

template<typename T>
CVec3<T> RotateVector(
    const CVec3<T> &vec0,
    const CVec3<T> &rot);

/**
 * @brief 3x3 Rotation matrix to rotate V into v with minimum rotation angle
 * @param[in] V rotation from
 * @param[in] v rotation to
 */
// TODO: consider making this function general to Eigen::Vector3x
//  using template and moving this functin to "geo_vec3.h"
template<typename REAL>
DFM2_INLINE std::array<REAL, 9> Mat3_MinimumRotation(
    const CVec3<REAL> &V,
    const CVec3<REAL> &v);

// TODO: consider making this function general to Eigen::Vector3x
//  using template and moving this functin to "geo_vec3.h"
template<typename REAL>
DFM2_INLINE std::array<REAL, 9> Mat3_ParallelTransport(
    const CVec3<REAL> &p0,
    const CVec3<REAL> &p1,
    const CVec3<REAL> &q0,
    const CVec3<REAL> &q1);

} // namespace delfem2

#if defined(_MSC_VER)
#  pragma warning( pop )
#endif

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/vec3.cpp"
#endif

#endif // DFM2_VEC3_H
