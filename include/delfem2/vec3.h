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
#  pragma warning( disable : 4201 )  // because we use nameless union in CVec3
#endif

namespace delfem2 {

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
  CVec3(const std::array<S, 3> &v) : p{
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
  inline CVec3 operator-() { return {-x, -y, -z}; }
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
    return {x + b.x, y + b.y, z + b.z};
  }

  inline CVec3 operator-(const CVec3 &b) const {
    return {x - b.x, y - b.y, z - b.z};
  }

  inline CVec3 operator*(T b) const {
    return {x * b, y * b, z * b};
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

  template<typename T1, typename std::enable_if_t<std::is_scalar_v<T1>> * = nullptr>
  CVec3 operator/(T1 d) const {
    if (fabs(d) < NEARLY_ZERO) { return CVec3(0, 0, 0); }
    return {x / d, y / d, z / d};
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
  // implicit cast
  operator std::array<T, 3>() const {
    return { p[0], p[1], p[2] };
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
  void normalize() {
    T invmag = 1 / norm();
    p[0] *= invmag;
    p[1] *= invmag;
    p[2] *= invmag;
  }

  //! @details named after Eigen library
  inline T norm() const {
    return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
  }

  //! @details named after Eigen library
  inline T squaredNorm() const {
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
  }

  //! @details named after Eigen library
  void setZero() {
    p[0] = 0;
    p[1] = 0;
    p[2] = 0;
  }

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
  using value_type = Scalar;
  static constexpr T NEARLY_ZERO = static_cast<T>(1.e-16);

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

template<typename T0, typename T1, typename std::enable_if_t<std::is_scalar_v<T1>> * = nullptr>
CVec3<T0> operator*(T1 d, const CVec3<T0> &rhs) {
  return {
      static_cast<T0>(rhs.x * d),
      static_cast<T0>(rhs.y * d),
      static_cast<T0>(rhs.z * d)};
}

template<typename T>
std::ostream &operator<<(
    std::ostream &output,
    const CVec3<T> &v) {
  output.setf(std::ios::scientific);
  output << v.x << " " << v.y << " " << v.z;
  return output;
}

template<typename T>
std::istream &operator>>(
    std::istream &input,
    CVec3<T> &v) {
  input >> v.x >> v.y >> v.z;
  return input;
}

template<typename T>
bool operator==(const CVec3<T> &lhs, const CVec3<T> &rhs) {
  if (fabs(lhs.p[0] - rhs.p[0]) < CVec3<T>::NEARLY_ZERO
      && fabs(lhs.p[1] - rhs.p[1]) < CVec3<T>::NEARLY_ZERO
      && fabs(lhs.p[2] - rhs.p[2]) < CVec3<T>::NEARLY_ZERO) { return true; }
  else { return false; }
}

template<typename T>
bool operator!=(const CVec3<T> &lhs, const CVec3<T> &rhs) {
  return !(lhs == rhs);
}

} // namespace delfem2

#if defined(_MSC_VER)
#  pragma warning( pop )
#endif

#endif // DFM2_VEC3_H
