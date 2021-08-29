/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/


#ifndef DFM2_VEC3_H
#define DFM2_VEC3_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delfem2/dfm2_inline.h"

#if defined(_MSC_VER)
#  pragma warning( push )
#  pragma warning( disable : 4201 )
#endif

#define NEARLY_ZERO 1.e-16

namespace delfem2 {

template<typename T>
DFM2_INLINE T Distance3(const T p0[3], const T p1[3]);

template<typename T>
DFM2_INLINE T SquareDistance3(const T p0[3], const T p1[3]);

template<typename T>
DFM2_INLINE T Length3(const T v[3]);

template<typename T>
DFM2_INLINE T SquareLength3(const T v[3]);

template<typename T>
DFM2_INLINE T Dot3(const T a[3], const T b[3]);

template<typename T>
DFM2_INLINE void Normalize3(T v[3]);

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

// above: functions general for any dimensions
// ------------------------------------------------
// below: functions specific to 3 dimension

template<typename T0, typename T1, typename T2>
void Cross3(
    T0 r[3],
    const T1 v1[3],
    const T2 v2[3]);

template<typename T>
T Area_Tri3(const T v1[3], const T v2[3], const T v3[3]);

template<typename T>
T ScalarTripleProduct3(
    const T a[], const T b[], const T c[]);

template<typename T>
void NormalTri3(T n[3],
                const T v1[3], const T v2[3], const T v3[3]);

template<typename REAL>
void UnitNormalAreaTri3(
    REAL n[3], REAL &a,
    const REAL v1[3], const REAL v2[3], const REAL v3[3]);

// --------------------------------

DFM2_INLINE void GetVertical2Vector3D(
    const double vec_n[3],
    double vec_x[3],
    double vec_y[3]);

// -------------------------------------------------------------

template<typename T>
class CVec3;

template<typename T0, typename T1>
CVec3<T0> operator*(T1 d, const CVec3<T0> &rhs);

template<typename T>
CVec3<T> operator/(const CVec3<T> &vec, T d);

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
  CVec3(const CVec3 &rhs) {
    p[0] = rhs.p[0];
    p[1] = rhs.p[1];
    p[2] = rhs.p[2];
  }

  template<typename S>
  explicit CVec3(const S *v) : p{v[0], v[1], v[2]} {}

  template<typename S>
  explicit CVec3(const std::vector<S> &v)  : p(v[0], v[1], v[2]) {}

  virtual ~CVec3() = default;

  // above: constructor / destructor
  // -------------------------------
  // below: operator

  inline CVec3 operator-() const { return ((T) (-1)) * (*this); }
  inline CVec3 operator+() const { return *this; }
  inline CVec3 operator+() { return *this; }
  inline CVec3 operator-() { return CVec3(-p[0], -p[1], -p[2]); }
  inline CVec3 &operator=(const CVec3 &b) {
    if (this != &b) {
      x = b.x;
      y = b.y;
      z = b.z;
    }
    return *this;
  }
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

  inline CVec3 operator^(const CVec3 &b) const {
    return CVec3(y * b.z - z * b.y,
                 z * b.x - x * b.z,
                 x * b.y - y * b.x);
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
  //
  template<typename DIST, typename ENG>
  void SetRandom(DIST &dist, ENG &eng) {
    p[0] = dist(eng);
    p[1] = dist(eng);
    p[2] = dist(eng);
  }
  //
  template<typename DIST, typename ENG>
  static CVec3 Random(DIST &dist, ENG &eng) {
    return CVec3(dist(eng), dist(eng), dist(eng));
  }
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

// ------------------------------

template<typename T>
T Dot(const CVec3<T> &arg1, const CVec3<T> &arg2);

template<typename T>
CVec3<T> Cross(const CVec3<T> &arg1, const CVec3<T> &arg2);


// --------------------------------------------------------------------------
// rule about naming, the method starts "Set" change it self (not const)

template<typename T>
CVec3<T> mult_GlAffineMatrix(
    const float *m,
    const CVec3<T> &p);

template<typename T>
CVec3<T> solve_GlAffineMatrix(
    const float *m,
    const CVec3<T> &p);

template<typename T>
CVec3<T> solve_GlAffineMatrixDirection(
    const float *m,
    const CVec3<T> &v);

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
CVec3<T> screenProjection(
    const CVec3<T> &v,
    const float *mMV,
    const float *mPj);

// opposite to the screen normal
template<typename T>
CVec3<T> screenDepthDirection(
    const CVec3<T> &v,
    const float *mMV,
    const float *mPj);

template<typename T>
CVec3<T> screenUnProjection(
    const CVec3<T> &v,
    const float *mMV,
    const float *mPj);

template<typename T>
CVec3<T> screenUnProjectionDirection(
    const CVec3<T> &v,
    const float *mMV,
    const float *mPj);

template<typename T>
T ScalarTripleProduct(
    const CVec3<T> &a, const CVec3<T> &b, const CVec3<T> &c);

template<typename T>
bool operator==(const CVec3<T> &lhs, const CVec3<T> &rhs);

template<typename T>
bool operator!=(const CVec3<T> &lhs, const CVec3<T> &rhs);

template<typename T>
void GetVertical2Vector(
    const CVec3<T> &vec_n, CVec3<T> &vec_x, CVec3<T> &vec_y);

// ---------------------------------------------

template<typename T>
void Cross(CVec3<T> &lhs, const CVec3<T> &v1, const CVec3<T> &v2);

template<typename T>
T Area_Tri(const CVec3<T> &v1, const CVec3<T> &v2, const CVec3<T> &v3);

template<typename T>
T SquareTriArea(const CVec3<T> &v1, const CVec3<T> &v2, const CVec3<T> &v3);

template<typename T>
double SquareDistance(
    const CVec3<T> &ipo0,
    const CVec3<T> &ipo1);

template<typename T>
double SquareLength(const CVec3<T> &point);

template<typename T>
double Length(
    const CVec3<T> &point);

template<typename T>
double Distance(
    const CVec3<T> &p0,
    const CVec3<T> &p1);

template<typename T>
double SqareLongestEdgeLength(
    const CVec3<T> &ipo0,
    const CVec3<T> &ipo1,
    const CVec3<T> &ipo2,
    const CVec3<T> &ipo3);

template<typename T>
double LongestEdgeLength(
    const CVec3<T> &ipo0,
    const CVec3<T> &ipo1,
    const CVec3<T> &ipo2,
    const CVec3<T> &ipo3);

template<typename T>
double SqareShortestEdgeLength(
    const CVec3<T> &ipo0,
    const CVec3<T> &ipo1,
    const CVec3<T> &ipo2,
    const CVec3<T> &ipo3);

template<typename T>
double ShortestEdgeLength(
    const CVec3<T> &ipo0,
    const CVec3<T> &ipo1,
    const CVec3<T> &ipo2,
    const CVec3<T> &ipo3);

template<typename T>
void Normal(
    CVec3<T> &vnorm,
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3);

template<typename T>
CVec3<T> Normal(
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3);

template<typename T>
void UnitNormal(
    CVec3<T> &vnorm,
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3);

template<typename T>
CVec3<T> UnitNormal(
    const CVec3<T> &v1,
    const CVec3<T> &v2,
    const CVec3<T> &v3);

template<typename T>
CVec3<T> RotateVector(const CVec3<T> &vec0, const CVec3<T> &rot);

template<typename T>
CVec3<T> RandVector();

template<typename T>
CVec3<T> RandUnitVector();

template<typename T>
CVec3<T> RandGaussVector();

// ----------------------------------------------------------
// here starts std::vector<CVector3>

template<typename T>
double Volume_Tet(
    int iv1,
    int iv2,
    int iv3,
    int iv4,
    const std::vector<CVec3<T> > &aPoint);

template<typename T>
double Area_Tri(
    int iv1,
    int iv2,
    int iv3,
    const std::vector<CVec3<T> > &aPoint);

template<typename T>
CVec3<T> CG_Tri3(
    unsigned int itri,
    const std::vector<unsigned int> &aTri,
    const std::vector<T> &aXYZ) {
  const unsigned int i0 = aTri[itri * 3 + 0];
  const unsigned int i1 = aTri[itri * 3 + 1];
  const unsigned int i2 = aTri[itri * 3 + 2];
  CVec3<T> p;
  p.p[0] = (aXYZ[i0 * 3 + 0] + aXYZ[i1 * 3 + 0] + aXYZ[i2 * 3 + 0]) / 3.0;
  p.p[1] = (aXYZ[i0 * 3 + 1] + aXYZ[i1 * 3 + 1] + aXYZ[i2 * 3 + 1]) / 3.0;
  p.p[2] = (aXYZ[i0 * 3 + 2] + aXYZ[i1 * 3 + 2] + aXYZ[i2 * 3 + 2]) / 3.0;
  return p;
}

template<typename T>
inline CVec3<T> Normal_Tri3(
    unsigned int itri,
    const std::vector<unsigned int> &aTri,
    const std::vector<T> &aXYZ) {
  const unsigned int i0 = aTri[itri * 3 + 0];
  const unsigned int i1 = aTri[itri * 3 + 1];
  const unsigned int i2 = aTri[itri * 3 + 2];
  const CVec3<T> p0(aXYZ[i0 * 3 + 0], aXYZ[i0 * 3 + 1], aXYZ[i0 * 3 + 2]);
  const CVec3<T> p1(aXYZ[i1 * 3 + 0], aXYZ[i1 * 3 + 1], aXYZ[i1 * 3 + 2]);
  const CVec3<T> p2(aXYZ[i2 * 3 + 0], aXYZ[i2 * 3 + 1], aXYZ[i2 * 3 + 2]);
  return (p1 - p0) ^ (p2 - p0);
}

} // namespace delfem2

#if defined(_MSC_VER)
#  pragma warning( pop )
#endif

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/vec3.cpp"
#endif

#endif // DFM2_VEC3_H
