/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/geo3_v23m34q.h"

#include <cmath>

#include "delfem2/vec2.h"
#include "delfem2/vec3.h"
#include "delfem2/vec3_funcs.h"
#include "delfem2/mat3.h"
#include "delfem2/mat3_funcs.h"
#include "delfem2/mat4.h"
#include "delfem2/quat.h"
#include "delfem2/geo_tri.h"

// ----------------------------------------

template<typename T>
DFM2_INLINE delfem2::CVec3<T>
delfem2::MatVec(
    const CMat3<T> &m,
    const CVec3<T> &vec0) {
  CVec3<T> vec1;
  delfem2::MatVec3(vec1.p, m.p_, vec0.p);
  return vec1;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::MatVec(const CMat3d& m, const CVec3d& vec0);
template delfem2::CVec3f delfem2::MatVec(const CMat3f& m, const CVec3f& vec0);
#endif

template<typename T>
DFM2_INLINE delfem2::CVec3<T> delfem2::MatVecTrans(
    const CMat3<T> &m,
    const CVec3<T> &vec0) {
  CVec3<T> vec1;
  MatTVec3(vec1.p, m.p_, vec0.p);
  return vec1;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CVec3d delfem2::MatVecTrans(const CMat3d&, const CVec3d&);
template delfem2::CVec3f delfem2::MatVecTrans(const CMat3f&, const CVec3f&);
#endif

// ---------------------------------------------------------------------

DFM2_INLINE void delfem2::SetDiag(
    CMat3d &m,
    const CVec3d &d) {
  double *mat = m.p_;
  mat[0 * 3 + 0] = d.x;
  mat[1 * 3 + 1] = d.y;
  mat[2 * 3 + 2] = d.z;
}


// -----------------------------------------------------
// below: rotational inertia

DFM2_INLINE void delfem2::Mat4_MatTransl(
    double m[16],
    const CMat3d &mat,
    const CVec3d &trans) {
  mat.AffineMatrixTrans(m);
  m[3 * 4 + 0] = trans.x;
  m[3 * 4 + 1] = trans.y;
  m[3 * 4 + 2] = trans.z;
}

DFM2_INLINE void delfem2::Mat4_ScaleMatTransl(
    double m[16],
    double scale,
    const CMat3d &mat,
    const CVec3d &trans) {
  mat.AffineMatrixTrans(m);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      m[i * 4 + j] *= scale;
    }
  }
  m[3 * 4 + 0] = trans.x;
  m[3 * 4 + 1] = trans.y;
  m[3 * 4 + 2] = trans.z;
}

// ----------------------------------------------------
// quaternion

namespace delfem2 {

template<typename REAL>
DFM2_INLINE CVec3<REAL> operator*(
    const CQuat<REAL> &q,
    const CVec3<REAL> &v) {
  CVec3<REAL> p;
  QuatVec(p.p,
          q.p, v.p);
  return p;
}
#ifdef DFM2_STATIC_LIBRARY
template CVec3f operator* (const CQuatf& q, const CVec3f& v);
template CVec3d operator* (const CQuatd& q, const CVec3d& v);
#endif

}

// ------------

DFM2_INLINE delfem2::CQuatd delfem2::Quat_CartesianAngle(const CVec3d &p) {
  CQuatd q;
  Quat_CartesianAngle(q.p, p.p);
  return q;
}

bool delfem2::Distortion_MappingTriangleFrom2To3Dim(
    double thresA,
    double thresE,
    unsigned int it0,
    const std::vector<unsigned int> &aTri,
    const std::vector<double> &aXYZ,
    const std::vector<double> &aTexP) { // check the distortion
  const unsigned int i0 = aTri[it0 * 3 + 0];
  const unsigned int i1 = aTri[it0 * 3 + 1];
  const unsigned int i2 = aTri[it0 * 3 + 2];
  const double area2 = Area_Tri2(aTexP.data() + i0 * 2, aTexP.data() + i1 * 2, aTexP.data() + i2 * 2);
  const double area3 = Area_Tri3(aXYZ.data() + i0 * 3, aXYZ.data() + i1 * 3, aXYZ.data() + i2 * 3);
  const double scoreArea = 0.5 * (area2 / area3 + area3 / area2);
  if (scoreArea < 0 || scoreArea > thresA) { return true; }
  const double len12 = Distance2(aTexP.data() + i1 * 2, aTexP.data() + i2 * 2);
  const double len20 = Distance2(aTexP.data() + i2 * 2, aTexP.data() + i0 * 2);
  const double len01 = Distance2(aTexP.data() + i0 * 2, aTexP.data() + i1 * 2);
  const double Len12 = Distance3(aXYZ.data() + i1 * 3, aXYZ.data() + i2 * 3);
  const double Len20 = Distance3(aXYZ.data() + i2 * 3, aXYZ.data() + i0 * 3);
  const double Len01 = Distance3(aXYZ.data() + i0 * 3, aXYZ.data() + i1 * 3);
  if (0.5 * (len12 / Len12 + Len12 / len12) > thresE) { return true; }
  if (0.5 * (len20 / Len20 + Len20 / len20) > thresE) { return true; }
  if (0.5 * (len01 / Len01 + Len01 / len01) > thresE) { return true; }
  return false;
}
