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
#include "delfem2/mat3.h"
#include "delfem2/mat4.h"
#include "delfem2/quat.h"


// ----------------------------------------

DFM2_INLINE delfem2::CVec2d delfem2::screenXYProjection(
    const CVec3d &v,
    const float *mMV,
    const float *mPj) {
  CVec3d sp0 = screenProjection(v, mMV, mPj);
  return delfem2::CVec2d(sp0.x, sp0.y);
}

DFM2_INLINE delfem2::CVec3d delfem2::GetCartesianRotationVector(
    const CMat3d &m) {
  const double *mat = m.p_;
  CVec3d a;
  a.p[0] = mat[7] - mat[5];
  a.p[1] = mat[2] - mat[6];
  a.p[2] = mat[3] - mat[1];
  double act = (m.Trace() - 1) * 0.5;
  if (act > +1) { act = +1; }
  if (act < -1) { act = -1; }
  double theta = acos(act);
  if (myIsNAN_Matrix3(theta)) { return a; }
  if (fabs(theta) < 1.0e-5) { return a * 0.5; }
  double mag = 0.5 * theta / sin(theta);
  a *= mag;
  return a;
}

DFM2_INLINE delfem2::CVec3d delfem2::GetSpinVector(
    const CMat3d &m) {
  const double *mat = m.p_;
  CVec3d r;
  r.p[0] = (mat[7] - mat[5]) * 0.5;
  r.p[1] = (mat[2] - mat[6]) * 0.5;
  r.p[2] = (mat[3] - mat[1]) * 0.5;
  return r;
}

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

DFM2_INLINE void delfem2::SetRotMatrix_Cartesian(
    CMat3d &m,
    const CVec3d &v) {
  m.SetRotMatrix_Cartesian(v.p);
}

DFM2_INLINE void delfem2::SetSpinTensor(CMat3d &m, const CVec3d &vec0) {
  Mat3_Spin(m.p_, vec0.p);
}

DFM2_INLINE void delfem2::SetOuterProduct(
    CMat3d &m,
    const CVec3d &vec0,
    const CVec3d &vec1) {
  double *mat = m.p_;
  mat[0] = vec0.x * vec1.x;
  mat[1] = vec0.x * vec1.y;
  mat[2] = vec0.x * vec1.z;
  mat[3] = vec0.y * vec1.x;
  mat[4] = vec0.y * vec1.y;
  mat[5] = vec0.y * vec1.z;
  mat[6] = vec0.z * vec1.x;
  mat[7] = vec0.z * vec1.y;
  mat[8] = vec0.z * vec1.z;
}

DFM2_INLINE void delfem2::SetProjection(CMat3d &m, const CVec3d &vec0) {
  double *mat = m.p_;
  const CVec3d &u = vec0.normalized();
  mat[0] = 1 - u.x * u.x;
  mat[1] = 0 - u.x * u.y;
  mat[2] = 0 - u.x * u.z;
  mat[3] = 0 - u.y * u.x;
  mat[4] = 1 - u.y * u.y;
  mat[5] = 0 - u.y * u.z;
  mat[6] = 0 - u.z * u.x;
  mat[7] = 0 - u.z * u.y;
  mat[8] = 1 - u.z * u.z;
}

// ----------------------------

DFM2_INLINE delfem2::CMat3d delfem2::Mirror(const CVec3d &n) {
  CVec3d N = n;
  N.normalize();
  return CMat3d::Identity() - 2 * delfem2::Mat3_OuterProduct(N, N);
}

DFM2_INLINE delfem2::CMat3d delfem2::Mat3_CrossCross(const CVec3d &v) {
  return Mat3(v) * Mat3(v);
}

DFM2_INLINE delfem2::CMat3d delfem2::RotMatrix_Cartesian(const CVec3d &v) {
  CMat3d m;
  SetRotMatrix_Cartesian(m, v);
  return m;
}

DFM2_INLINE delfem2::CMat3d delfem2::Mat3(const CVec3d &vec0) {
  CMat3d m;
  SetSpinTensor(m, vec0);
  return m;
}

DFM2_INLINE delfem2::CMat3d delfem2::Mat3(const CVec3d &vec0, const CVec3d &vec1) {
  CMat3d m;
  SetOuterProduct(m, vec0, vec1);
  return m;
}

DFM2_INLINE delfem2::CMat3d delfem2::Mat3(
    const CVec3d &vec0,
    const CVec3d &vec1,
    const CVec3d &vec2) {
  CMat3d m;
  double *mat = m.p_;
  mat[0 * 3 + 0] = vec0.x;
  mat[0 * 3 + 1] = vec1.x;
  mat[0 * 3 + 2] = vec2.x;
  mat[1 * 3 + 0] = vec0.y;
  mat[1 * 3 + 1] = vec1.y;
  mat[1 * 3 + 2] = vec2.y;
  mat[2 * 3 + 0] = vec0.z;
  mat[2 * 3 + 1] = vec1.z;
  mat[2 * 3 + 2] = vec2.z;
  return m;
}

DFM2_INLINE delfem2::CMat3d delfem2::Mat3_Spin(const CVec3d &vec0) {
  CMat3d m;
  ::delfem2::Mat3_Spin(m.p_, vec0.p);
  return m;
}

DFM2_INLINE delfem2::CMat3d delfem2::Mat3_OuterProduct(const CVec3d &vec0, const CVec3d &vec1) {
  CMat3d m;
  SetOuterProduct(m, vec0, vec1);
  return m;
}

DFM2_INLINE delfem2::CMat3d delfem2::Mat3_RotCartesian(const CVec3d &vec0) {
  CMat3d m;
  m.SetRotMatrix_Cartesian(vec0.x, vec0.y, vec0.z);
  return m;
}

// ------------------

namespace delfem2 {

template <typename T>
DFM2_INLINE CVec3<T> operator*(
    const CVec3<T> &v,
    const CMat3<T> &m) {
  return MatVecTrans(m, v);
}
#ifdef DFM2_STATIC_LIBRARY
template CVec3d operator*(const CVec3d &, const CMat3d &);
template CVec3f operator*(const CVec3f &, const CMat3f &);
#endif

// -------------------------------------------

template <typename T>
DFM2_INLINE CVec3<T> operator*(
    const CMat3<T> &m,
    const CVec3<T> &v) {
  return MatVec(m, v);
}
#ifdef DFM2_STATIC_LIBRARY
template CVec3d operator*(const CMat3d &, const CVec3d &);
template CVec3f operator*(const CMat3f &, const CVec3f &);
#endif

}

// ------------------------------

template<typename REAL>
delfem2::CMat3<REAL> delfem2::Mat3_MinimumRotation(
    const CVec3<REAL> &V,
    const CVec3<REAL> &v) {
  CVec3<REAL> ep = V.normalized();
  CVec3<REAL> eq = v.normalized();
  CVec3<REAL> n = ep ^ eq;
  const double st2 = n.dot(n);
  CMat3<REAL> m;
  if (st2 < 1.0e-4f) {
    m.p_[0] = 1.f + 0.5f * (n.x * n.x - st2);
    m.p_[1] = -n.z + 0.5f * (n.x * n.y);
    m.p_[2] = +n.y + 0.5f * (n.x * n.z);
    m.p_[3] = +n.z + 0.5f * (n.y * n.x);
    m.p_[4] = 1.f + 0.5f * (n.y * n.y - st2);
    m.p_[5] = -n.x + 0.5f * (n.y * n.z);
    m.p_[6] = -n.y + 0.5f * (n.z * n.x);
    m.p_[7] = +n.x + 0.5f * (n.z * n.y);
    m.p_[8] = 1.f + 0.5f * (n.z * n.z - st2);
    return m;
  }
  const double st = sqrt(st2);
  const double ct = ep.dot(eq);
  n.normalize();
  m.p_[0] = ct + (1.f - ct) * n.x * n.x;
  m.p_[1] = -n.z * st + (1.f - ct) * n.x * n.y;
  m.p_[2] = +n.y * st + (1.f - ct) * n.x * n.z;
  m.p_[3] = +n.z * st + (1.f - ct) * n.y * n.x;
  m.p_[4] = ct + (1.f - ct) * n.y * n.y;
  m.p_[5] = -n.x * st + (1.f - ct) * n.y * n.z;
  m.p_[6] = -n.y * st + (1.f - ct) * n.z * n.x;
  m.p_[7] = +n.x * st + (1.f - ct) * n.z * n.y;
  m.p_[8] = ct + (1.f - ct) * n.z * n.z;
  return m;
}
#ifdef DFM2_STATIC_LIBRARY
template delfem2::CMat3d delfem2::Mat3_MinimumRotation(const CVec3d& V, const CVec3d& v);
#endif


// --------------------------

DFM2_INLINE delfem2::CMat3d delfem2::Mat3_ParallelTransport(
    const CVec3d &p0,
    const CVec3d &p1,
    const CVec3d &q0,
    const CVec3d &q1) {
  return Mat3_MinimumRotation(p1 - p0, q1 - q0);
}

// -----------------------------------------------------
// below: rotational inertia

// moment of inertia around origin triangle vtx (origin,d0,d1,d2) the area_density=1
// see http://www.dcs.warwick.ac.uk/~rahil/files/RigidBodySimulation.pdf
DFM2_INLINE delfem2::CMat3d delfem2::Mat3_IrotTri(
    const CVec3d &d0,
    const CVec3d &d1,
    const CVec3d &d2) {

  CVec3d dv = d0 + d1 + d2;
  CMat3d I0 =
      Mat3_OuterProduct(d0, d0) +
      Mat3_OuterProduct(d1, d1) +
      Mat3_OuterProduct(d2, d2) +
      Mat3_OuterProduct(dv, dv);
  double tr0 = I0.Trace();
  CMat3d I = tr0 * CMat3d::Identity() - I0;

  double darea = ((d1 - d0) ^ (d2 - d0)).norm();
  I *= darea / 24.0;
  return I;
}

// moment of inertia triangle pyramid with vtx (origin,d0,d1,d2) volume_density = 1
// see http://www.dcs.warwick.ac.uk/~rahil/files/RigidBodySimulation.pdf
DFM2_INLINE delfem2::CMat3d delfem2::Mat3_IrotTriSolid(
    const CVec3d &d0,
    const CVec3d &d1,
    const CVec3d &d2) {
  CVec3d dv = d0 + d1 + d2;
  CMat3d I0 =
      Mat3_OuterProduct(d0, d0) +
      Mat3_OuterProduct(d1, d1) +
      Mat3_OuterProduct(d2, d2) +
      Mat3_OuterProduct(dv, dv);
  double tr0 = I0.Trace();
  CMat3d I = tr0 * CMat3d::Identity() - I0;

  double darea = d0.dot(d1 ^ d2);
  I *= darea / 120.0;
  return I;
}

DFM2_INLINE delfem2::CMat3d delfem2::Mat3_IrotLineSeg(
    const CVec3d &d0,
    const CVec3d &d1) {
  CVec3d dv = d1 - d0;
  double l = dv.norm();
  CMat3d I;
  {
    I = dv.squaredNorm() * CMat3d::Identity() - Mat3_OuterProduct(dv, dv);
    I *= l / 12.0;
  }
  CVec3d p = (d0 + d1) * 0.5;
  I += l * (p.squaredNorm() * CMat3d::Identity() - Mat3_OuterProduct(p, p));
  return I;
}

DFM2_INLINE delfem2::CMat3d delfem2::Mat3_IrotPoint(
    const CVec3d &d0) {
  return (d0.squaredNorm() * CMat3d::Identity() - Mat3_OuterProduct(d0, d0));
}


// above: rotational inertia
// ---------------------------------------------------------------------


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

// ---------------------------------------------------------------------------------------


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

DFM2_INLINE void delfem2::UpdateRotationsByMatchingCluster_Linear(
    std::vector<double> &aQuat1,
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aXYZ1,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup) {
  const size_t np = aXYZ0.size() / 3;
  for (unsigned int ip = 0; ip < np; ++ip) {
    const CVec3d Pi(aXYZ0.data() + ip * 3);
    const CVec3d pi(aXYZ1.data() + ip * 3);
    const CQuatd Qi(aQuat1.data() + ip * 4);
    CMat3d Mat;
    Mat.setZero();
    CVec3d rhs;
    rhs.setZero();
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      const unsigned int jp = psup[ipsup];
      const CVec3d v0 = Qi * (CVec3d(aXYZ0.data() + jp * 3) - Pi);
      const CVec3d d01 = CVec3d(aXYZ1.data() + jp * 3) - pi - v0;
      Mat += Mat3_CrossCross(v0);
      rhs += d01 ^ v0;
    }
    CVec3d sol = Mat.Inverse() * rhs;
    CQuatd q0 = Quat_CartesianAngle(sol);
    CQuatd q1 = q0 * Qi;
    q1.CopyTo(aQuat1.data() + ip * 4);
  }
}

DFM2_INLINE void delfem2::UpdateRotationsByMatchingCluster_SVD(
    std::vector<double> &aQuat1,
    unsigned int ip,
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aXYZ1,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup) {
  const CVec3d Pi(aXYZ0.data() + ip * 3);
  const CVec3d pi(aXYZ1.data() + ip * 3);
  const CMat3d R0i = CMat3d::Quat(aQuat1.data() + ip * 4);
  double A[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
    const unsigned int jp = psup[ipsup];
    const CVec3d &v0 = R0i * (CVec3d(aXYZ0.data() + jp * 3) - Pi);
    const CVec3d v1 = CVec3d(aXYZ1.data() + jp * 3) - pi;
    const double *dp = v1.p;
    const double *dq = v0.p;
    A[0 * 3 + 0] += dp[0] * dq[0];
    A[0 * 3 + 1] += dp[0] * dq[1];
    A[0 * 3 + 2] += dp[0] * dq[2];
    A[1 * 3 + 0] += dp[1] * dq[0];
    A[1 * 3 + 1] += dp[1] * dq[1];
    A[1 * 3 + 2] += dp[1] * dq[2];
    A[2 * 3 + 0] += dp[2] * dq[0];
    A[2 * 3 + 1] += dp[2] * dq[1];
    A[2 * 3 + 2] += dp[2] * dq[2];
  }
  CMat3d dRi;
  GetRotPolarDecomp(dRi.p_, A, 40);
  CMat3d R1 = dRi * R0i;
  CQuatd q1;
  R1.GetQuat_RotMatrix(q1.p);
  q1.CopyTo(aQuat1.data() + ip * 4);
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
