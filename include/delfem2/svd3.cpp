//
// Created by Nobuyuki Umetani on 2022/01/30.
//

#include "svd3.h"

#include <cmath>

#include "delfem2/dfm2_inline.h"

// ---------------------
// local functions

namespace delfem2::svd3 {

template<typename T>
DFM2_INLINE T SquareNormFrobenius_SymMat3(
    const T sm[6]) {
  return sm[0] * sm[0] + sm[1] * sm[1] + sm[2] * sm[2] + 2 * (sm[3] * sm[3] + sm[4] * sm[4] + sm[5] * sm[5]);
}

template<typename T>
DFM2_INLINE T Det_Mat3(const T U[9]) {
  return +U[0] * U[4] * U[8] + U[3] * U[7] * U[2] + U[6] * U[1] * U[5]
      - U[0] * U[7] * U[5] - U[6] * U[4] * U[2] - U[3] * U[1] * U[8];
}

DFM2_INLINE void Cross3(double r[3], const double v1[3], const double v2[3]) {
  r[0] = v1[1] * v2[2] - v2[1] * v1[2];
  r[1] = v1[2] * v2[0] - v2[2] * v1[0];
  r[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

DFM2_INLINE void Normalize3(double v[3]) {
  double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

DFM2_INLINE double SqLength3(const double v[3]) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

template<typename T>
DFM2_INLINE void MatMatT3(
    T *ABt,
    const T *A,
    const T *B) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ABt[i * 3 + j] =
          A[i * 3 + 0] * B[j * 3 + 0] +
              A[i * 3 + 1] * B[j * 3 + 1] +
              A[i * 3 + 2] * B[j * 3 + 2];
    }
  }
}

DFM2_INLINE void SortEigen3(
    double G[3],
    double V[9]) {
  double t;
  if (G[1] > G[0]) {  // swap 01
    t = G[0];
    G[0] = G[1];
    G[1] = t;
    t = V[0];
    V[0] = V[1];
    V[1] = t;
    t = V[3];
    V[3] = V[4];
    V[4] = t;
    t = V[6];
    V[6] = V[7];
    V[7] = t;
  }
  if (G[2] > G[1]) {
    t = G[1];
    G[1] = G[2];
    G[2] = t;
    t = V[1];
    V[1] = V[2];
    V[2] = t;
    t = V[4];
    V[4] = V[5];
    V[5] = t;
    t = V[7];
    V[7] = V[8];
    V[8] = t;
  }
  if (G[1] > G[0]) { // swap 01
    t = G[0];
    G[0] = G[1];
    G[1] = t;
    t = V[0];
    V[0] = V[1];
    V[1] = t;
    t = V[3];
    V[3] = V[4];
    V[4] = t;
    t = V[6];
    V[6] = V[7];
    V[7] = t;
  }
}

}  // namespace delfem2::svd3


// -------------------
// implementation

DFM2_INLINE bool delfem2::EigenSym3(
    double u[9],
    double l[3],
    const double sm[6],
    int nitr) {
  // initialize u as identity matrix
  u[0] = u[4] = u[8] = 1.0;
  u[1] = u[2] = u[3] = u[5] = u[6] = u[7] = 0.0;
  l[0] = l[1] = l[2] = 0.0;
  double dnrm = svd3::SquareNormFrobenius_SymMat3(sm);
  if (dnrm < 1.0e-30) { return false; } // this matrix is too small
  const double scale = sqrt(dnrm);
  const double invscl = 1.0 / scale;
  double sms[6] = {sm[0] * invscl, sm[1] * invscl, sm[2] * invscl, sm[3] * invscl, sm[4] * invscl, sm[5] * invscl};
  for (int itr = 0; itr < nitr; itr++) {
    const double m[6] = {sms[0], sms[1], sms[2], sms[3], sms[4], sms[5]};
    const double v[9] = {u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8]};
    const double a12 = fabs(sms[3]);
    const double a20 = fabs(sms[4]);
    const double a01 = fabs(sms[5]);
    if (a12 >= a20 && a12 >= a01) {
      // when a12 sms[3] is the biggest
      const double t = 0.5 * atan2(2 * m[3], m[2] - m[1]);
      const double ct = cos(t);
      const double st = sin(t);
      sms[1] = ct * ct * m[1] + st * st * m[2] - 2 * st * ct * m[3];
      sms[2] = ct * ct * m[2] + st * st * m[1] + 2 * st * ct * m[3];
      sms[3] = 0; // (ct*ct-st*st)*m[3]+st*ct*(m[1]-m[2]);
      sms[4] = st * m[5] + ct * m[4];
      sms[5] = ct * m[5] - st * m[4];
      //
      u[1] = +ct * v[1] - st * v[2];
      u[2] = +st * v[1] + ct * v[2];
      u[4] = +ct * v[4] - st * v[5];
      u[5] = +st * v[4] + ct * v[5];
      u[7] = +ct * v[7] - st * v[8];
      u[8] = +st * v[7] + ct * v[8];
    } else if (a20 >= a01 && a20 >= a12) {
      // when a20 sms[4] is the biggest
      // the above condition statement shoud pass exactly once for each iteration.
      const double t = 0.5 * atan2(2 * m[4], m[2] - m[0]);
      const double ct = cos(t);
      const double st = sin(t);
      sms[0] = ct * ct * m[0] + st * st * m[2] - 2 * st * ct * m[4];
      sms[2] = ct * ct * m[2] + st * st * m[0] + 2 * st * ct * m[4];
      sms[3] = st * m[5] + ct * m[3];
      sms[4] = 0; // (ct*ct-st*st)*m[4]+st*ct*(m[0]-m[2]);
      sms[5] = ct * m[5] - st * m[3];
      //
      u[0] = +ct * v[0] - st * v[2];
      u[2] = +st * v[0] + ct * v[2];
      u[3] = +ct * v[3] - st * v[5];
      u[5] = +st * v[3] + ct * v[5];
      u[6] = +ct * v[6] - st * v[8];
      u[8] = +st * v[6] + ct * v[8];
    } else {
      // when a01 sms[5] is the biggest
      // the condition statement shoud pass exactly once for each iteration.
      const double t = 0.5 * atan2(2 * m[5], m[1] - m[0]);
      const double ct = cos(t);
      const double st = sin(t);
      sms[0] = ct * ct * m[0] + st * st * m[1] - 2 * st * ct * m[5];
      sms[1] = ct * ct * m[1] + st * st * m[0] + 2 * st * ct * m[5];
      sms[3] = st * m[4] + ct * m[3];
      sms[4] = ct * m[4] - st * m[3];
      sms[5] = 0; // (ct*ct-st*st)*m[5]+st*ct*(m[0]-m[1]);
      //
      u[0] = +ct * v[0] - st * v[1];
      u[1] = +st * v[0] + ct * v[1];
      u[3] = +ct * v[3] - st * v[4];
      u[4] = +st * v[3] + ct * v[4];
      u[6] = +ct * v[6] - st * v[7];
      u[7] = +st * v[6] + ct * v[7];
    }
  }
  l[0] = scale * sms[0];
  l[1] = scale * sms[1];
  l[2] = scale * sms[2];
  return true;
}

// m = UGV^T
DFM2_INLINE void delfem2::Svd3(
    double U[9],
    double G[3],
    double V[9],
    const double m[9],
    int nitr) {
  // M^TM = VGGV^T
  const double mtm[6] = {
      m[0] * m[0] + m[3] * m[3] + m[6] * m[6],
      m[1] * m[1] + m[4] * m[4] + m[7] * m[7],
      m[2] * m[2] + m[5] * m[5] + m[8] * m[8],
      m[1] * m[2] + m[4] * m[5] + m[7] * m[8],
      m[2] * m[0] + m[5] * m[3] + m[8] * m[6],
      m[0] * m[1] + m[3] * m[4] + m[6] * m[7]};
  double lv[3];
  EigenSym3(V, lv,
            mtm, nitr);
  svd3::SortEigen3(lv, V);
  if (lv[0] < 0) { lv[0] = 0.0; }
  if (lv[1] < 0) { lv[1] = 0.0; }
  if (lv[2] < 0) { lv[2] = 0.0; }
  G[0] = sqrt(lv[0]);
  G[1] = sqrt(lv[1]);
  G[2] = sqrt(lv[2]);

  double u0[3] = {
      m[0] * V[0] + m[1] * V[3] + m[2] * V[6],
      m[3] * V[0] + m[4] * V[3] + m[5] * V[6],
      m[6] * V[0] + m[7] * V[3] + m[8] * V[6]};
  double u1[3] = {
      m[0] * V[1] + m[1] * V[4] + m[2] * V[7],
      m[3] * V[1] + m[4] * V[4] + m[5] * V[7],
      m[6] * V[1] + m[7] * V[4] + m[8] * V[7]};
  double u2[3] = {
      m[0] * V[2] + m[1] * V[5] + m[2] * V[8],
      m[3] * V[2] + m[4] * V[5] + m[5] * V[8],
      m[6] * V[2] + m[7] * V[5] + m[8] * V[8]};

  if (svd3::Det_Mat3(V) < 0) {  // making right hand coordinate
    V[0 * 3 + 2] *= -1;
    V[1 * 3 + 2] *= -1;
    V[2 * 3 + 2] *= -1;
    G[2] *= -1;
    U[0 * 3 + 2] *= -1;
    U[1 * 3 + 2] *= -1;
    U[2 * 3 + 2] *= -1;
  }

  const double sql0 = svd3::SqLength3(u0);
  if (sql0 > 1.0e-20) {
    svd3::Normalize3(u0);
    const double sql1 = svd3::SqLength3(u1);
    if (sql1 < 1.0e-20) {
      u1[0] = 1.0 - fabs(u0[0]);
      u1[1] = 1.0 - fabs(u0[1]);
      u1[2] = 1.0 - fabs(u0[2]);
    } else {
      svd3::Normalize3(u1);
    }
    const double d01 = u0[0] * u1[0] + u0[1] * u1[1] + u0[2] * u1[2];
    u1[0] -= d01 * u0[0];
    u1[1] -= d01 * u0[1];
    u1[2] -= d01 * u0[2];
    svd3::Normalize3(u1);
    double s2[3];
    svd3::Cross3(s2, u0, u1);
    const double d22 = u2[0] * s2[0] + u2[1] * s2[1] + u2[2] * s2[2];
    u2[0] = s2[0];
    u2[1] = s2[1];
    u2[2] = s2[2];
    if (d22 < 0) {
      G[2] *= -1;
    }
  } else {
    u0[0] = 1;
    u0[1] = 0;
    u0[2] = 0;
    u1[0] = 0;
    u1[1] = 1;
    u1[2] = 0;
    u2[0] = 0;
    u2[1] = 0;
    u2[2] = 1;
  }
  U[0] = u0[0];
  U[1] = u1[0];
  U[2] = u2[0];
  U[3] = u0[1];
  U[4] = u1[1];
  U[5] = u2[1];
  U[6] = u0[2];
  U[7] = u1[2];
  U[8] = u2[2];
}

DFM2_INLINE void delfem2::GetRotPolarDecomp(
    double R[9],
    const double am[9],
    int nitr) {
  double U[9], G[3], V[9];
  // UGV^T = am
  Svd3(
      U, G, V,
      am, nitr);
  // R = UV^T
  svd3::MatMatT3(
      R,
      U, V);
}