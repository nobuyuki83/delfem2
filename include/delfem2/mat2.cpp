/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mat2.h"

namespace delfem2 {
namespace mat2 {

template<typename T>
T SquareLength2(const T v[2]) {
  return v[0] * v[0] + v[1] * v[1];
}

template<typename T>
T Length2(const T v[2]) {
  return sqrt(v[0] * v[0] + v[1] * v[1]);
}

template<typename T>
void Normalize2(T w[2]) {
  T l = Length2(w);
  T invl = static_cast<T>(1) / l;
  w[0] *= invl;
  w[1] *= invl;
}

}
}

// -------------------------------------

template<typename T>
void delfem2::MatVec2(
    T w[2],
    const T A[4],
    const T v[2]) {
  w[0] = A[0] * v[0] + A[1] * v[1];
  w[1] = A[2] * v[0] + A[3] * v[1];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MatVec2(float w[2], const float A[4], const float v[2]);
template void delfem2::MatVec2(double w[2], const double A[4], const double v[2]);
#endif

// -------------------------------------

template<typename T>
void delfem2::MatMat2(T AB[4], const T A[4], const T B[4]) {
  AB[0 * 2 + 0] = A[0 * 2 + 0] * B[0 * 2 + 0] + A[0 * 2 + 1] * B[1 * 2 + 0];
  AB[0 * 2 + 1] = A[0 * 2 + 0] * B[0 * 2 + 1] + A[0 * 2 + 1] * B[1 * 2 + 1];
  AB[1 * 2 + 0] = A[1 * 2 + 0] * B[0 * 2 + 0] + A[1 * 2 + 1] * B[1 * 2 + 0];
  AB[1 * 2 + 1] = A[1 * 2 + 0] * B[0 * 2 + 1] + A[1 * 2 + 1] * B[1 * 2 + 1];
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::MatMat2(float AB[4], const float A[4], const float B[4]);
template void delfem2::MatMat2(double AB[4], const double A[4], const double B[4]);
#endif

// --------------------------------------------------------------

DFM2_INLINE bool delfem2::InverseMat2(
    double invB[4],
    const double B[4]) {
  double det = B[0] * B[3] - B[1] * B[2];
  if (fabs(det) < 1.0e-10) return false;
  double invdet = 1.0 / det;
  invB[0] = +invdet * B[3];
  invB[1] = -invdet * B[1];
  invB[2] = -invdet * B[2];
  invB[3] = +invdet * B[0];
  return true;
}

DFM2_INLINE void delfem2::gramian2
    (double AtA[3],
     const double A[4]) {
  AtA[0] = A[0 * 2 + 0] * A[0 * 2 + 0] + A[1 * 2 + 0] * A[1 * 2 + 0];
  AtA[1] = A[0 * 2 + 0] * A[0 * 2 + 1] + A[1 * 2 + 0] * A[1 * 2 + 1];
  AtA[2] = AtA[1];
  AtA[3] = A[0 * 2 + 1] * A[0 * 2 + 1] + A[1 * 2 + 1] * A[1 * 2 + 1];
}

DFM2_INLINE void delfem2::VLVt2
    (double A[4],
     double l0,
     double l1,
     const double V[4]) {
  A[0] = l0 * V[0] * V[0] + l1 * V[1] * V[1];
  A[1] = l0 * V[2] * V[0] + l1 * V[3] * V[1];
  A[2] = l0 * V[0] * V[2] + l1 * V[1] * V[3];
  A[3] = l0 * V[2] * V[2] + l1 * V[3] * V[3];
}

DFM2_INLINE void delfem2::RotationalComponentOfMatrix2(
    double R[4],
    const double M[4]) {
  namespace lcl = delfem2::mat2;
  const double eps = 1.0e-20;
  double A[4];
  {
    double s = fabs(M[0]) + fabs(M[1]) + fabs(M[2]) + fabs(M[3]);
    if (s < 1.0e-10) {
      R[0] = 1.0;
      R[1] = 0.0;
      R[2] = 0.0;
      R[3] = 1.0;
      return;
    }
    double invs = 1.0 / s;
    A[0] = invs * M[0];
    A[1] = invs * M[1];
    A[2] = invs * M[2];
    A[3] = invs * M[3];
  }
  double G[4];
  gramian2(G, A);
  double l0, l1;
  double v0[2], v1[2];
  {
    double b = G[0] + G[3];
    double c = G[0] * G[3] - G[1] * G[2];
    double d = b * b - 4 * c;
    if (d < eps) {
      l0 = 0.5 * b;
      l1 = 0.5 * b;
      v0[0] = 0;
      v0[1] = 1;
      v1[0] = 1;
      v1[1] = 0;
    } else {
      d = sqrt(d);
      l0 = 0.5 * (b + d);
      l1 = 0.5 * (b - d);
      v0[0] = G[1];
      v0[1] = G[3] - l1;
      if (lcl::SquareLength2(v0) > eps) { lcl::Normalize2(v0); }
      v1[0] = G[0] - l0;
      v1[1] = G[2];
      if (lcl::SquareLength2(v1) > eps) { lcl::Normalize2(v1); }
    }
  }
  double V[4] = {v0[0], v1[0], v0[1], v1[1]};
  if (l0 < eps) { l0 = 1; }
  if (l1 < eps) { l1 = 1; }
  double il0 = 1.0 / sqrt(l0);
  double il1 = 1.0 / sqrt(l1);
  double invS[4];
  VLVt2(invS, il0, il1, V);
  MatMat2(R, A, invS);
}


// ------------------------------

namespace delfem2 {

template<typename T>
CMat2<T> operator*(
    T d,
    const CMat2<T> &rhs) {
  CMat2<T> temp = rhs;
  temp *= d;
  return temp;
}

}

namespace delfem2 {

template<typename T>
CMat2<T> operator+(const CMat2<T> &lhs, const CMat2<T> &rhs) {
  CMat2<T> q;
  for (int i = 0; i < 4; ++i) { q.p[i] = lhs.p[i] + rhs.p[i]; }
  return q;
}
#ifdef DFM2_STATIC_LIBRARY
template CMat2d operator + (const CMat2d& lhs, const CMat2d& rhs);
template CMat2f operator + (const CMat2f& lhs, const CMat2f& rhs);
#endif

}

// --------

namespace delfem2 {

template<typename T>
CMat2<T> operator-(
    const CMat2<T> &lhs,
    const CMat2<T> &rhs) {
  return CMat2<T>(
      lhs.p[0] - rhs.p[0],
      lhs.p[1] - rhs.p[1],
      lhs.p[2] - rhs.p[2],
      lhs.p[3] - rhs.p[3]);
}
#ifdef DFM2_STATIC_LIBRARY
template CMat2d operator - (const CMat2d& lhs, const CMat2d& rhs);
template CMat2f operator - (const CMat2f& lhs, const CMat2f& rhs);
#endif

}

// ---------------------------

namespace delfem2 {

template<typename T>
CMat2<T> operator*(
    const CMat2<T> &lhs,
    const CMat2<T> &rhs) {
  CMat2<T> m;
  MatMat2(m.p, lhs.p, rhs.p);
  return m;
}
#ifdef DFM2_STATIC_LIBRARY
template CMat2d operator * (const CMat2d& lhs, const CMat2d& rhs);
template CMat2f operator * (const CMat2f& lhs, const CMat2f& rhs);
#endif

}

// ---------------------------

template<typename T>
void delfem2::polar_decomposition(
    CMat2<T> &R,
    CMat2<T> &S,
    const CMat2<T> &m) {
  auto x = m(0, 0) + m(1, 1);
  auto y = m(1, 0) - m(0, 1);
  auto scale = static_cast<T>(1) / std::sqrt(x * x + y * y);
  auto c = x * scale;
  auto s = y * scale;
  R(0, 0) = c;
  R(0, 1) = -s;
  R(1, 0) = s;
  R(1, 1) = c;
  S = R.transpose() * m;
}

template<typename T>
void delfem2::svd(CMat2<T> &U,
                  CMat2<T> &sig,
                  CMat2<T> &V,
                  const CMat2<T> &m) {
  CMat2<T> S;
  polar_decomposition(U, S, m);
  T c, s;
  if (std::abs(S(0, 1)) < 1e-6f) {
    sig = S;
    c = 1;
    s = 0;
  } else {
    auto tao = 0.5f * (S(0, 0) - S(1, 1));
    auto w = std::sqrt(tao * tao + S(0, 1) * S(0, 1));
    auto t = tao > 0 ? S(0, 1) / (tao + w) : S(0, 1) / (tao - w);
    c = 1.0f / std::sqrt(t * t + 1);
    s = -t * c;
    sig(0, 0) = (c * c) * S(0, 0) - 2 * c * s * S(0, 1) + (s * s) * S(1, 1);
    sig(1, 1) = (s * s) * S(0, 0) + 2 * c * s * S(0, 1) + (c * c) * S(1, 1);
  }
  if (sig(0, 0) < sig(1, 1)) {
    std::swap(sig(0, 0), sig(1, 1));
    V(0, 0) = -s;
    V(0, 1) = -c;
    V(1, 0) = c;
    V(1, 1) = -s;
  } else {
    V(0, 0) = c;
    V(0, 1) = -s;
    V(1, 0) = s;
    V(1, 1) = c;
  }
  V.transposeInPlace();
  U = U * V;
}
