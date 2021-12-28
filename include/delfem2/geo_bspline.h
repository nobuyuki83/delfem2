/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @detail The order of dependency in delfem2:
 * line < ray < edge < polyline < quadratic < cubic < bspline << plane < tri < quad
 */

#ifndef DFM2_CURVE_BSPLINE_H
#define DFM2_CURVE_BSPLINE_H

#include <cassert>
#include <cstdio>
#include <vector>
#include <algorithm>

namespace delfem2 {

/**
 *
 * @tparam SCALAR float or double
 * @param[out] coeff array of coefficients
 * @param[in] idx_segment
 * @param[in] num_segment
 * @detail (i-th point's weight) = coeff[i][0] + coeff[i][1] * t + coeff[i][2] * t * t
 */
template<typename SCALAR>
void CoefficientsOfOpenUniformBSpline_Quadratic(
  SCALAR coeff[3][3],
  int idx_segment,
  int num_segment) {

  const int k0 = std::clamp<int>(idx_segment - 1, 0, num_segment) - idx_segment;
  const int k1 = std::clamp<int>(idx_segment + 0, 0, num_segment) - idx_segment;
  const int k2 = std::clamp<int>(idx_segment + 1, 0, num_segment) - idx_segment;
  const int k3 = std::clamp<int>(idx_segment + 2, 0, num_segment) - idx_segment;
  assert(-1 <= k0 && k0 <= k1 && k1 <= k2 && k2 <= k3 && k3 <= 2);

  const SCALAR c00 = (k2 == k1) ? 0.0 : 1. / static_cast<SCALAR>(k2 - k1);
  const SCALAR c10 = (k2 == k0) ? 0.0 : 1. / static_cast<SCALAR>(k2 - k0);
  const SCALAR c11 = (k3 == k1) ? 0.0 : 1. / static_cast<SCALAR>(k3 - k1);

  const SCALAR d22 = c00 * c10;
  const SCALAR d02 = c10 * c00;
  const SCALAR d13 = c11 * c00;
  const SCALAR d11 = c00 * c11;

  coeff[0][0] = k2 * k2 * d22;
  coeff[0][1] = -2 * k2 * d22;
  coeff[0][2] = d22;
  //
  coeff[1][0] = -k0 * k2 * d02 - k1 * k3 * d13;
  coeff[1][1] = (k0 + k2) * d02 + (k1 + k3) * d13;
  coeff[1][2] = -(d02 + d13);
  //
  coeff[2][0] = k1 * k1 * d11;
  coeff[2][1] = -2 * k1 * d11;
  coeff[2][2] = d11;
}

/**
 * Quadratic B-Spline with "open and uniform knot vector"
 * knot vector = [0,0,0,1,2,3,...,N-1,N,N,N] where N is the nubmer of segments, which is poly.size()-2
 * @param[in] t parameter of curve that takes [0,N]
 * @param[in] poly position of the the control points
 * @return sampled point
 */
template<typename VEC>
VEC Sample_QuadraticBsplineCurve(
  typename VEC::Scalar t,
  const std::vector<VEC> &poly) {
  using SCALAR = typename VEC::Scalar;

  const int num_segment = poly.size() - 2;
  const int idx_segment = static_cast<int>(t) + (t == num_segment ? -1 : 0);
  assert(idx_segment >= 0 && idx_segment < num_segment);

  t -= idx_segment;
  assert(t >= 0 && t <= 1);

  SCALAR coeff[3][3];
  CoefficientsOfOpenUniformBSpline_Quadratic(coeff, idx_segment, num_segment);

  const SCALAR w0 = coeff[0][0] + coeff[0][1] * t + coeff[0][2] * t * t;
  const SCALAR w1 = coeff[1][0] + coeff[1][1] * t + coeff[1][2] * t * t;
  const SCALAR w2 = coeff[2][0] + coeff[2][1] * t + coeff[2][2] * t * t;

  assert(fabs(w0 + w1 + w2 - 1.) < 1.0e-10);
  assert(w0 >= 0 && w1 >= 0 && w2 >= 0);
  return poly[idx_segment] * w0 + poly[idx_segment + 1] * w1 + poly[idx_segment + 2] * w2;
}

/**
 * at parameter "t" the position is evaluated as
 * const SCALAR w0 = coeff[0][0] + coeff[0][1] * t + coeff[0][2] * t * t;
 * const SCALAR w1 = coeff[1][0] + coeff[1][1] * t + coeff[1][2] * t * t;
 * const SCALAR w2 = coeff[2][0] + coeff[2][1] * t + coeff[2][2] * t * t;
 * return points[0] * w0 + points[1] * w1 + points[2] * w2;
 */
template<typename VEC>
double Length_ParametricCurve_Quadratic(
  const double coeff[3][3],
  const VEC points[3]) {
  const VEC p0 = coeff[0][1] * points[0] + coeff[1][1] * points[1] + coeff[2][1] * points[2];
  const VEC p1 = coeff[0][2] * points[0] + coeff[1][2] * points[1] + coeff[2][2] * points[2];
  // dl = sqrt(coe[0] + coe[1]*t + coe[2]*t^2) dt
  const double coe[3] = {
    p0.squaredNorm(),
    4*p0.dot(p1),
    4*p1.squaredNorm()};
  auto intt = [&coe](double t) -> double {
    const double tmp0 = (coe[1] + 2 * coe[2] * t) * std::sqrt(coe[0] + t * (coe[1] + coe[2] * t));
    const double tmp3 = sqrt(coe[2]);
    const double tmp1 = coe[1] + 2 * coe[2] * t + 2 * tmp3 * sqrt(coe[0] + t * (coe[1] + coe[2] * t));
    const double tmp2 = (coe[1] * coe[1] - 4 * coe[0] * coe[2]) * std::log(tmp1);
    return tmp0 / (4. * coe[2]) - tmp2 / (8. * coe[2] * tmp3);
  };
  return intt(1) - intt(0);
}

// above: quadratic
// ==========================================
// below: cubic

/**
 *
 * @tparam SCALAR float or double
 * @param[out] coeff array of coefficients
 * @param[in] idx_segment
 * @param[in] num_segment
 * @detail (i-th point's weight) = coeff[i][0] + coeff[i][1] * t + coeff[i][2] * t * t + coeff[i][3] * t * t * t
 */
template<typename SCALAR>
void CoefficientsOfOpenUniformBSpline_Cubic(
  SCALAR coeff[4][4],
  int idx_segment,
  int num_segment) {

  // knot vector for this segement
  const int k0 = std::clamp<int>(idx_segment - 2, 0, num_segment) - idx_segment;
  const int k1 = std::clamp<int>(idx_segment - 1, 0, num_segment) - idx_segment;
  const int k2 = std::clamp<int>(idx_segment + 0, 0, num_segment) - idx_segment;
  const int k3 = std::clamp<int>(idx_segment + 1, 0, num_segment) - idx_segment;
  const int k4 = std::clamp<int>(idx_segment + 2, 0, num_segment) - idx_segment;
  const int k5 = std::clamp<int>(idx_segment + 3, 0, num_segment) - idx_segment;
  assert(-2 <= k0 && k0 <= k1 && k1 <= k2 && k2 <= k3 && k3 <= k4 && k4 <= k5 && k5 <= 3);

  const SCALAR c32 = (k3 == k2) ? 0 : 1 / static_cast<SCALAR>(k3 - k2);
  const SCALAR c31 = (k3 == k1) ? 0 : 1 / static_cast<SCALAR>(k3 - k1);
  const SCALAR c42 = (k4 == k2) ? 0 : 1 / static_cast<SCALAR>(k4 - k2);
  const SCALAR c30 = (k3 == k0) ? 0 : 1 / static_cast<SCALAR>(k3 - k0);
  const SCALAR c41 = (k4 == k1) ? 0 : 1 / static_cast<SCALAR>(k4 - k1);
  const SCALAR c52 = (k5 == k2) ? 0 : 1 / static_cast<SCALAR>(k5 - k2);

  {
    // (k3-t) * (k3-t) * (k3-t)
    const SCALAR d333 = c30 * c31 * c32;
    coeff[0][0] = k3 * k3 * k3 * d333;
    coeff[0][1] = -3 * k3 * k3 * d333;
    coeff[0][2] = +3 * k3 * d333;
    coeff[0][3] = -d333;
  }
  {
    // (t-k0) * (k3-t) * (k3-t)
    // (k4-t) * (t-k1) * (k3-t)
    // (k4-t) * (k4-t) * (t-k2)
    const SCALAR d033 = c30 * c31 * c32;
    const SCALAR d413 = c41 * c31 * c32;
    const SCALAR d442 = c41 * c42 * c32;
    coeff[1][0] = -k0 * k3 * k3 * d033 - k4 * k1 * k3 * d413 - k4 * k4 * k2 * d442;
    coeff[1][1] =
      +(2 * k0 * k3 + k3 * k3) * d033
        + (k4 * k1 + k1 * k3 + k3 * k4) * d413
        + (k4 * k4 + 2 * k4 * k2) * d442;
    coeff[1][2] = -(k0 + 2 * k3) * d033 - (k4 + k1 + k3) * d413 - (2 * k4 + k2) * d442;
    coeff[1][3] = d033 + d413 + d442;
  }
  {
    // (t-k1) * (t-k1) * (k3-t) / (k4-k1) / (k3-k1) / (k3-k2)
    // (t-k1) * (k4-t) * (t-k2) / (k4-k1) / (k4-k2) / (k3-k2)
    // (k5-t) * (t-k2) * (t-k2) / (k5-k2) / (k4-k2) / (k3-k2)
    const SCALAR d113 = c41 * c31 * c32;
    const SCALAR d142 = c41 * c42 * c32;
    const SCALAR d522 = c52 * c42 * c32;
    coeff[2][0] = k1 * k1 * k3 * d113 + k1 * k4 * k2 * d142 + k5 * k2 * k2 * d522;
    coeff[2][1] =
      -(2 * k1 * k3 + k1 * k1) * d113
        - (k1 * k4 + k4 * k2 + k2 * k1) * d142
        - (2 * k5 * k2 + k2 * k2) * d522;
    coeff[2][2] = (2 * k1 + k3) * d113 + (k1 + k4 + k2) * d142 + (k5 + 2 * k2) * d522;
    coeff[2][3] = -d113 - d142 - d522;
  }
  {
    // (t-k2) * (t-k2) * (t-k2)
    const SCALAR d222 = c52 * c42 * c32;
    coeff[3][0] = -k2 * k2 * k2 * d222;
    coeff[3][1] = +3 * k2 * k2 * d222;
    coeff[3][2] = -3 * k2 * d222;
    coeff[3][3] = d222;
  }

  /*
  const double w00 = safe_divide(k3 - t, k3 - k2);
  const double w10 = safe_divide(k3 - t, k3 - k1);
  const double w11 = safe_divide(k4 - t, k4 - k2);
  const double w20 = safe_divide(k3 - t, k3 - k0);
  const double w21 = safe_divide(k4 - t, k4 - k1);
  const double w22 = safe_divide(k5 - t, k5 - k2);

  const double w0 = w20 * w10 * w00;
  const double w1 = (1 - w20) * w10 * w00 + w21 * (1 - w10) * w00 + w21 * w11 * (1 - w00);
  const double w2 = (1 - w21) * (1 - w10) * w00 + (1 - w21) * w11 * (1 - w00) + w22 * (1 - w11) * (1 - w00);
  const double w3 = (1 - w22) * (1 - w11) * (1 - w00);
   */
}

/**
 * Quadratic B-Spline with "open and uniform knot vector"
 * knot vector = [0,0,0,1,2,3,...,N-1,N,N,N] where N is poly.size()-2
 * @param t parameter of curve that takes [0,1]
 * @param poly position of the the control points
 * @return sampled point
 */
template<typename VEC>
VEC Sample_CubicBsplineCurve(
  typename VEC::Scalar t,
  const std::vector<VEC> &poly) {
  using SCALAR = typename VEC::Scalar;

  const int num_segment = poly.size() - 3;
  const int idx_segment = static_cast<int>(t) + (t == num_segment ? -1 : 0);
  assert(idx_segment >= 0 && idx_segment < num_segment);

  t -= idx_segment;
  assert(t >= 0 && t <= 1);

  SCALAR coeff[4][4];
  CoefficientsOfOpenUniformBSpline_Cubic(coeff, idx_segment, num_segment);

  SCALAR v0 = coeff[0][0] + coeff[0][1] * t + coeff[0][2] * t * t + coeff[0][3] * t * t * t;
  SCALAR v1 = coeff[1][0] + coeff[1][1] * t + coeff[1][2] * t * t + coeff[1][3] * t * t * t;
  SCALAR v2 = coeff[2][0] + coeff[2][1] * t + coeff[2][2] * t * t + coeff[2][3] * t * t * t;
  SCALAR v3 = coeff[3][0] + coeff[3][1] * t + coeff[3][2] * t * t + coeff[3][3] * t * t * t;

  assert(fabs(v0 + v1 + v2 + v3 - 1.) < 1.0e-10);
  assert(v0 >= -1.0e-10 && v1 >= -1.0e-10 && v2 >= -1.0e-10 && v3 >= -1.0e-10);
  return poly[idx_segment] * v0 + poly[idx_segment + 1] * v1 + poly[idx_segment + 2] * v2 + poly[idx_segment + 3] * v3;
}


// above: cubic
// ====================================
// below: general

void FlatKnot(
  std::vector<double> &aKnotFlat,
  const std::vector<int> &aKnotMulti,
  const std::vector<double> &aKnot) {
  assert(aKnot.size() == aKnotMulti.size());
  aKnotFlat.clear();
  for (size_t ik = 0; ik < aKnot.size(); ++ik) {
    for (int iik = 0; iik < aKnotMulti[ik]; ++iik) {
      aKnotFlat.push_back(aKnot[ik]);
    }
  }
}

template<typename T>
T DeBoorBSpline(
  double u,
  int ndegree,
  const std::vector<T> &aCP,
  const std::vector<double> &aKnot) {
  assert(ndegree > 0);
  assert(aKnot.size() == aCP.size() + ndegree + 1);
  const double eps = 1.0e-10;
  //
  unsigned int iks;
  {
    for (iks = ndegree; iks < aKnot.size() - ndegree; ++iks) {
      double u0 = aKnot[iks];
      double u1 = aKnot[iks + 1];
      if (u >= u0 - eps && u <= u1 + eps) { break; }
    }
  }
  std::vector<T> aTmp;
  {
    aTmp.reserve(ndegree + 1);
    for (unsigned int ik = iks - ndegree; ik < iks + 1; ++ik) {
      assert(ik >= 0);
      aTmp.push_back(aCP[ik]);
    }
  }
  for (int r = 0; r < ndegree; ++r) {
    for (int j = 0; j < ndegree - r; ++j) {
      double u0 = aKnot[j + iks - ndegree + 1 + r];
      double u1 = aKnot[j + iks + 1];
//      assert(u>=u0-eps && u<u1+eps);
      double a = (u - u0) / (u1 - u0);
      aTmp[j] = (1 - a) * aTmp[j] + a * aTmp[j + 1];
    }
  }
  return aTmp[0];
}

template<typename T>
void SampleBSpline(
  std::vector<T> &polyline0,
  const int nsmpl,
  const int ndegree,
  const std::vector<double> &aKnotFlat,
  const std::vector<T> &aCtrlPoint) {
  polyline0.clear();
  double u0 = aKnotFlat[0];
  double u1 = aKnotFlat[aKnotFlat.size() - 1];
  for (int i = 0; i < nsmpl + 1; ++i) {
    double u = (double) i * (u1 - u0) / nsmpl + u0;
    T p = DeBoorBSpline(u, ndegree, aCtrlPoint, aKnotFlat);
    polyline0.push_back(p);
  }
}

template<typename VEC, int norderplus1>
VEC Sample_BsplineCurveSegment(
  double t,
  const int knot[norderplus1 * 2],
  const VEC polyloc[norderplus1]) {

  assert(t >= 0 && t <= 1);

  const auto safe_divide = [](double a, int b) {
    return (b == 0) ? 0. : a / static_cast<double>(b);
  };

  double coefficient[norderplus1 + 1][norderplus1];
  for (int i = 0; i <= norderplus1; i++) {
    for (int m = 0; m < norderplus1; m++) {
      coefficient[i][m] = 0.;
    }
  }
  coefficient[norderplus1 - 1][0] = 1.;

  for (int i = 2; i < norderplus1 + 1; i++) {
    for (int j = 0; j < norderplus1; j++) {
      const double coe[4] = {
        -safe_divide(knot[j], knot[j + i - 1] - knot[j]),
        +safe_divide(knot[j + i], knot[j + i] - knot[j + 1]),
        +safe_divide(1, knot[j + i - 1] - knot[j]),
        -safe_divide(1, knot[j + i] - knot[j + 1])};
      for (int m = norderplus1 - 1; m > 0; m--) {
        coefficient[j][m]
          = coe[0] * coefficient[j][m]
          + coe[1] * coefficient[j + 1][m]
          + coe[2] * coefficient[j][m - 1]
          + coe[3] * coefficient[j + 1][m - 1];
      }
      coefficient[j][0]
        = coe[0] * coefficient[j][0]
        + coe[1] * coefficient[j + 1][0];
    }
  }

  double weight[norderplus1 + 1];
  for (int i = 0; i < norderplus1 + 1; i++) {
    weight[i] = coefficient[i][norderplus1 - 1];
    for (int j = 1; j < norderplus1; ++j) {
      weight[i] = weight[i] * t + coefficient[i][norderplus1 - 1 - j];
    }
  }

  VEC p(0, 0);
  for (int i = 0; i < norderplus1; i++) {
    p += polyloc[i] * weight[i];
  }
  return p;
}

template<typename VEC, int norderplus1>
VEC Sample_BsplineCurve(
  double t,
  const std::vector<VEC> &poly) {

  const int N = poly.size() + 1 - norderplus1; // N=max{knot vector}
  const int index = static_cast<int>(t) + (t == N ? -1 : 0);
  assert(index >= 0 && index < poly.size() + norderplus1 - 1);

  int knot[norderplus1 * 2];
  for (int i = 0; i < norderplus1 * 2; i++) {
    knot[i] = std::clamp<int>(index + i - (norderplus1 - 1), 0, N) - index;
  }
  //
  VEC polyloc[norderplus1];
  for (int i = 0; i < norderplus1; i++) {
    polyloc[i] = poly[index + i];
  }
  //
  return Sample_BsplineCurveSegment<VEC, norderplus1>(t - index, knot, polyloc);
}

template<typename VEC, int norderplus1>
VEC Sample_BsplineDerivative(
  double t,
  const std::vector<VEC> &poly) {
  const int N = poly.size() + 1 - norderplus1;

  const auto knot_generator = [&](int i) -> double {
    return std::clamp(i - norderplus1 + 1, 0, N) / static_cast<double>(N);
  };
  const int index = t + (t == N ? -1 : 0);
  //
  VEC Q_poly[norderplus1 - 1];
  for (int i = index; i < index + norderplus1 - 1; i++) {
    double w = (norderplus1 - 1) / (knot_generator(i + norderplus1) - knot_generator(i + 1));
    Q_poly[i - index] = w * (poly[i + 1] - poly[i]);
  }
  //
  int knot[norderplus1 * 2 - 2];
  for (int i = 0; i < norderplus1 * 2 - 2; i++) {
    knot[i] = std::clamp(index + i - norderplus1 + 2, 0, N) - index;
  }
  //
  return Sample_BsplineCurveSegment<VEC, norderplus1 - 1>(t - index, knot, Q_poly) / N;
}

}

#endif /* DFM2_CURVE_BSPLINE_H */
