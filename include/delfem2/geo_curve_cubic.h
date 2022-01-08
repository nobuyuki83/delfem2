/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @detail The order of dependency in delfem2:
 * aabb ->
 * line -> ray -> edge -> polyline ->
 * curve_quadratic -> curve_cubic -> curve_ndegree ->
 * plane -> tri -> quad
 */

#ifndef DFM2_CURVE_CUBIC_BEZIER_H
#define DFM2_CURVE_CUBIC_BEZIER_H

#include <cassert>
#include <cstdio>
#include <vector>
#include <algorithm>

#include "quadrature.h"
#include "polynomial_root.h"

namespace delfem2 {

template<typename VEC>
VEC PointOnCubicBezierCurve(
  typename VEC::Scalar t,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3,
  const VEC &p4) {
  auto tp = 1 - t;
  return t * t * t * p4
    + 3 * t * t * tp * p3
    + 3 * t * tp * tp * p2
    + tp * tp * tp * p1;
}

template<typename T>
T Tangent_CubicBezierCurve(
  double t,
  const T &p1, const T &p2, const T &p3, const T &p4) {
  double tp = 1.0 - t;
  return 3 * t * t * p4
    + 3 * t * (2 - 3 * t) * p3
    + 3 * tp * (1 - 3 * t) * p2
    - 3 * tp * tp * p1;
}

template<class VEC>
void Polyline_BezierCubic(
  std::vector<VEC> &aP,
  const unsigned int n,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3,
  const VEC &p4) {
  using SCALAR = typename VEC::Scalar;
  aP.resize(n);
  for (unsigned int i = 0; i < n; ++i) {
    const SCALAR t = static_cast<SCALAR>(i) / static_cast<SCALAR>(n - 1);
    aP[i] = PointOnCubicBezierCurve(
      t,
      p1, p2, p3, p4);
  }
}

/*
template<class VEC>
void Polyline_CubicBezierCurve(
  std::vector<VEC>& aP,
  const int n,
  const std::vector<VEC>& aCP) {
  int ns = (int) (aCP.size() / 3);
  aP.resize(ns * n + 1);
  for (int is = 0; is < ns; is++) {
    for (int i = 0; i < n; i++) {
      double t = (double) i / n;
      aP[is * n + i] = PointOnCubicBezierCurve(
        t,
        aCP[is * 3 + 0],
        aCP[is * 3 + 1],
        aCP[is * 3 + 2],
        aCP[is * 3 + 3]);
    }
  }
  aP[ns * n] = aCP[ns * 3];
}
 */

// -----------------------------

template<typename T>
bool getParameterCubicBezier_IntersectionWithPlane(
  double &t,
  const T &org,
  const T &nrm,
  const T &p1,
  const T &p2,
  const T &p3,
  const T &p4) {
  double h1 = (p1 - org).dot(nrm);
  double h2 = (p2 - org).dot(nrm);
  double h3 = (p3 - org).dot(nrm);
  double h4 = (p4 - org).dot(nrm);
  double ref = fabs(h1) + fabs(h2) + fabs(h3) + fabs(h4);
  double eps = 1.0e-5;
  if (fabs(h1) < ref * eps && fabs(h4) < ref * eps) return false;
  if (h1 * h4 > 0) return false;
  if (fabs(h1) < ref * eps) {
    t = 0.0;
    return true;
  }
  if (fabs(h4) < ref * eps) {
    t = 1.0;
    return true;
  }
  t = 0.5;
  for (int itr = 0; itr < 10; ++itr) {
    double tp = 1.0 - t;
    double h = t * t * t * h4
      + 3 * t * t * tp * h3
      + 3 * t * tp * tp * h2
      + tp * tp * tp * h1;
    if (fabs(h) < 1.0e-6 * ref) return true;
    double dh = 3 * t * t * h4
      + 3 * t * (2 - 3 * t) * h3
      + 3 * tp * (1 - 3 * t) * h2
      - 3 * tp * tp * h1;
    t -= (h / dh);
  }
  return false;
}

template<typename T>
void getCubicBezierCurve(
  const int n,
  std::vector<T> &aP,
  const std::vector<T> &aCP) {
  int ns = (int) (aCP.size() / 3);
  aP.resize(ns * n + 1);
  for (int is = 0; is < ns; is++) {
    for (int i = 0; i < n; i++) {
      double t = (double) i / n;
      aP[is * n + i] = getPointCubicBezierCurve(
        t,
        aCP[is * 3 + 0],
        aCP[is * 3 + 1],
        aCP[is * 3 + 2],
        aCP[is * 3 + 3]);
    }
  }
  aP[ns * n] = aCP[ns * 3];
}

template<typename VEC>
typename VEC::Scalar Length_CubicBezierCurve_Quadrature(
  const VEC &p0,  // end point
  const VEC &p1,  // the control point next tp p0
  const VEC &p2,  // the control point next to p1
  const VEC &p3,  // another end point
  int gauss_order)  // order of Gaussian quadrature to use
{
  assert(gauss_order < 6);
  using SCALAR = typename VEC::Scalar;
  SCALAR totalLength = 0;
  const unsigned int iw0 = kNumIntegrationPoint_GaussianQuadrature[gauss_order];
  const unsigned int iw1 = kNumIntegrationPoint_GaussianQuadrature[gauss_order + 1];
  for (unsigned int i = iw0; i < iw1; i++) {
    double t = (kPositionWeight_GaussianQuadrature<double>[i][0] + 1) / 2;
    double w = kPositionWeight_GaussianQuadrature<double>[i][1];
    const VEC dt = 3 * (1 - t) * (1 - t) * (p1 - p0)
      + 6 * (1 - t) * t * (p2 - p1)
      + 3 * t * t * (p3 - p2);
    totalLength += dt.norm() * w;
  }
  return totalLength / 2;
}

// ======================================

template<typename VEC>
typename VEC::Scalar Nearest_CubicBezierCurve(
  const VEC &q,
  const VEC &p0,  // end point
  const VEC &p1,  // the control point
  const VEC &p2,
  const VEC &p3,
  unsigned int num_samples,
  unsigned int num_newton_itr) {   // another end point

  using SCALAR = typename VEC::Scalar;

  const VEC a = -p0 + 3 * p1 - 3 * p2 + p3;
  const VEC b = 3 * p0 - 6 * p1 + 3 * p2;
  const VEC c = -3 * p0 + 3 * p1;
  const VEC d = p0 - q;

  SCALAR t = 0, dist_min = (p0 - q).norm();
  for (unsigned int i = 1; i < num_samples + 1; i++) {
    SCALAR t0 = static_cast<SCALAR>(i) / static_cast<SCALAR>(num_samples);
    SCALAR dist0 = (a * (t0 * t0 * t0) + b * (t0 * t0) + c * t0 + d).norm();
    if (dist0 < dist_min) {
      dist_min = dist0;
      t = t0;
    }
  }

  const SCALAR s0 = 2 * c.dot(d);
  const SCALAR s1 = 2 * c.squaredNorm() + 4 * b.dot(d);
  const SCALAR s2 = 6 * b.dot(c) + 6 * a.dot(d);
  const SCALAR s3 = 4 * b.squaredNorm() + 8 * a.dot(c);
  const SCALAR s4 = 10 * a.dot(b);
  const SCALAR s5 = 6 * a.squaredNorm();
  const SCALAR u0 = s1;
  const SCALAR u1 = 2 * s2;
  const SCALAR u2 = 3 * s3;
  const SCALAR u3 = 4 * s4;
  const SCALAR u4 = 5 * s5;

  for (unsigned int itr = 0; itr < num_newton_itr; ++itr) {
    SCALAR t0 = t;
    SCALAR dw = s0 + t * (s1 + t * (s2 + t * (s3 + t * (s4 + t * s5))));
    SCALAR ddw = u0 + t * (u1 + t * (u2 + t * (u3 + t * u4)));
    t -= dw / ddw;
    t = (t < 0) ? 0 : t;
    t = (t > 1) ? 1 : t;
    SCALAR dist0 = (a * (t * t * t) + b * (t * t) + c * t + d).norm();
    if (dist0 > dist_min) {
      t = t0;
      break;
    }   // winding back
    dist_min = dist0;
  }
  return t;
}

// ------------------------------

/**
 * @return the parameter of cubic bezier curve nearest to q
 */
template<typename VEC>
typename VEC::Scalar Nearest_CubicBezierCurve_Strum(
  const VEC &q,
  const VEC &p0,  // end point
  const VEC &p1,  // the control point
  const VEC &p2,
  const VEC &p3,
  unsigned int num_bisection = 15) {
  // Precompute coefficients
  // p = at^3 + bt^2 + ct + d
  const VEC a = -p0 + 3 * p1 - 3 * p2 + p3;
  const VEC b = 3 * p0 - 6 * p1 + 3 * p2;
  const VEC c = -3 * p0 + 3 * p1;
  const VEC d = p0 - q;

  // Derivative of squared distance function
  // We use the squared distance because it is easier to find its derivative
  const double coe[6] = {
    c.dot(d),
    c.squaredNorm() + 2 * b.dot(d),
    3 * (b.dot(c) + a.dot(d)),
    2 * b.squaredNorm() + 4 * a.dot(c),
    5 * a.dot(b),
    3 * a.squaredNorm()};

  auto roots = RootsOfPolynomial<6>(coe, num_bisection);
  roots.push_back(1); // check t=1
  double best_t = 0., best_dist = d.squaredNorm(); // check t=0
  for (double t: roots) {
    double dist0 = (d + t * (c + t * (b + t * a))).squaredNorm();
    if (dist0 > best_dist) { continue; }
    best_t = t;
    best_dist = dist0;
  }
  return best_t;
}

// ----------------------------------------

template<typename VEC>
auto Area_CubicBezierCurve2(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3) -> typename VEC::Scalar {
  using T = typename VEC::Scalar;
  const T tmp0 = -3 * p2[0] * p0[1] - p3[0] * p0[1]
    - 3 * p2[0] * p1[1]
    - 3 * p3[0] * p1[1]
    - 6 * p3[0] * p2[1]
    + 6 * p2[0] * p3[1]
    + 10 * p3[0] * p3[1]
    + 3 * p1[0] * (-2 * p0[1] + p2[1] + p3[1])
    + p0[0] * (-10 * p0[1] + 6 * p1[1] + 3 * p2[1] + p3[1]);
  return (p0[0] * p0[1]) / 2 - (p3[0] * p3[1]) / 2 + tmp0 / 20;
}

/**
 * @tparam ndim dimension of the geometry
 * @return min_x, min_y, (min_z), max_x, max_y, (max_z)
 */
template<typename VEC>
auto AABB_CubicBezierCurve(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3) -> std::array<typename VEC::Scalar, VEC::SizeAtCompileTime * 2> {
  using SCALAR = typename VEC::Scalar;
  constexpr int ndim = VEC::SizeAtCompileTime;

  const double EPSILON = 1.0e-10;

  const VEC a = -p0 + 3 * p1 - 3 * p2 + p3;
  const VEC b = 3 * p0 - 6 * p1 + 3 * p2;
  const VEC c = -3 * p0 + 3 * p1;
  const VEC d = p0;

  std::array<SCALAR, ndim * 2> res;
  for (int i = 0; i < ndim; i++) {
    if (fabs(a[i]) < EPSILON) {
      SCALAR r = c[i] / (-2 * b[i]);
      r = std::clamp<SCALAR>(r, 0, 1);
      const SCALAR e = d[i] + r * (c[i] + r * (b[i] + r * a[i]));
      res[i] = std::min({p0[i], p3[i], e});
      res[i + ndim] = std::max({p0[i], p3[i], e});
    } else {
      const SCALAR det = std::sqrt(b[i] * b[i] - 3 * a[i] * c[i]);
      SCALAR r0 = (-b[i] + det) / (3 * a[i]);
      SCALAR r1 = (-b[i] - det) / (3 * a[i]);
      r0 = std::clamp<SCALAR>(r0, 0, 1);
      r1 = std::clamp<SCALAR>(r1, 0, 1);
      const SCALAR e0 = d[i] + r0 * (c[i] + r0 * (b[i] + r0 * a[i]));
      const SCALAR e1 = d[i] + r1 * (c[i] + r1 * (b[i] + r1 * a[i]));
      res[i] = std::min({p0[i], p3[i], e0, e1});
      res[i + ndim] = std::max({p0[i], p3[i], e0, e1});
    }
  }
  return res;
}

/**
 * @return 4 control points of the first segment followed by 4 control points of the second segment (in original sequence)
 */
template<typename VEC>
std::array<VEC, 8> Split_CubicBezierCurve(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3,
  typename VEC::Scalar t) {
  // p_i^j = p_i^{j-1} * (1 - t0) + p_{i+1}^{j-1} * t0
  auto mix = [&t](const VEC &q0, const VEC &q1) {
    return q0 * (1 - t) + q1 * t;
  };

  std::array<VEC, 8> res;
  res[0] = p0;
  res[1] = mix(p0, p1); // p01
  VEC p11 = mix(p1, p2);
  res[2] = mix(res[1], p11); // p02
  res[6] = mix(p2, p3); // p21
  res[7] = p3;
  res[5] = mix(p11, res[6]); // p12
  res[3] = res[4] = mix(res[2], res[5]); // p03
  return res;
}

template<typename VEC>
double Length_CubicBezierCurve_QuadratureSubdivision(
  const VEC &p0,  // end point
  const VEC &p1,  // the control point next tp p0
  const VEC &p2,  // the control point next to p1
  const VEC &p3,  // another end point
  double tolerance = 1E-8,
  int maxDepth = 8) {
  using SCALAR = typename VEC::Scalar;

  // derivative of cubic Bezier
  auto cubicBezierdt = [&](SCALAR t) {
    return 3 * (1 - t) * (1 - t) * (p1 - p0) + 6 * (1 - t) * t * (p2 - p1) + 3 * t * t * (p3 - p2);
  };

  constexpr unsigned int np_gauss = 7;
  constexpr unsigned int np_kronrod = 15;
  constexpr unsigned int noff = kNumIntegrationPoint_GaussianQuadrature[np_gauss - 1];

  // integrate using gaussian and gauss-kronrod quadrature at the same time
  SCALAR length_gauss = 0.;
  SCALAR length_kronrod = 0.;
  for (unsigned int iw = 0; iw < np_gauss; ++iw) { // shared part
    const double x0 = kPositionWeight_GaussianQuadrature<double>[noff + iw][0];
    const double w0 = kPositionWeight_GaussianQuadrature<double>[noff + iw][1];
    const double w1 = kPositionWeight_GaussKronrodQuadrature<double>[iw][1];
    assert(fabs(x0 - kPositionWeight_GaussKronrodQuadrature<double>[iw][0]) < 1.0e-8);
    const SCALAR nodeValue = cubicBezierdt((x0 + 1) / 2).norm();
    length_gauss += nodeValue * w0 / 2;
    length_kronrod += nodeValue * w1 / 2;
  }

  // kronrod-only terms
  for (unsigned int iw = np_gauss; iw < np_kronrod; ++iw) {
    const double x0 = kPositionWeight_GaussKronrodQuadrature<double>[iw][0];
    const double w1 = kPositionWeight_GaussKronrodQuadrature<double>[iw][1];
    const SCALAR nodeValue = cubicBezierdt((x0 + 1) / 2).norm();
    length_kronrod += nodeValue * w1 / 2;
  }

  if (std::abs(length_kronrod - length_gauss) < std::abs(length_kronrod) * tolerance) {
    return length_kronrod;
  } else {
    if (maxDepth == 1) {
      return length_kronrod;
    } else {
      std::array<VEC, 8> subdiv = Split_CubicBezierCurve(p0, p1, p2, p3, 0.5);
      double len0 = Length_CubicBezierCurve_QuadratureSubdivision<VEC>(
        subdiv[0], subdiv[1], subdiv[2], subdiv[3], tolerance, maxDepth - 1);
      double len1 = Length_CubicBezierCurve_QuadratureSubdivision<VEC>(
        subdiv[4], subdiv[5], subdiv[6], subdiv[7], tolerance, maxDepth - 1);
      return len0 + len1;
    }
  }

  return -1; // suppress compiler warning
}

// above: Bezier
// ===================================================
// below: bspline

/**
 *
 * @tparam SCALAR float or double
 * @param[out] coeff array of coefficients
 * @param[in] idx_segment
 * @param[in] num_segment
 * @detail (i-th point's weight) = coeff[i][0] + coeff[i][1] * t + coeff[i][2] * t * t + coeff[i][3] * t * t * t
 */
template<typename SCALAR>
void CoefficientsOfOpenUniformCubicBSpline(
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
  CoefficientsOfOpenUniformCubicBSpline(coeff, idx_segment, num_segment);

  SCALAR v0 = coeff[0][0] + coeff[0][1] * t + coeff[0][2] * t * t + coeff[0][3] * t * t * t;
  SCALAR v1 = coeff[1][0] + coeff[1][1] * t + coeff[1][2] * t * t + coeff[1][3] * t * t * t;
  SCALAR v2 = coeff[2][0] + coeff[2][1] * t + coeff[2][2] * t * t + coeff[2][3] * t * t * t;
  SCALAR v3 = coeff[3][0] + coeff[3][1] * t + coeff[3][2] * t * t + coeff[3][3] * t * t * t;

  assert(fabs(v0 + v1 + v2 + v3 - 1.) < 1.0e-10);
  assert(v0 >= -1.0e-10 && v1 >= -1.0e-10 && v2 >= -1.0e-10 && v3 >= -1.0e-10);
  return poly[idx_segment] * v0 + poly[idx_segment + 1] * v1 + poly[idx_segment + 2] * v2 + poly[idx_segment + 3] * v3;
}

template<typename VEC>
void Polyline_CubicBSplineCurve(
    std::vector<VEC>& sample,
    unsigned int num_point,
    const std::vector<VEC>& cps){
  sample.clear();
  for(unsigned int ip=0;ip<num_point;++ip) {
    double t = static_cast<double>(ip) / static_cast<double>(num_point-1) * (cps.size() - 3);
    VEC p = delfem2::Sample_CubicBsplineCurve(t, cps);
    sample.push_back(p);
  }
}

}

#endif /* DFM2_CURVE_CUBIC_BEZIER_H */
