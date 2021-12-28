/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_CURVE_QUADRATIC_BEZIER_H
#define DFM2_CURVE_QUADRATIC_BEZIER_H

#include <cassert>
#include <cstdio>
#include <vector>
#include <algorithm>

#include "quadrature.h"
#include "polynomial_root.h"

namespace delfem2 {

template<typename VEC>
VEC PointOnQuadraticBezierCurve(
  typename VEC::Scalar t,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3) {
  auto tp = 1 - t;
  return (t * t) * p3 + (2 * t * tp) * p2 + (tp * tp) * p1;
}

template<typename VEC>
double Length_QuadraticBezierCurve_Quadrature(
  const VEC &p0,  // end point
  const VEC &p1,  // the control point
  const VEC &p2,  // another end point
  int gauss_order) // order of Gaussian quadrature to use
{
  using SCALAR = typename VEC::Scalar;
  assert(gauss_order < 4);
  SCALAR totalLength = 0;
  for (unsigned int i = 0; i < NIntLineGauss[gauss_order]; i++) {
    const SCALAR t = (LineGauss<SCALAR>[gauss_order][i][0] + 1) / 2;
    const SCALAR w = LineGauss<SCALAR>[gauss_order][i][1];
    const VEC dt = 2 * (1 - t) * (p1 - p0) + 2 * t * (p2 - p1);
    totalLength += dt.norm() * w;
  }
  return totalLength / 2;
}

template<typename VEC>
typename VEC::Scalar Length_QuadraticBezierCurve_Analytic(
  const VEC &p0,  // end point
  const VEC &p1,  // the control point
  const VEC &p2)   // another end point
{
  using SCALAR = typename VEC::Scalar;

  // closed form solution from https://math.stackexchange.com/questions/12186/arc-length-of-b%c3%a9zier-curves
  const SCALAR a = (p2 - 2 * p1 + p0).squaredNorm();

  if (a < 1.0e-5) {
    return Length_QuadraticBezierCurve_Quadrature<VEC>(p0, p1, p2, 3);
  }

  const SCALAR b = (p1 - p0).dot(p2 - p1 - p1 + p0);
  const SCALAR c = (p1 - p0).squaredNorm();
  const SCALAR tmp1 = std::sqrt(a * (c + 2 * b + a));
  const SCALAR tmp2 = std::sqrt(a * c);
  const SCALAR tmp3 = (-b * b + a * c);
  const SCALAR t1 = tmp1 * (b + a) + tmp3 * std::log(b + a + tmp1);
  const SCALAR t0 = tmp2 * b + tmp3 * std::log(b + tmp2);
  return (t1 - t0) / (a * std::sqrt(a));
}

/**
 * @return the parameter of quadratic bezier curve nearest to q
 */
template<typename VEC>
typename VEC::Scalar Nearest_QuadraticBezierCurve(
  const VEC &q,
  const VEC &p0,  // end point
  const VEC &p1,  // the control point
  const VEC &p2,
  unsigned int num_samples,
  unsigned int num_newton_itr) {   // another end point

  using SCALAR = typename VEC::Scalar;
  const VEC a = p0 - 2 * p1 + p2;
  const VEC b = -2 * p0 + 2 * p1;
  const VEC c = p0 - q;

  SCALAR t = 0, dist_min = c.norm();
  for (unsigned int i = 1; i < num_samples + 1; i++) {
    typename VEC::Scalar t0 = static_cast<SCALAR>(i) / static_cast<SCALAR>(num_samples);
    SCALAR dist0 = (a * (t0 * t0) + b * t0 + c).norm();
    if (dist0 < dist_min) {
      dist_min = dist0;
      t = t0;
    }
  }

  const SCALAR s0 = 2 * b.dot(c);
  const SCALAR s1 = 2 * b.squaredNorm() + 4 * a.dot(c);
  const SCALAR s2 = 6 * a.dot(b);
  const SCALAR s3 = 4 * a.squaredNorm();
  const SCALAR u0 = s1;
  const SCALAR u1 = 2 * s2;
  const SCALAR u2 = 3 * s3;

  for (unsigned int itr = 0; itr < num_newton_itr; ++itr) {
    SCALAR t0 = t;
    SCALAR dw = s0 + t * (s1 + t * (s2 + t * s3));
    SCALAR ddw = u0 + t * (u1 + t * u2);
    t -= dw / ddw;
    t = (t < 0) ? 0 : t;
    t = (t > 1) ? 1 : t;
    SCALAR dist0 = (a * (t * t) + b * t + c).norm();
    if (dist0 > dist_min) {
      t = t0;
      break;
    }   // winding back
    dist_min = dist0;
  }
  return t;
}

/**
 * @return the parameter of quadratic bezier curve nearest to q
 */
template<typename VEC>
typename VEC::Scalar Nearest_QuadraticBezierCurve_Strum(
  const VEC &q,
  const VEC &p0,  // end point
  const VEC &p1,  // the control point
  const VEC &p2,
  unsigned int num_bisection) {
  // Precompute coefficients
  // p = at^2 + bt + c
  const VEC a = p0 - 2 * p1 + p2;
  const VEC b = -2 * p0 + 2 * p1;
  const VEC c = p0 - q;

  // Derivative of squared distance function
  // We use the squared distance because it is easier to find its derivative
  const double coe[4] = {
    b.dot(c),
    b.squaredNorm() + 2 * a.dot(c),
    3 * a.dot(b),
    2 * a.squaredNorm()};

  auto roots = RootsOfPolynomial<4>(coe, num_bisection);
  roots.push_back(1); // check t = 1
  double best_t = 0., best_dist = c.squaredNorm(); // check t = 0
  for (double t: roots) {
    double dist0 = (c + t * (b + t * a)).squaredNorm();
    if (dist0 > best_dist) { continue; }
    best_t = t;
    best_dist = dist0;
  }
  return best_t;
}

template<typename VEC>
auto Area_QuadraticBezierCurve2(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2) -> typename VEC::Scalar {
  using T = typename VEC::Scalar;
  const T tmp0 = p1[0] * (p2[1] - p0[1]) + p1[1] * (p0[0] - p2[0]);
  const T tmp1 = p0[0] * p2[1] - p2[0] * p0[1];
  return tmp0 / 3 + tmp1 / 6;
}

/**
 * @tparam ndim dimension of the geometry
 * @return min_x, min_y, (min_z), max_x, max_y, (max_z)
 */
template<int ndim, typename VEC>
auto AABB_QuadraticBezierCurve(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2) -> std::array<typename VEC::Scalar, ndim * 2> {
  using SCALAR = typename VEC::Scalar;

  // bezier expression: p(t) = at^2 + bt + c
  // p'(t) = 2at + b
  const VEC a = p0 - 2 * p1 + p2;
  const VEC b = -2 * p0 + 2 * p1;
  const VEC c = p0;

  std::array<SCALAR, ndim * 2> res;
  for (int i = 0; i < ndim; i++) {
    SCALAR t = b[i] / (-2 * a[i]);
    t = std::clamp<SCALAR>(t, 0, 1);
    SCALAR extreme = c[i] + t * (b[i] + (t * a[i]));
    res[i] = std::min({p0[i], p2[i], extreme});
    res[i + ndim] = std::max({p0[i], p2[i], extreme});
  }

  return res;
}

/**
 * @return 3 control points of the first segment followed by 3 control points of the second segment (in original sequence)
 */
template<typename VEC>
std::array<VEC, 6> Split_QuadraticBezierCurve(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  typename VEC::Scalar t) {
  // p_i^j = p_i^{j-1} * (1 - t0) + p_{i+1}^{j-1} * t0
  auto mix = [&t](const VEC &q0, const VEC &q1) {
    return q0 * (1 - t) + q1 * t;
  };

  std::array<VEC, 6> res;
  res[0] = p0;
  res[1] = mix(p0, p1); // p01
  res[4] = mix(p1, p2); // p11
  res[2] = res[3] = mix(res[1], res[4]); // p02
  res[5] = p2;
  return res;
}

// =====================


template<typename VEC, unsigned gauss_order>
typename VEC::Scalar Length_QuadraticBezierCurve_QuadratureSubdiv(
  const VEC &p0,  // end point
  const VEC &p1,  // the control point next tp p0
  const VEC &p2,  // another end point
  double tolerance = 1E-8,
  int maxDepth = 8) {
  using SCALAR = typename VEC::Scalar;

  // derivative of cubic Bezier
  auto quadraticBezierdt = [&](SCALAR t) {
    return 2 * (1 - t) * (p1 - p0) + 2 * t * (p2 - p1);
  };

  constexpr auto gdw = gauss_detail<double, gauss_order>::weights();
  constexpr auto gda = gauss_detail<double, gauss_order>::abscissa();
  constexpr auto krw = gauss_kronrod_detail<double, gauss_order * 2 + 1>::weights();
  constexpr auto kra = gauss_kronrod_detail<double, gauss_order * 2 + 1>::abscissa();

  // integrate using gaussian and gauss-kronrod quadrature at the same time
  SCALAR length_gauss = 0;
  SCALAR length_gauss_kronrod = 0;

  { // 0
    SCALAR nodeValue = quadraticBezierdt(0.5).norm();
    length_gauss += nodeValue * gdw[0];
    length_gauss_kronrod += nodeValue * krw[0];
  }
  for (int gauss_node = 1; gauss_node < (gauss_order + 1) / 2; gauss_node++) {
    SCALAR a0 = gda[gauss_node];
    {
      SCALAR nodeValue = quadraticBezierdt((1 + a0) / 2).norm();
      length_gauss += nodeValue * gdw[gauss_node];
      length_gauss_kronrod += nodeValue * krw[gauss_node * 2];
    }
    {
      SCALAR nodeValue = quadraticBezierdt((1 - a0) / 2).norm();
      length_gauss += nodeValue * gdw[gauss_node];
      length_gauss_kronrod += nodeValue * krw[gauss_node * 2];
    }
  }

  // kronrod-only terms
  for (int kronrod_node = 1; kronrod_node <= gauss_order; kronrod_node += 2) {
    double a0 = kra[kronrod_node];
    {
      const SCALAR nodeValue = quadraticBezierdt((1 + a0) / 2).norm();
      length_gauss_kronrod += nodeValue * krw[kronrod_node];
    }
    {
      const SCALAR nodeValue = quadraticBezierdt((1 - a0) / 2).norm();
      length_gauss_kronrod += nodeValue * krw[kronrod_node];
    }
  }

  length_gauss /= 2;
  length_gauss_kronrod /= 2;

  if (std::abs(length_gauss_kronrod - length_gauss) < std::abs(length_gauss_kronrod) * tolerance) {
    return length_gauss_kronrod;
  } else {
    if (maxDepth == 1) {
      // std::cout << "Warning: Max depth reached, current estimated error = " << std::abs(length_gauss_kronrod - length_gauss) / std::abs(length_gauss_kronrod) << std::endl;
      return length_gauss_kronrod;
    } else { // split
      std::array<VEC, 6> subdiv = Split_QuadraticBezierCurve(p0, p1, p2, 0.5);
      const SCALAR len0 = Length_QuadraticBezierCurve_QuadratureSubdiv<VEC, gauss_order>(
        subdiv[0], subdiv[1], subdiv[2], tolerance, maxDepth - 1);
      const SCALAR len1 = Length_QuadraticBezierCurve_QuadratureSubdiv<VEC, gauss_order>(
        subdiv[3], subdiv[4], subdiv[5], tolerance, maxDepth - 1);
      return len0 + len1;
    }
  }

  return -1; // suppress compiler warning
}

}

#endif /* DFM2_CURVE_QUADRATIC_BEZIER_H */
