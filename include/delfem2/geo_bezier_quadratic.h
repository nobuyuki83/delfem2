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

#include "delfem2/dfm2_inline.h"

namespace delfem2::geo_bezier_quadratic {

constexpr static unsigned int NIntLineGauss[4] = {
  1, 2, 3, 4
};

template<typename T>
constexpr static T LineGauss[4][4][2] =
  {
    {
      {0.0, 2.0},
      {0.0, 0.0},
      {0.0, 0.0},
      {0.0, 0.0},
    },
    {
      {-0.577350269189626, 1.0},
      {0.577350269189626, 1.0},
      {0.0, 0.0},
      {0.0, 0.0},
    },
    {
      {-0.774596669241483, 0.555555555555556},
      {0.0, 0.888888888888889},
      {0.774596669241483, 0.555555555555556},
      {0.0, 0.0},
    },
    {
      {-0.861136311594053, 0.347854845137454},
      {-0.339981043584856, 0.652145154862546},
      {0.339981043584856, 0.652145154862546},
      {0.861136311594053, 0.347854845137454},
    }
  };

} // delfem2::pgeo

namespace delfem2 {

template<typename VEC>
VEC PointOnQuadraticBezierCurve(
  double t,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3) {
  double tp = 1.0 - t;
  return (t * t) * p3 + (2 * t * tp) * p2 + (tp * tp) * p1;
}

template<typename VEC>
double Length_QuadraticBezierCurve_Quadrature(
  const VEC &p0,  // end point
  const VEC &p1,  // the control point
  const VEC &p2,  // another end point
  int gauss_order) // order of Gaussian quadrature to use
{
  namespace lcl = ::delfem2::geo_bezier_quadratic;
  using SCALAR = typename VEC::Scalar;

  assert(gauss_order < 4);
  SCALAR totalLength = 0;
  for (unsigned int i = 0; i < lcl::NIntLineGauss[gauss_order]; i++) {
    const double t = (lcl::LineGauss<double>[gauss_order][i][0] + 1) / 2;
    const double w = lcl::LineGauss<double>[gauss_order][i][1];
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
    if( dist0 > dist_min ){ t = t0; break; }   // winding back
    dist_min = dist0;
  }
  return t;
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

}

#endif /* DFM2_CURVE_QUADRATIC_BEZIER_H */
