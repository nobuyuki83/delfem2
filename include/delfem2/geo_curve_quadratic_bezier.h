/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_BEZIER_QUADRATIC_H
#define DFM2_BEZIER_QUADRATIC_H

#include <cassert>
#include <cstdio>
#include <vector>

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
  using SCALAR = decltype(p0[0]);
    
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
double Length_QuadraticBezierCurve_Analytic(
  const VEC &p0,  // end point
  const VEC &p1,  // the control point
  const VEC &p2)   // another end point
{
  using SCALAR = decltype(p0[0]);

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
double Nearest_QuadraticBezierCurve(
  const VEC &q,
  const VEC &p0,  // end point
  const VEC &p1,  // the control point
  const VEC &p2,
  unsigned int num_samples,
  unsigned int num_newton_itr) {   // another end point

  // Precompute coefficients
  // p = at^2 + bt + c
  const VEC a = p0 - 2 * p1 + p2;
  const VEC b = -2 * p0 + 2 * p1;
  const VEC c = p0;

  double t = 0, dist_min = (p0 - q).norm();
  for (unsigned int i = 1; i < num_samples + 1; i++) {
    double t0 = static_cast<double>(i) / static_cast<double>(num_samples);
    double dist0 = (a * (t0 * t0) + b * t0 + c - q).norm();
    if (dist0 < dist_min) {
      dist_min = dist0;
      t = t0;
    }
  }

  const double s0 = 2 * b.dot(c - q);
  const double s1 = 2 * b.squaredNorm() + 4 * a.dot(c - q);
  const double s2 = 6 * a.dot(b);
  const double s3 = 4 * a.squaredNorm();
  const double u0 = s1;
  const double u1 = 2*s2;
  const double u2 = 3*s3;

  for (unsigned int itr = 0; itr < num_newton_itr; ++itr) {
    double dw = s0 + t * (s1 + t * (s2 + t * s3));
    double ddw = u0 + t * (u1 + t * u2);
    t -= dw / ddw;
    t = (t < 0) ? 0 : t;
    t = (t > 1) ? 1 : t;
    // t = std::clamp(t, 0., 1.);
  }
  return t;
}

template <typename VEC>
auto Area_QuadraticBezierCurve(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2)  -> decltype(p0[0])
{
  using T = decltype(p0[0]);
  const T tmp0 = p1[0] * (p2[1] - p0[1]) + p1[1] * (p0[0] - p2[0]);
  const T tmp1 = p0[0] * p2[1] - p2[0] * p0[1];
  return tmp0 / 3 +  tmp1 / 6;
}

}

#endif /* DFM2_PGEO_H */
