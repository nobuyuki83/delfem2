/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_CURVE_CUBIC_BEZIER_H
#define DFM2_CURVE_CUBIC_BEZIER_H

#include <cassert>
#include <cstdio>
#include <vector>

namespace delfem2::geo_bezier_cubic {

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
VEC PointOnCubicBezierCurve(
  double t,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3,
  const VEC &p4) {
  double tp = 1.0 - t;
  return t * t * t * p4
    + 3 * t * t * tp * p3
    + 3 * t * tp * tp * p2
    + tp * tp * tp * p1;
}

template<typename T>
T getTangentCubicBezierCurve(
  double t,
  const T &p1, const T &p2, const T &p3, const T &p4) {
  double tp = 1.0 - t;
  return 3 * t * t * p4
    + 3 * t * (2 - 3 * t) * p3
    + 3 * tp * (1 - 3 * t) * p2
    - 3 * tp * tp * p1;
}

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
double Length_CubicBezierCurve_Quadrature(
  const VEC &p0,  // end point
  const VEC &p1,  // the control point next tp p0
  const VEC &p2,  // the control point next to p1
  const VEC &p3,  // another end point
  int gauss_order)  // order of Gaussian quadrature to use
{
  namespace lcl = ::delfem2::geo_bezier_cubic;
  assert(gauss_order < 4);
  using SCALAR = decltype(p0[0]);
  SCALAR totalLength = 0;
  for (unsigned int i = 0; i < lcl::NIntLineGauss[gauss_order]; i++) {
    double t = (lcl::LineGauss<double>[gauss_order][i][0] + 1) / 2;
    double w = lcl::LineGauss<double>[gauss_order][i][1];
    const VEC dt = 3 * (1 - t) * (1 - t) * (p1 - p0)
      + 6 * (1 - t) * t * (p2 - p1)
      + 3 * t * t * (p3 - p2);
    totalLength += dt.norm() * w;
  }
  return totalLength / 2;
}

template<typename VEC>
double Nearest_CubicBezierCurve(
  const VEC &q,
  const VEC &p0,  // end point
  const VEC &p1,  // the control point
  const VEC &p2,
  const VEC &p3,
  unsigned int num_samples,
  unsigned int num_newton_itr) {   // another end point

  const VEC a = -p0 + 3 * p1 - 3 * p2 + p3;
  const VEC b = 3 * p0 - 6 * p1 + 3 * p2;
  const VEC c = -3 * p0 + 3 * p1;
  const VEC d = p0 - q;

  double t = 0, dist_min = (p0 - q).norm();
  for (unsigned int i = 1; i < num_samples + 1; i++) {
    double t0 = static_cast<double>(i) / static_cast<double>(num_samples);
    double dist0 = (a * (t0 * t0 * t0) + b * (t0 * t0) + c * t0 + d).norm();
    if (dist0 < dist_min) {
      dist_min = dist0;
      t = t0;
    }
  }

  const double s0 = 2 * c.dot(d);
  const double s1 = 2 * c.squaredNorm() + 4 * b.dot(d);
  const double s2 = 6 * b.dot(c) + 6 * a.dot(d);
  const double s3 = 4 * b.squaredNorm() + 8 * a.dot(c);
  const double s4 = 10 * a.dot(b);
  const double s5 = 6 * a.squaredNorm();
  const double u0 = s1;
  const double u1 = 2 * s2;
  const double u2 = 3 * s3;
  const double u3 = 4 * s4;
  const double u4 = 5 * s5;

  for (unsigned int itr = 0; itr < num_newton_itr; ++itr) {
    double dw = s0 + t * (s1 + t * (s2 + t * (s3 + t * (s4 + t * s5))));
    double ddw = u0 + t * (u1 + t * (u2 + t * (u3 + t * u4)));
    t -= dw / ddw;
    t = (t < 0) ? 0 : t;
    t = (t > 1) ? 1 : t;
  }
  return t;
}

template <typename VEC>
auto Area_CubicBezierCurve(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3) -> decltype(p0[0])
{
  using T = decltype(p0[0]);
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
template <int ndim, typename VEC>
auto AABB_CubicBezierCurve(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3) -> std::array<decltype(p0[0]),ndim*2>
{
    using SCALAR = decltype(p0[0]);
    std::array<SCALAR, ndim*2> res;
    /*
     write something here
     */
    return res;
}

}

#endif /* DFM2_PGEO_H */
