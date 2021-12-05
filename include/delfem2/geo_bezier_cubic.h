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
#include <algorithm>

namespace delfem2::geo_bezier_cubic {


constexpr static unsigned int NIntLineGauss[6] = {
  1, 2, 3, 4, 5, 6
};

// https://pomax.github.io/bezierinfo/legendre-gauss.html
template<typename T>
constexpr static T LineGauss[6][6][2] =
  {
    { // 0
      {0.0, 2.0},
    },
    { // 1
      {-0.577350269189626, 1.0},
      {0.577350269189626, 1.0},
    },
    {  // 2
      {-0.774596669241483, 0.555555555555556},
      {0.0, 0.888888888888889},
      {0.774596669241483, 0.555555555555556},
    },
    {  // 3
      {-0.861136311594053, 0.347854845137454},
      {-0.339981043584856, 0.652145154862546},
      {0.339981043584856, 0.652145154862546},
      {0.861136311594053, 0.347854845137454},
    },
    {  // 4
      {0.0000000000000000, 0.5688888888888889},
      {-0.5384693101056831, 0.4786286704993665},
      {0.5384693101056831, 0.4786286704993665},
      {-0.9061798459386640, 0.2369268850561891},
      {0.9061798459386640, 0.2369268850561891},
    },
    {  // 5
      {0.6612093864662645, 0.3607615730481386},
      {-0.6612093864662645, 0.3607615730481386},
      {-0.2386191860831969, 0.4679139345726910},
      {0.2386191860831969, 0.4679139345726910},
      {-0.9324695142031521, 0.1713244923791704},
      {0.9324695142031521, 0.1713244923791704}
    }
  };

} // delfem2::pgeo

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
typename VEC::Scalar Length_CubicBezierCurve_Quadrature(
  const VEC &p0,  // end point
  const VEC &p1,  // the control point next tp p0
  const VEC &p2,  // the control point next to p1
  const VEC &p3,  // another end point
  int gauss_order)  // order of Gaussian quadrature to use
{
  namespace lcl = ::delfem2::geo_bezier_cubic;
  assert(gauss_order < 6);
  using SCALAR = typename VEC::Scalar;
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
    if( dist0 > dist_min ){ t = t0; break; }   // winding back
    dist_min = dist0;
  }
  return t;
}

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
template<int ndim, typename VEC>
auto AABB_CubicBezierCurve(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3) -> std::array<typename VEC::Scalar, ndim * 2> {
  using SCALAR = typename VEC::Scalar;

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
template <typename VEC>
std::array<VEC,8> Split_CubicBezierCurve(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3,
  typename VEC::Scalar t)
{
    // p_i^j = p_i^{j-1} * (1 - t0) + p_{i+1}^{j-1} * t0
    auto mix = [&t](const VEC &q0, const VEC &q1) {
        return q0 * (1 - t) + q1 * t;
    };

    std::array<VEC,8> res;
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


}

#endif /* DFM2_CURVE_CUBIC_BEZIER_H */
