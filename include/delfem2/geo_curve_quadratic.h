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

#ifndef DFM2_CURVE_QUADRATIC_H
#define DFM2_CURVE_QUADRATIC_H

#include <cassert>
#include <cstdio>
#include <vector>
#include <algorithm>

#include "quadrature.h"
#include "polynomial_root.h"
#include "geo_aabb.h"

namespace delfem2 {

/**
 * @param ngauss number of Gausssian quadrature points minus one
 *
 * at parameter "t" the position is evaluated as \n
 * return q0 + q1 * t + q2 * t * t;
 */
template<typename VEC>
double Length_QuadraticCurve_Gauss(
  [[maybe_unused]] const VEC &q0,
  const VEC &q1,
  const VEC &q2,
  int ngauss) {
  const double coe[3] = {q1.squaredNorm(), 4 * q1.dot(q2), 4 * q2.squaredNorm()};

  const unsigned int nw0 = kNumIntegrationPoint_GaussianQuadrature[ngauss];
  const unsigned int nw1 = kNumIntegrationPoint_GaussianQuadrature[ngauss + 1];

  double sum = 0.;
  for (unsigned int i = nw0; i < nw1; i++) {
    const double t = (1. + kPositionWeight_GaussianQuadrature<double>[i][0]) / 2;
    const double w = kPositionWeight_GaussianQuadrature<double>[i][1] / 2;
    const double v = std::sqrt(coe[0] + coe[1] * t + coe[2] * t * t);
    sum += w * v;
  }
  return sum;
}

/**
 * @details at parameter "t" the position is evaluated as \n
 * return q0 + q1 * t + q2 * t * t;
 * https://math.stackexchange.com/questions/12186/arc-length-of-b%c3%a9zier-curves
 */
template<typename VEC>
double Length_QuadraticCurve(
  const VEC &q0,
  const VEC &q1,
  const VEC &q2) {
  if (q2.squaredNorm() < 1.0e-5) {
    return Length_QuadraticCurve_Gauss(q0, q1, q2, 3);
  }
  // dl = sqrt(coe[0] + coe[1]*t + coe[2]*t^2) dt
  const double coe[3] = {
    q1.squaredNorm(),
    4 * q1.dot(q2),
    4 * q2.squaredNorm()};
  auto intt = [&coe](double t) -> double {
    const double tmp0 = (coe[1] + 2 * coe[2] * t) * std::sqrt(coe[0] + t * (coe[1] + coe[2] * t));
    const double tmp3 = sqrt(coe[2]);
    const double tmp1 = coe[1] + 2 * coe[2] * t + 2 * tmp3 * sqrt(coe[0] + t * (coe[1] + coe[2] * t));
    const double tmp2 = (coe[1] * coe[1] - 4 * coe[0] * coe[2]) * std::log(tmp1);
    return tmp0 / (4. * coe[2]) - tmp2 / (8. * coe[2] * tmp3);
  };
  return intt(1) - intt(0);
}


// c + b*t + a*t^2
template<typename VEC>
double Nearest_Origin_QuadraticCurve(
  const VEC &c,
  const VEC &b,
  const VEC &a,
  int num_bisection) {
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

/**
 * @tparam ndim dimension of the geometry
 * @return min_x, min_y, (min_z), max_x, max_y, (max_z)
 */
template<typename VEC>
auto AABB_QuadraticCurve(
  const VEC &c,
  const VEC &b,
  const VEC &a) -> std::array<typename VEC::Scalar, VEC::SizeAtCompileTime * 2> {
  using SCALAR = typename VEC::Scalar;
  constexpr int ndim = VEC::SizeAtCompileTime;
  const VEC p2 = a + b + c;
  std::array<SCALAR, ndim * 2> res;
  for (int i = 0; i < ndim; i++) {
    SCALAR t = b[i] / (-2 * a[i]);
    t = std::clamp<SCALAR>(t, 0, 1);
    SCALAR extreme = c[i] + t * (b[i] + (t * a[i]));
    res[i] = std::min({c[i], p2[i], extreme});
    res[i + ndim] = std::max({c[i], p2[i], extreme});
  }
  return res;
}


// above: quadratic curve
// ===============================================
// below: quadratic Bezier curve

template<typename VEC>
VEC PointOnQuadraticBezierCurve(
  typename VEC::Scalar t,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3) {
  auto tp = 1 - t;
  return (t * t) * p3 + (2 * t * tp) * p2 + (tp * tp) * p1;
}

template<class VEC>
void Polyline_BezierQuadratic(
  std::vector<VEC> &aP,
  const unsigned int n,
  const VEC &p1,
  const VEC &p2,
  const VEC &p3) {
  using SCALAR = typename VEC::Scalar;
  aP.resize(n);
  for (unsigned int i = 0; i < n; ++i) {
    const double t = static_cast<SCALAR>(i) / (static_cast<SCALAR>(n) - 1);
    aP[i] = PointOnQuadraticBezierCurve(
      t,
      p1, p2, p3);
  }
}

template<typename VEC>
double Length_QuadraticBezierCurve_Quadrature(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2,
  int ngauss) {
  const VEC q2 = p0 - 2 * p1 + p2;
  const VEC q1 = -2 * p0 + 2 * p1;
  const VEC q0 = p0;
  return Length_QuadraticCurve_Gauss(q0, q1, q2, ngauss);
}

template<typename VEC>
typename VEC::Scalar Length_QuadraticBezierCurve_Analytic(
  const VEC &p0,  // end point
  const VEC &p1,  // the control point
  const VEC &p2)   // another end point
{
  const VEC q2 = p0 - 2 * p1 + p2;
  const VEC q1 = -2 * p0 + 2 * p1;
  const VEC q0 = p0;
  return Length_QuadraticCurve(q0, q1, q2);
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

  //  SCALAR disttance = (at^2 + bt + c).norm();
  const VEC a = p0 - 2 * p1 + p2;
  const VEC b = -2 * p0 + 2 * p1;
  const VEC c = p0 - q;

  SCALAR t = 0, dist_min = c.norm();  // check t=0
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
typename VEC::Scalar Nearest_QuadraticBezierCurve_Sturm(
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
  return Nearest_Origin_QuadraticCurve(c, b, a, num_bisection);
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
template<typename VEC>
auto AABB_QuadraticBezierCurve(
  const VEC &p0,
  const VEC &p1,
  const VEC &p2) -> std::array<typename VEC::Scalar, VEC::SizeAtCompileTime * 2> {
  // bezier expression: p(t) = at^2 + bt + c
  const VEC a = p0 - 2 * p1 + p2;
  const VEC b = -2 * p0 + 2 * p1;
  const VEC c = p0;
  return AABB_QuadraticCurve(c, b, a);
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

template<typename VEC>
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

  unsigned int np_gauss = 7;
  unsigned int np_kronrod = 15;
  unsigned int noff = kNumIntegrationPoint_GaussianQuadrature[np_gauss - 1];

  // integrate using gaussian and gauss-kronrod quadrature at the same time
  SCALAR length_gauss = 0.;
  SCALAR length_kronrod = 0.;
  for (unsigned int iw = 0; iw < np_gauss; ++iw) { // shared part
    const SCALAR x0 = kPositionWeight_GaussianQuadrature<SCALAR>[noff + iw][0];
    const SCALAR w0 = kPositionWeight_GaussianQuadrature<SCALAR>[noff + iw][1];
    const SCALAR w1 = kPositionWeight_GaussKronrodQuadrature<SCALAR>[iw][1];
    assert(fabs(x0 - kPositionWeight_GaussKronrodQuadrature<SCALAR>[iw][0]) < 1.0e-8);
    const SCALAR nodeValue = quadraticBezierdt((x0 + 1) / 2).norm();
    length_gauss += nodeValue * w0 / 2;
    length_kronrod += nodeValue * w1 / 2;
  }

  // kronrod-only terms
  for (unsigned int iw = np_gauss; iw < np_kronrod; ++iw) {
    const SCALAR x0 = kPositionWeight_GaussKronrodQuadrature<SCALAR>[iw][0];
    const SCALAR w1 = kPositionWeight_GaussKronrodQuadrature<SCALAR>[iw][1];
    const SCALAR nodeValue = quadraticBezierdt((x0 + 1) / 2).norm();
    length_kronrod += nodeValue * w1 / 2;
  }

  if (std::abs(length_kronrod - length_gauss) < std::abs(length_kronrod) * tolerance) {
    return length_kronrod;
  } else {
    if (maxDepth == 1) {
      return length_kronrod;
    } else { // split
      std::array<VEC, 6> subdiv = Split_QuadraticBezierCurve(p0, p1, p2, 0.5);
      const SCALAR len0 = Length_QuadraticBezierCurve_QuadratureSubdiv<VEC>(
        subdiv[0], subdiv[1], subdiv[2], tolerance, maxDepth - 1);
      const SCALAR len1 = Length_QuadraticBezierCurve_QuadratureSubdiv<VEC>(
        subdiv[3], subdiv[4], subdiv[5], tolerance, maxDepth - 1);
      return len0 + len1;
    }
  }

  return -1; // suppress compiler warning
}

// above: quadratic Bezier
// ================================================
// below: quadratic BSpline

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
  unsigned int num_segment) {

  const int k0 = std::clamp<int>(idx_segment - 1, 0, static_cast<int>(num_segment)) - idx_segment;
  const int k1 = std::clamp<int>(idx_segment + 0, 0, static_cast<int>(num_segment)) - idx_segment;
  const int k2 = std::clamp<int>(idx_segment + 1, 0, static_cast<int>(num_segment)) - idx_segment;
  const int k3 = std::clamp<int>(idx_segment + 2, 0, static_cast<int>(num_segment)) - idx_segment;
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

  const double s = t - idx_segment;  // parameter in this segment
  assert(s >= 0 && s <= 1);

  SCALAR coeff[3][3];
  CoefficientsOfOpenUniformBSpline_Quadratic(coeff, idx_segment, num_segment);

  const SCALAR w0 = coeff[0][0] + coeff[0][1] * s + coeff[0][2] * s * s;
  const SCALAR w1 = coeff[1][0] + coeff[1][1] * s + coeff[1][2] * s * s;
  const SCALAR w2 = coeff[2][0] + coeff[2][1] * s + coeff[2][2] * s * s;

  assert(fabs(w0 + w1 + w2 - 1.) < 1.0e-5);
  assert(w0 >= 0 && w1 >= 0 && w2 >= 0);
  return poly[idx_segment] * w0 + poly[idx_segment + 1] * w1 + poly[idx_segment + 2] * w2;
}

/**
 * Quadratic B-Spline with "open and uniform knot vector"
 * knot vector = [0,0,0,1,2,3,...,N-1,N,N,N] where N is the nubmer of segments, which is poly.size()-2
 * @param[in] t parameter of curve that takes [0,N]
 * @param[in] poly position of the the control points
 * @return sampled point
 */
template<typename VEC>
VEC Tangent_QuadraticBsplineCurve(
  typename VEC::Scalar t,
  const std::vector<VEC> &poly) {
  using SCALAR = typename VEC::Scalar;

  const int num_segment = poly.size() - 2;
  const int idx_segment = static_cast<int>(t) + (t == num_segment ? -1 : 0);
  assert(idx_segment >= 0 && idx_segment < num_segment);
  const double s = t - idx_segment;  // parameter in this segment
  assert(s >= 0 && s <= 1);

  SCALAR coeff[3][3];
  CoefficientsOfOpenUniformBSpline_Quadratic(coeff, idx_segment, num_segment);

  const SCALAR w0 = coeff[0][1] + 2 * coeff[0][2] * s;
  const SCALAR w1 = coeff[1][1] + 2 * coeff[1][2] * s;
  const SCALAR w2 = coeff[2][1] + 2 * coeff[2][2] * s;
  return poly[idx_segment] * w0 + poly[idx_segment + 1] * w1 + poly[idx_segment + 2] * w2;
}

template<class VEC>
typename VEC::Scalar Nearest_QuadraticBSplineCurve(
  const std::vector<VEC> &poly,
  const VEC &scr) {
  assert( poly.size() > 2 );
  using SCALAR = typename VEC::Scalar;
  const unsigned int num_segment = poly.size() - 2;
  SCALAR dist_best = (poly[0] - scr).norm();
  SCALAR t_best = 0.;
  for (unsigned int iseg = 0; iseg < num_segment; ++iseg) {
    const VEC p0 = poly[iseg] - scr;
    const VEC p1 = poly[iseg+1] - scr;
    const VEC p2 = poly[iseg+2] - scr;
    if( !IsContact_Orgin_AabbOfPoint3(p0,p1,p2,dist_best) ){ continue; }
    SCALAR coeff[3][3];
    CoefficientsOfOpenUniformBSpline_Quadratic(coeff, iseg, num_segment);
    const VEC q0 = coeff[0][0] * p0 + coeff[1][0] * p1 + coeff[2][0] * p2;
    const VEC q1 = coeff[0][1] * p0 + coeff[1][1] * p1 + coeff[2][1] * p2;
    const VEC q2 = coeff[0][2] * p0 + coeff[1][2] * p1 + coeff[2][2] * p2;
    double t0 = Nearest_Origin_QuadraticCurve(q0,q1,q2, 16);
    double dist0 = (q0 + q1*t0 + q2*t0*t0).norm();
    if( dist0 > dist_best ){ continue; }
    dist_best = dist0;
    t_best = iseg + t0;
    // std::cout << "  " << iseg << " " << num_segment << "  ---> " << t_best << " " << dist_best << std::endl;
  }
  return t_best;
}

/**
 * curve is represented with parameter t as p0 + p1 * t + p2 * t^2
 * @return winding number around origin
 */
template <typename VEC>
typename VEC::Scalar WindingNumber_QuadraticCurve2(
    const VEC& p0,
    const VEC& p1,
    const VEC& p2,
    unsigned int num_bisection,
    unsigned int num_sample){
  using SCALAR = typename VEC::Scalar;
  const SCALAR coe[4] = {
      p1.dot(p0),
      p1.squaredNorm() + 2 * p2.dot(p0),
      3 * p2.dot(p1),
      2 * p2.squaredNorm()};
  std::vector<SCALAR> roots = delfem2::RootsOfPolynomial<4>(coe, num_bisection);
  for(unsigned int i=0;i<num_sample+2;i++){
    roots.push_back( static_cast<SCALAR>(i)/static_cast<SCALAR>(num_sample+1) );
  }
  std::sort(roots.begin(), roots.end());
  SCALAR theta = 0;
  VEC q0 = p0;
  for(unsigned int ir0=0;ir0<roots.size()-1;++ir0) {
    SCALAR t1 = roots[ir0+1];
    const VEC q1 = p2 * t1 * t1 + p1 * t1 + p0;
    const SCALAR sn = q1[1]*q0[0] - q1[0]*q0[1];
    const SCALAR cs = q0[0] * q1[0] + q0[1] * q1[1];
    theta += std::atan2(sn,cs);
    q0 = q1;
  }
  return theta;
}

}

#endif /* DFM2_CURVE_QUADRATIC_H */
