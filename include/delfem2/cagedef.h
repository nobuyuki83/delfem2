//
// Created by Nobuyuki Umetani on 2021-08-29.
//

/**
 * implementation based on
 * - Tao Ju, Scott Schaefer, and Joe Warren. 2005.
 *   "Mean value coordinates for closed triangular meshes".
 *   ACM Trans. Graph. 24, 3 (July 2005), 561–566.
 */

#ifndef DFM2_CAGEDEF_H
#define DFM2_CAGEDEF_H

#include <cmath>

namespace delfem2::cagedef {

template <typename VEC>
double ScalarTripleProduct(
    const VEC &a,
    const VEC &b,
    const VEC &c) {
  return
  a[0] * (b[1] * c[2] - b[2] * c[1]) +
  a[1] * (b[2] * c[0] - b[0] * c[2]) +
  a[2] * (b[0] * c[1] - b[1] * c[0]);
}

}  // namespace delfem2::cagedef

namespace delfem2 {

/**
 * compute weight for the mean value coordinate.
 * @tparam VEC should work for "delfem2::CVec3" or "Eigen::Vector3"
 * @param w
 * @param v0
 * @param v1
 * @param v2
 */
template <class VEC>
void MeanValueCoordinate(
    double w[3],
    const VEC &v0,
    const VEC &v1,
    const VEC &v2) {
  namespace lcl = delfem2::cagedef;
  double eps = 1.0e-5;
  double d0 = v0.norm();
  double d1 = v1.norm();
  double d2 = v2.norm();
  const VEC u0 = v0 / d0;
  const VEC u1 = v1 / d1;
  const VEC u2 = v2 / d2;
  double l0 = (u1 - u2).norm();
  double l1 = (u2 - u0).norm();
  double l2 = (u0 - u1).norm();
  if (l0 < eps || l1 < eps || l2 < eps) {
    w[0] = 0;
    w[1] = 0;
    w[2] = 0;
    return;
  }
  double t0 = 2 * asin(l0 * 0.5);
  double t1 = 2 * asin(l1 * 0.5);
  double t2 = 2 * asin(l2 * 0.5);
  double h = (t0 + t1 + t2) * 0.5;
  double c0 = 2 * sin(h) * sin(h - t0) / (sin(t1) * sin(t2)) - 1;
  double c1 = 2 * sin(h) * sin(h - t1) / (sin(t2) * sin(t0)) - 1;
  double c2 = 2 * sin(h) * sin(h - t2) / (sin(t0) * sin(t1)) - 1;
  double vol012 = ScalarTripleProduct(u0, u1, u2);
  double sign = (vol012 > 0) ? 1 : -1;
  double s0 = sign * sqrt(1.0 - c0 * c0);
  double s1 = sign * sqrt(1.0 - c1 * c1);
  double s2 = sign * sqrt(1.0 - c2 * c2);
  if (isnan(s0) || isnan(s1) || isnan(s2)) {
    w[0] = 0;
    w[1] = 0;
    w[2] = 0;
    return;
  }
  if (fabs(d0 * sin(t1) * s2) < eps || fabs(d1 * sin(t2) * s0) < eps || fabs(d2 * sin(t0) * s1) < eps) {
    w[0] = 0;
    w[1] = 0;
    w[2] = 0;
    return;
  }
  w[0] = (t0 - c2 * t1 - c1 * t2) / (d0 * sin(t1) * s2);
  w[1] = (t1 - c0 * t2 - c2 * t0) / (d1 * sin(t2) * s0);
  w[2] = (t2 - c1 * t0 - c0 * t1) / (d2 * sin(t0) * s1);
}

}

#endif // DFM2_CAGEDEF_H
