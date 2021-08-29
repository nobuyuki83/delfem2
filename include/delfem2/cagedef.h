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
#include <vector>

namespace delfem2::cagedef {

template<typename VEC>
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
template<class VEC>
void MeanValueCoordinate_Triangle(
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
  if (std::isnan(s0) || std::isnan(s1) || std::isnan(s2)) {
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

// --------

template<class VEC>
void MeanValueCoordinate_Polygon2(
    double *aW,
    double px, double py,
    const double *aXY,
    unsigned int nv) {
  for (unsigned int iv = 0; iv < nv; ++iv) { aW[iv] = 0.0; }
  for (unsigned int iv = 0; iv < nv; ++iv) {
    VEC v0(aXY[iv * 2 + 0] - px, aXY[iv * 2 + 1] - py);
    if (v0.Length() > 1.0e-10) { continue; }
    aW[iv] = 1.0;
    return;
  }
  for (unsigned int ie = 0; ie < nv; ++ie) {
    unsigned int iv0 = (ie + 0) % nv;
    unsigned int iv1 = (ie + 1) % nv;
    VEC v0(aXY[iv0 * 2 + 0] - px, aXY[iv0 * 2 + 1] - py);
    VEC v1(aXY[iv1 * 2 + 0] - px, aXY[iv1 * 2 + 1] - py);
    const double l0 = v0.Length();
    const double l1 = v1.Length();
    if (fabs((v0 * v1) / (l0 * l1) + 1) > 1.0e-10) { continue; }
    aW[iv0] = l1 / (l0 + l1);
    aW[iv1] = l0 / (l0 + l1);
    return;
  }
  double sum = 0;
  for (unsigned int ie = 0; ie < nv; ++ie) {
    unsigned int iv0 = (ie + 0) % nv;
    unsigned int iv1 = (ie + 1) % nv;
    unsigned int iv2 = (ie + 2) % nv;
    VEC v0(aXY[iv0 * 2 + 0] - px, aXY[iv0 * 2 + 1] - py);
    VEC v1(aXY[iv1 * 2 + 0] - px, aXY[iv1 * 2 + 1] - py);
    VEC v2(aXY[iv2 * 2 + 0] - px, aXY[iv2 * 2 + 1] - py);
    double c01 = (v0 * v1) / (v0.Length() * v1.Length());
    double s01 = (Cross(v0, v1) > 0) ? 1 : -1;
    double c12 = (v1 * v2) / (v1.Length() * v2.Length());
    double s12 = (Cross(v1, v2) > 0) ? 1 : -1;
    double t01 = s01 * sqrt((1 - c01) / (1 + c01));
    double t12 = s12 * sqrt((1 - c12) / (1 + c12));
    double w1 = (t01 + t12) / v1.Length();
    aW[iv1] = w1;
    sum += w1;
  }
  for (unsigned int iv = 0; iv < nv; ++iv) {
    aW[iv] /= sum;
  }
}

// --------------

template<class VEC>
void MeanValueCoordinate_Polygon2(
    std::vector<double> &aW,
    VEC &p,
    std::vector<VEC> &aVtx) {
  const int nv = (int) aVtx.size();
  aW.assign(nv, 0.0);
  double sum = 0;
  for (int ie = 0; ie < nv; ++ie) {
    int iv0 = (ie + 0) % nv;
    int iv1 = (ie + 1) % nv;
    int iv2 = (ie + 2) % nv;
    VEC v0 = aVtx[iv0] - p;
    VEC v1 = aVtx[iv1] - p;
    VEC v2 = aVtx[iv2] - p;
    double c01 = (v0 * v1) / (v0.Length() * v1.Length());
    double c12 = (v1 * v2) / (v1.Length() * v2.Length());
    double t01 = sqrt((1 - c01) / (1 + c01));
    double t12 = sqrt((1 - c12) / (1 + c12));
    double w1 = (t01 + t12) / v1.Length();
    aW[iv1] = w1;
    sum += w1;
  }
  for (int iv = 0; iv < nv; ++iv) {
    aW[iv] /= sum;
  }
}

}

#endif // DFM2_CAGEDEF_H
