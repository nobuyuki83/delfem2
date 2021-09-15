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
#include <cassert>

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

template <typename VEC>
std::vector<double> ComputeMeanValueCoordinate(
    const std::vector<double> &vec_xyz_ini,
    const std::vector<double> &vec_xyz_cage0,
    const std::vector<unsigned int> &aTri_cage0) {
  const auto num_point = static_cast<unsigned int>(vec_xyz_ini.size() / 3);
  const auto num_point_cage0 = static_cast<unsigned int>(vec_xyz_cage0.size() / 3);
  std::vector<double> matrix0(num_point * num_point_cage0, 0.0);
  for (unsigned int iq = 0; iq < num_point; ++iq) {
    VEC q0(vec_xyz_ini.data() + iq * 3);
    for (unsigned int itc = 0; itc < aTri_cage0.size() / 3; ++itc) {
      const unsigned int ip0 = aTri_cage0[itc * 3 + 0];
      const unsigned int ip1 = aTri_cage0[itc * 3 + 1];
      const unsigned int ip2 = aTri_cage0[itc * 3 + 2];
      VEC p0 = VEC(vec_xyz_cage0.data() + ip0 * 3) - q0;
      VEC p1 = VEC(vec_xyz_cage0.data() + ip1 * 3) - q0;
      VEC p2 = VEC(vec_xyz_cage0.data() + ip2 * 3) - q0;
      double w[3];
      MeanValueCoordinate_Triangle<VEC>(w, p0, p1, p2);
      matrix0[iq * num_point_cage0 + ip0] += w[0];
      matrix0[iq * num_point_cage0 + ip1] += w[1];
      matrix0[iq * num_point_cage0 + ip2] += w[2];
    }
  }
  return matrix0;
}

/**
 * [np,ndof] = [np,np_cage] * [ndof, np_cage]^T
 * @tparam VEC delfem2::CVec3d
 * @param[in] vec_xyz_ini
 * @param[in] vec_xyz_cage0
 * @param[in] tri_vtxidx_cage0
 * @param[in] ndof
 * @param[in] dof_dofcage
 * @return
 */
template <typename VEC>
std::vector<double> ComputeMeanValueCoordinateReduced(
    const std::vector<double> &vec_xyz,
    const std::vector<double> &vec_xyz_cage0,
    const std::vector<unsigned int> &tri_vtxidx_cage0,
    unsigned int ndof,
    const std::vector<double>& dof_dofcage) {
  const size_t num_vtx = vec_xyz.size() / 3;
  const size_t num_vtx_cage = vec_xyz_cage0.size() / 3;
  const size_t num_tri_cage = tri_vtxidx_cage0.size() / 3;
  assert( dof_dofcage.size() == ndof*num_vtx_cage );
  std::vector<double> mat_res(num_vtx * ndof, 0.0);
  std::vector<double> mat_tmp;
  for (unsigned int iq = 0; iq < num_vtx; ++iq) {
    mat_tmp.assign(num_vtx_cage,0.0);
    const VEC q0(vec_xyz.data() + iq * 3);
    for (unsigned int itc = 0; itc < num_tri_cage; ++itc) {
      const unsigned int ip0 = tri_vtxidx_cage0[itc * 3 + 0];
      const unsigned int ip1 = tri_vtxidx_cage0[itc * 3 + 1];
      const unsigned int ip2 = tri_vtxidx_cage0[itc * 3 + 2];
      const VEC p0 = VEC(vec_xyz_cage0.data() + ip0 * 3) - q0;
      const VEC p1 = VEC(vec_xyz_cage0.data() + ip1 * 3) - q0;
      const VEC p2 = VEC(vec_xyz_cage0.data() + ip2 * 3) - q0;
      double w[3];
      MeanValueCoordinate_Triangle<VEC>(w, p0, p1, p2);
      mat_tmp[ip0] += w[0];
      mat_tmp[ip1] += w[1];
      mat_tmp[ip2] += w[2];
    }
    for(unsigned int idof=0;idof<ndof;++idof) {
      for(unsigned int ip=0;ip<num_vtx_cage;++ip) {
        mat_res[iq * ndof + idof] += dof_dofcage[idof * num_vtx_cage + ip] * mat_tmp[ip];
      }
    }
  }
  return mat_res;
}

// ---------------------------------------

template<class VEC>
void MeanValueCoordinate_Polygon2(
    double *aW,
    double px,
    double py,
    const double *aXY,
    size_t nv) {
  for (unsigned int iv = 0; iv < nv; ++iv) { aW[iv] = 0.0; }
  for (unsigned int iv = 0; iv < nv; ++iv) {
    VEC v0(aXY[iv * 2 + 0] - px, aXY[iv * 2 + 1] - py);
    if (v0.norm() > 1.0e-10) { continue; }
    aW[iv] = 1.0;
    return;
  }
  for (unsigned int ie = 0; ie < nv; ++ie) {
    unsigned int iv0 = (ie + 0) % nv;
    unsigned int iv1 = (ie + 1) % nv;
    VEC v0(aXY[iv0 * 2 + 0] - px, aXY[iv0 * 2 + 1] - py);
    VEC v1(aXY[iv1 * 2 + 0] - px, aXY[iv1 * 2 + 1] - py);
    const double l0 = v0.norm();
    const double l1 = v1.norm();
    if (fabs((v0.dot(v1)) / (l0 * l1) + 1) > 1.0e-10) { continue; }
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
    double c01 = (v0.dot(v1)) / (v0.norm() * v1.norm());
    double s01 = (Cross(v0, v1) > 0) ? 1 : -1;
    double c12 = (v1.dot(v2)) / (v1.norm() * v2.norm());
    double s12 = (Cross(v1, v2) > 0) ? 1 : -1;
    double t01 = s01 * sqrt((1 - c01) / (1 + c01));
    double t12 = s12 * sqrt((1 - c12) / (1 + c12));
    double w1 = (t01 + t12) / v1.norm();
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
