//
// Created by Nobuyuki Umetani on 2021/11/25.
//

#ifndef DFM2_SAMPLER_TRIMESH_H_
#define DFM2_SAMPLER_TRIMESH_H_

#include <random>

#include "delfem2/vec3.h"

namespace delfem2 {

class RandomSamplingOnMeshTri3 {
 public:
  RandomSamplingOnMeshTri3(
    const std::vector<double> &vtx_xyz_,
    const std::vector<unsigned int> &tri_vtx_)
    : vtx_xyz(vtx_xyz_), tri_vtx(tri_vtx_) {
    unsigned int ntri = tri_vtx.size() / 3;
    cumulative_area_sum.reserve(ntri + 1);
    cumulative_area_sum.push_back(0.);
    for (unsigned int itri = 0; itri < ntri; ++itri) {
      const double a0 = delfem2::Area_Tri3(
        vtx_xyz.data() + tri_vtx[itri * 3 + 0] * 3,
        vtx_xyz.data() + tri_vtx[itri * 3 + 1] * 3,
        vtx_xyz.data() + tri_vtx[itri * 3 + 2] * 3);
      const double t0 = cumulative_area_sum[cumulative_area_sum.size() - 1];
      cumulative_area_sum.push_back(a0 + t0);
    }
    rndeng.seed(std::random_device{}());
    dist_01 = std::uniform_real_distribution<double>(0, 1);
  }
  std::tuple<unsigned int, double, double> Sample() {
    double val01 = dist_01(rndeng);
    unsigned int ntri = tri_vtx.size() / 3;
    assert(cumulative_area_sum.size() == ntri + 1);
    double a0 = val01 * cumulative_area_sum[ntri];
    unsigned int itri_l = 0, itri_u = ntri;
    for (;;) {
      assert(cumulative_area_sum[itri_l] < a0);
      assert(a0 < cumulative_area_sum[itri_u]);
      unsigned int itri_h = (itri_u + itri_l) / 2;
      if (itri_u - itri_l == 1) { break; }
      if (cumulative_area_sum[itri_h] < a0) { itri_l = itri_h; }
      else { itri_u = itri_h; }
    }
    assert(cumulative_area_sum[itri_l] < a0);
    assert(a0 < cumulative_area_sum[itri_l + 1]);
    double r0 = (a0 - cumulative_area_sum[itri_l]) / (cumulative_area_sum[itri_l + 1] - cumulative_area_sum[itri_l]);
    double r1 = dist_01(rndeng);
    if (r0 + r1 > 1) {
      double r0a = r0, r1a = r1;
      r0 = 1 - r1a;
      r1 = 1 - r0a;
    }
    return {itri_l, r0, r1};
  }
 public:
  const std::vector<double> &vtx_xyz;
  const std::vector<unsigned int> &tri_vtx;
  std::vector<double> cumulative_area_sum;
  std::mt19937 rndeng;
  std::uniform_real_distribution<double> dist_01;
};


class RandomSamplingOnMeshTri3Selective {
 public:
  RandomSamplingOnMeshTri3Selective(
    const std::vector<double> &vtx_xyz_,
    const std::vector<unsigned int> &tri_vtx_,
    const std::function<bool (unsigned int)> &fnc)
    : vtx_xyz(vtx_xyz_), tri_vtx(tri_vtx_) {
    unsigned int ntri = tri_vtx.size() / 3;
    cumulative_area_sum.reserve(ntri + 1);
    cumulative_area_sum.push_back(0.);
    for (unsigned int itri = 0; itri < ntri; ++itri) {
      if( !fnc(itri) ){ continue; }
      const double a0 = delfem2::Area_Tri3(
        vtx_xyz.data() + tri_vtx[itri * 3 + 0] * 3,
        vtx_xyz.data() + tri_vtx[itri * 3 + 1] * 3,
        vtx_xyz.data() + tri_vtx[itri * 3 + 2] * 3);
      const double t0 = cumulative_area_sum[cumulative_area_sum.size() - 1];
      flggedtri_tri.push_back(itri);
      cumulative_area_sum.push_back(a0 + t0);
    }
    rndeng.seed(std::random_device{}());
    dist_01 = std::uniform_real_distribution<double>(0, 1);
  }
  std::tuple<unsigned int, double, double> Sample() {
    double val01 = dist_01(rndeng);
    unsigned int num_flaggedtri = flggedtri_tri.size();
    assert(cumulative_area_sum.size() == flggedtri_tri.size() + 1);
    double a0 = val01 * cumulative_area_sum[num_flaggedtri];
    unsigned int itri_l = 0, itri_u = num_flaggedtri;
    for (;;) {
      assert(cumulative_area_sum[itri_l] < a0);
      assert(a0 < cumulative_area_sum[itri_u]);
      unsigned int itri_h = (itri_u + itri_l) / 2;
      if (itri_u - itri_l == 1) { break; }
      if (cumulative_area_sum[itri_h] < a0) { itri_l = itri_h; }
      else { itri_u = itri_h; }
    }
    assert(cumulative_area_sum[itri_l] < a0);
    assert(a0 < cumulative_area_sum[itri_l + 1]);
    double r0 = (a0 - cumulative_area_sum[itri_l]) / (cumulative_area_sum[itri_l + 1] - cumulative_area_sum[itri_l]);
    double r1 = dist_01(rndeng);
    if (r0 + r1 > 1) {
      double r0a = r0, r1a = r1;
      r0 = 1 - r1a;
      r1 = 1 - r0a;
    }
    return {flggedtri_tri[itri_l], r0, r1};
  }
 public:
  const std::vector<double> &vtx_xyz;
  const std::vector<unsigned int> &tri_vtx;
  std::vector<double> cumulative_area_sum;
  std::vector<unsigned int> flggedtri_tri;
  std::mt19937 rndeng;
  std::uniform_real_distribution<double> dist_01;
};

}


#endif // DFM2_SAMPLER_TRIMESH_H_
