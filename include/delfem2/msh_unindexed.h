//
// Created by Nobuyuki Umetani on 2021/11/29.
//

#ifndef DFM2_MSH_UNINDEXED_H_
#define DFM2_MSH_UNINDEXED_H_

#include <vector>

#include "delfem2/vec3.h"

namespace delfem2 {

void UnidexedVertexDataTriMesh(
  std::vector<double> &tri_xyz,
  const std::vector<double> &vtx_xyz,
  const std::vector<unsigned int> &tri_vtx) {
  const unsigned int ntri = tri_vtx.size() / 3;
  tri_xyz.resize(ntri * 9);
  for (unsigned int it = 0; it < tri_vtx.size() / 3; ++it) {
    for (unsigned int ino = 0; ino < 3; ++ino) {
      const double *p0 = vtx_xyz.data() + tri_vtx[it * 3 + ino] * 3;
      tri_xyz[it * 9 + ino * 3 + 0] = p0[0];
      tri_xyz[it * 9 + ino * 3 + 1] = p0[1];
      tri_xyz[it * 9 + ino * 3 + 2] = p0[2];
    }
  }
}

void UnindexedNormalTriMesh3(
  std::vector<double> &tri_nrm,
  const std::vector<double> &vtx_xyz,
  const std::vector<unsigned int> &tri_vtx) {
  const unsigned int ntri = tri_vtx.size() / 3;
  tri_nrm.resize(ntri * 9);
  for (unsigned int it = 0; it < tri_vtx.size() / 3; ++it) {
    double n0[3];
    {
      const double *p0 = vtx_xyz.data() + tri_vtx[it * 3 + 0] * 3;
      const double *p1 = vtx_xyz.data() + tri_vtx[it * 3 + 1] * 3;
      const double *p2 = vtx_xyz.data() + tri_vtx[it * 3 + 2] * 3;
      double area;
      delfem2::UnitNormalAreaTri3(n0, area, p0, p1, p2);
    }
    for (unsigned int ino = 0; ino < 3; ++ino) {
      tri_nrm[it * 9 + ino * 3 + 0] = n0[0];
      tri_nrm[it * 9 + ino * 3 + 1] = n0[1];
      tri_nrm[it * 9 + ino * 3 + 2] = n0[2];
    }
  }
}

void UnindexedColorTriMesh3(
  std::vector<double> &tri_rgb,
  const double flg_rgb[][3],
  const std::vector<int> &tri_flg,
  const std::vector<unsigned int> &tri_vtx) {
  const unsigned int ntri = tri_vtx.size() / 3;
  tri_rgb.resize(ntri * 9);
  for (unsigned int it = 0; it < tri_vtx.size() / 3; ++it) {
    unsigned int iflg = tri_flg[it];
    for (unsigned int ino = 0; ino < 3; ++ino) {
      tri_rgb[it * 9 + ino * 3 + 0] = flg_rgb[iflg][0];
      tri_rgb[it * 9 + ino * 3 + 1] = flg_rgb[iflg][1];
      tri_rgb[it * 9 + ino * 3 + 2] = flg_rgb[iflg][2];
    }
  }
}

void UnindexedEdgeMesh3_UnindexedTrimesh3
 (std::vector<double> &edge_xyz,
  const std::vector<double> &tri_xyz) {
  unsigned int ntri = tri_xyz.size()/9;
  edge_xyz.resize(ntri * 18);
  for (unsigned int itri = 0; itri < ntri; itri++) {
    for (unsigned int i = 0; i < 3; i++) {
      edge_xyz[itri * 18 + 0 + i] = tri_xyz[itri * 9 + 0 + i];
      edge_xyz[itri * 18 + 3 + i] = tri_xyz[itri * 9 + 3 + i];
      edge_xyz[itri * 18 + 6 + i] = tri_xyz[itri * 9 + 3 + i];
      edge_xyz[itri * 18 + 9 + i] = tri_xyz[itri * 9 + 6 + i];
      edge_xyz[itri * 18 + 12 + i] = tri_xyz[itri * 9 + 6 + i];
      edge_xyz[itri * 18 + 15 + i] = tri_xyz[itri * 9 + 0 + i];
    }
  }
}

}

#endif //DFM2_MSH_UNINDEXED_H_
