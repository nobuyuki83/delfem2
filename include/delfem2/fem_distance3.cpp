/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/fem_distance3.h"

#include "delfem2/vec2.h"
#include "delfem2/mat3.h"
#include "delfem2/mat3_funcs.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/vec3_funcs.h"

// -------------------------------------------------

DFM2_INLINE double delfem2::WdWddW_SquareLengthLineseg3D(
    CVec3d dw_dp[2],
    CMat3d ddw_ddp[2][2],
    double stiffness,
    const CVec3d vtx_xyz_ini[2],
    double edge_length_ini) {
  const double l = sqrt(
      (vtx_xyz_ini[0][0] - vtx_xyz_ini[1][0]) * (vtx_xyz_ini[0][0] - vtx_xyz_ini[1][0]) +
          (vtx_xyz_ini[0][1] - vtx_xyz_ini[1][1]) * (vtx_xyz_ini[0][1] - vtx_xyz_ini[1][1]) +
          (vtx_xyz_ini[0][2] - vtx_xyz_ini[1][2]) * (vtx_xyz_ini[0][2] - vtx_xyz_ini[1][2]));
  const CVec3d v = vtx_xyz_ini[0] - vtx_xyz_ini[1];
  const double R = edge_length_ini - l;
  dw_dp[0] = (-R * stiffness / l) * v;
  dw_dp[1] = (+R * stiffness / l) * v;
  const CMat3d mvv = Mat3_OuterProduct(v, v);
  double t0 = stiffness * edge_length_ini / (l * l * l);
  double t1 = stiffness * (l - edge_length_ini) / l;
  const CMat3d m = t0 * mvv + t1 * CMat3d::Identity(1.0);
  ddw_ddp[0][0] = +m;
  ddw_ddp[0][1] = -m;
  ddw_ddp[1][0] = -m;
  ddw_ddp[1][1] = +m;
  return 0.5 * stiffness * R * R;
}

template<typename T>
DFM2_INLINE void delfem2::CdC_SquareLengthLineseg3D(
    T &c,
    T dc_dpos[2][3],
    const T pos_xyz[2][3],
    T length_ini) {
  const T vx = pos_xyz[0][0] - pos_xyz[1][0];
  const T vy = pos_xyz[0][1] - pos_xyz[1][1];
  const T vz = pos_xyz[0][2] - pos_xyz[1][2];
  const T l = std::sqrt(vx * vx + vy * vy + vz * vz);
  c = l - length_ini;
  dc_dpos[0][0] = vx / l;
  dc_dpos[0][1] = vy / l;
  dc_dpos[0][2] = vz / l;
  dc_dpos[1][0] = -dc_dpos[0][0];
  dc_dpos[1][1] = -dc_dpos[0][1];
  dc_dpos[1][2] = -dc_dpos[0][2];
}



DFM2_INLINE void delfem2::CdC_LengthTriEdges23(
  double C[3],
  double dCdp[3][9],
  const double P[3][2], // (in) undeformed triangle vertex positions
  const double p[3][3] // (in) deformed triangle vertex positions
) {
  const double L12 = Distance2(P[1], P[2]);
  const double L20 = Distance2(P[2], P[0]);
  const double L01 = Distance2(P[0], P[1]);
  CVec3d v12(p[1][0] - p[2][0], p[1][1] - p[2][1], p[1][2] - p[2][2]);
  CVec3d v20(p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]);
  CVec3d v01(p[0][0] - p[1][0], p[0][1] - p[1][1], p[0][2] - p[1][2]);
  const double l12 = v12.norm();
  const double l20 = v20.norm();
  const double l01 = v01.norm();
  C[0] = l12 - L12;
  C[1] = l20 - L20;
  C[2] = l01 - L01;
  v12 /= l12;
  v20 /= l20;
  v01 /= l01;
  for (int i = 0; i < 27; ++i) { (&dCdp[0][0])[i] = 0.0; }
  v12.CopyTo(dCdp[0] + 3 * 1);
  v12.CopyToScale(dCdp[0] + 3 * 2, -1.0);
  v20.CopyTo(dCdp[1] + 3 * 2);
  v20.CopyToScale(dCdp[1] + 3 * 0, -1.0);
  v01.CopyTo(dCdp[2] + 3 * 0);
  v01.CopyToScale(dCdp[2] + 3 * 1, -1.0);
}

/**
 *
 * @param C
 * @param dCdp
 * @param[in] P undeformed triangle vertex positions
 * @param[in] p deformed triangle vertex positions
 * @param L01
 */
DFM2_INLINE void delfem2::CdC_LengthTetEdges(
  double C[6],
  double dCdp[6][12],
  const double P[4][3],
  const double p[4][3]) {
  const double L01 = Distance3(P[0], P[1]);
  const double L02 = Distance3(P[0], P[2]);
  const double L03 = Distance3(P[0], P[3]);
  const double L12 = Distance3(P[1], P[2]);
  const double L13 = Distance3(P[1], P[3]);
  const double L23 = Distance3(P[2], P[3]);
  CVec3d v01(p[0][0] - p[1][0], p[0][1] - p[1][1], p[0][2] - p[1][2]);
  CVec3d v02(p[0][0] - p[2][0], p[0][1] - p[2][1], p[0][2] - p[2][2]);
  CVec3d v03(p[0][0] - p[3][0], p[0][1] - p[3][1], p[0][2] - p[3][2]);
  CVec3d v12(p[1][0] - p[2][0], p[1][1] - p[2][1], p[1][2] - p[2][2]);
  CVec3d v13(p[1][0] - p[3][0], p[1][1] - p[3][1], p[1][2] - p[3][2]);
  CVec3d v23(p[2][0] - p[3][0], p[2][1] - p[3][1], p[2][2] - p[3][2]);
  const double l01 = v01.norm();
  const double l02 = v02.norm();
  const double l03 = v03.norm();
  const double l12 = v12.norm();
  const double l13 = v13.norm();
  const double l23 = v23.norm();
  C[0] = l01 - L01;
  C[1] = l02 - L02;
  C[2] = l03 - L03;
  C[3] = l12 - L12;
  C[4] = l13 - L13;
  C[5] = l23 - L23;
  // ----
  v01 /= l01;
  v02 /= l02;
  v03 /= l03;
  v12 /= l12;
  v13 /= l13;
  v23 /= l23;
  // ----
  for (int i = 0; i < 6 * 3 * 4; ++i) { (&dCdp[0][0])[i] = 0.0; }
  v01.CopyTo(dCdp[0] + 0 * 3);
  v01.CopyToScale(dCdp[0] + 1 * 3, -1.0);
  v02.CopyTo(dCdp[1] + 0 * 3);
  v02.CopyToScale(dCdp[1] + 2 * 3, -1.0);
  v03.CopyTo(dCdp[2] + 0 * 3);
  v03.CopyToScale(dCdp[2] + 3 * 3, -1.0);
  v12.CopyTo(dCdp[3] + 1 * 3);
  v12.CopyToScale(dCdp[3] + 2 * 3, -1.0);
  v13.CopyTo(dCdp[4] + 1 * 3);
  v13.CopyToScale(dCdp[4] + 3 * 3, -1.0);
  v23.CopyTo(dCdp[5] + 2 * 3);
  v23.CopyToScale(dCdp[5] + 3 * 3, -1.0);
}

