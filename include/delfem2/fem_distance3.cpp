/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/fem_distance3.h"

#include "delfem2/geo3_v23m34q.h"
#include "delfem2/mat3.h"
#include "delfem2/mshuni.h"

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
  const CMat3d m = (stiffness * edge_length_ini / (l * l * l)) * Mat3_OuterProduct(v, v)
      + (stiffness * (l - edge_length_ini) / l) * Mat3_Identity(1.0);
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

