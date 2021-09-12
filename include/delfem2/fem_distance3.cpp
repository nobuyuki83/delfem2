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
    CVec3d dW_dP[2],
    CMat3d ddW_ddP[2][2],
    const double stiff,
    const CVec3d p[2],
    double L0) {
  const double l = sqrt(
      +(p[0][0] - p[1][0]) * (p[0][0] - p[1][0])
          + (p[0][1] - p[1][1]) * (p[0][1] - p[1][1])
          + (p[0][2] - p[1][2]) * (p[0][2] - p[1][2]));
  const CVec3d v = p[0] - p[1];
  const double R = L0 - l;
  dW_dP[0] = (-R * stiff / l) * v;
  dW_dP[1] = (+R * stiff / l) * v;
  const CMat3d m = (stiff * L0 / (l * l * l)) * Mat3_OuterProduct(v, v) + (stiff * (l - L0) / l) * Mat3_Identity(1.0);
  ddW_ddP[0][0] = +m;
  ddW_ddP[0][1] = -m;
  ddW_ddP[1][0] = -m;
  ddW_ddP[1][1] = +m;
  return 0.5 * stiff * R * R;
}

template<typename T>
DFM2_INLINE T delfem2::WdW_SquareLengthLineseg3D(
    T dW_dP[2][3],
    const T stiff,
    const T p[2][3],
    T L0) {
  const double l = sqrt(
      (p[0][0] - p[1][0]) * (p[0][0] - p[1][0]) +
          (p[0][1] - p[1][1]) * (p[0][1] - p[1][1]) +
          (p[0][2] - p[1][2]) * (p[0][2] - p[1][2]));
  const double R = L0 - l;
  dW_dP[0][0] = (-R * stiff / l) * (p[0][0] - p[1][0]);
  dW_dP[0][1] = (-R * stiff / l) * (p[0][1] - p[1][1]);
  dW_dP[0][2] = (-R * stiff / l) * (p[0][2] - p[1][2]);
  dW_dP[1][0] = -dW_dP[0][0];
  dW_dP[1][1] = -dW_dP[0][1];
  dW_dP[1][2] = -dW_dP[0][2];
  return 0.5 * stiff * R * R;
}

