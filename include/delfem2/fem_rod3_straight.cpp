/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/fem_rod3_straight.h"

#include "delfem2/geo3_v23m34q.h"
#include "delfem2/mat3.h"

// -------------------------------------------------

template <typename T>
void delfem2::CdC_Rod3BendStraight(
    T C[3],
    T dC_dP[3][3][3],
    const T vtx_xyz[3][3]) {
  namespace dfm2 = delfem2;
  using V3 = delfem2::CVec3<T>;
  using M3 = delfem2::CMat3<T>;
  const V3 v0 = V3(vtx_xyz[1]) - V3(vtx_xyz[0]);
  const V3 v1 = V3(vtx_xyz[2]) - V3(vtx_xyz[1]);
  const T l0 = v0.norm();
  const T l2 = v1.norm();
  const V3 u0 = v0 / l0;
  const V3 u2 = v1 / l2;
  const M3 du0dp0 = -(M3::Identity() - M3::OuterProduct(u0.p, u0.p)) / l0;
  const M3 du2dp2 = +(M3::Identity() - M3::OuterProduct(u2.p, u2.p)) / l2;
  const T cos0 = 1 + u0.dot(u2);
  const V3 tan0 = u0.cross(u2) / cos0;
  const M3 dtdu0 = (-M3::Spin(u2.p) - M3::OuterProduct(tan0.p, u2.p)) / cos0;
  const M3 dtdu2 = (+M3::Spin(u0.p) - M3::OuterProduct(tan0.p, u0.p)) / cos0;
  const M3 dtdp0 = dtdu0 * du0dp0;
  const M3 dtdp2 = dtdu2 * du2dp2;
  tan0.CopyTo(C);
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j) {
      dC_dP[i][0][j] = dtdp0(i,j);
      dC_dP[i][2][j] = dtdp2(i,j);
      dC_dP[i][1][j] = -dC_dP[i][0][j] - dC_dP[i][2][j];
    }
  }
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::CdC_Rod3BendStraight(
    float C[3],
    float dC_dP[3][3][3],
    const float vtx_xyz[3][3]);
template void delfem2::CdC_Rod3BendStraight(
    double C[3],
    double dC_dP[3][3][3],
    const double vtx_xyz[3][3]);
#endif

// --------------------------------------------

double delfem2::WdWddW_Rod3BendStraight(
    delfem2::CVec3d dW_dP[3],
    delfem2::CMat3d ddW_ddP[3][3],
    const delfem2::CVec3d vec_pos[3]) {
  namespace dfm2 = delfem2;
  const double l0 = (vec_pos[1] - vec_pos[0]).norm();
  const double l2 = (vec_pos[2] - vec_pos[1]).norm();
  const dfm2::CVec3d u0 = (vec_pos[1] - vec_pos[0]) / l0;
  const dfm2::CVec3d u2 = (vec_pos[2] - vec_pos[1]) / l2;
  const dfm2::CMat3d du0 = -(dfm2::CMat3d::Identity() - dfm2::CMat3d::OuterProduct(u0.p, u0.p)) / l0;
  const dfm2::CMat3d du2 = +(dfm2::CMat3d::Identity() - dfm2::CMat3d::OuterProduct(u2.p, u2.p)) / l2;
  const double c = 1 + u0.dot(u2);
  const dfm2::CVec3d t = u0.cross(u2) / c;
  const double w = 0.5 * t.squaredNorm();
  const dfm2::CMat3d dtdu0 =
      (-dfm2::CMat3d::Spin(u2.p) - dfm2::CMat3d::OuterProduct(t.p, u2.p)) / c;
  const dfm2::CMat3d dtdu2 =
      (+dfm2::CMat3d::Spin(u0.p) - dfm2::CMat3d::OuterProduct(t.p, u0.p)) / c;
  dW_dP[0] = +t * dtdu0 * du0;
  dW_dP[2] = +t * dtdu2 * du2;
  dW_dP[1] = -(dW_dP[0] + dW_dP[2]);
  {
    const dfm2::CMat3d tddtddu0 =
        -dfm2::CMat3d::Spin(u2.p) * dfm2::CMat3d::OuterProduct(t.p, u2.p) / (c * c) +
            dfm2::CMat3d::OuterProduct(u2.p, t.p) * dfm2::CMat3d::Spin(u2.p) / (c * c) +
            2 * dfm2::CMat3d::OuterProduct(u2.p, u2.p) * t.squaredNorm() / (c * c);
    const dfm2::CVec3d x = +t * dtdu0;
    const dfm2::CMat3d xddu0 =
        (dfm2::CMat3d::OuterProduct(u0.p, u0.p) * 3. * u0.dot(x)
            - dfm2::CMat3d::OuterProduct(x.p, u0.p)
            - dfm2::CMat3d::OuterProduct(u0.p, x.p)
            - dfm2::CMat3d::Identity() * u0.dot(x)
        ) / (l0 * l0);
    ddW_ddP[0][0] = du0.transpose() * dtdu0.transpose() * dtdu0 * du0 +
        du0.transpose() * tddtddu0 * du0 +
        xddu0;
  }
  {
    const dfm2::CMat3d tddtdu2du0 =
        -dfm2::CMat3d::Spin(t.p) / c
            - (t.squaredNorm() / c) * dfm2::CMat3d::Identity()
            + dfm2::CMat3d::Spin(t.p) * dfm2::CMat3d::OuterProduct(u2.p, u0.p) / (c * c)
            + dfm2::CMat3d::OuterProduct(u2.p, u0.p) * dfm2::CMat3d::Spin(t.p) / (c * c)
            + 2 * dfm2::CMat3d::OuterProduct(u2.p, u0.p) * t.squaredNorm() / (c * c);
    ddW_ddP[2][0] = du2.transpose() * dtdu2.transpose() * dtdu0 * du0
        + du2.transpose() * tddtdu2du0.transpose() * du0;
    ddW_ddP[0][2] = ddW_ddP[2][0].transpose();
  }
  {
    const dfm2::CMat3d tddtddu2 =
        dfm2::CMat3d::Spin(u0.p) * dfm2::CMat3d::OuterProduct(t.p, u0.p) / (c * c) -
            dfm2::CMat3d::OuterProduct(u0.p, t.p) * dfm2::CMat3d::Spin(u0.p) / (c * c) +
            2 * dfm2::CMat3d::OuterProduct(u0.p, u0.p) * t.squaredNorm() / (c * c);
    const dfm2::CVec3d x = +t * dtdu2;
    const dfm2::CMat3d xddu2 =
        (dfm2::CMat3d::OuterProduct(u2.p, u2.p) * 3. * u2.dot(x)
            - dfm2::CMat3d::OuterProduct(x.p, u2.p)
            - dfm2::CMat3d::OuterProduct(u2.p, x.p)
            - dfm2::CMat3d::Identity() * u2.dot(x)
        ) / (l2 * l2);
    ddW_ddP[2][2] = du2.transpose() * dtdu2.transpose() * dtdu2 * du2 +
        du2.transpose() * tddtddu2 * du2 +
        xddu2;
  }
  ddW_ddP[0][1] = -(ddW_ddP[0][0] + ddW_ddP[0][2]);
  ddW_ddP[1][0] = ddW_ddP[0][1].transpose();
  ddW_ddP[2][1] = -(ddW_ddP[2][0] + ddW_ddP[2][2]);
  ddW_ddP[1][2] = ddW_ddP[2][1].transpose();
  ddW_ddP[1][1] = -(ddW_ddP[1][0] + ddW_ddP[1][2]);
  return w;
}
