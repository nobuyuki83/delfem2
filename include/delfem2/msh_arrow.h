//
// Created by Nobuyuki Umetani on 2021/12/16.
//

#ifndef MSH_ARROW_H_
#define MSH_ARROW_H_

#include "delfem2/mat4.h"
#include "delfem2/vec3.h"

namespace delfem2 {

template<typename REAL, typename INT>
DFM2_INLINE void Mesh_ArrowOcta(
  std::vector<REAL> &vtx_xyz,
  std::vector<INT> &tri_vtx,
  const CMat4<REAL> &affine_trans,
  const CVec3<REAL> &p0,
  const CVec3<REAL> &d,
  REAL rad_ratio,
  REAL node_ratio) {
  using CV3 = CVec3<REAL>;
  if (d.norm() < 1.0e-10) { return; }
  CV3 z = d;
  z.normalize();
  CV3 x, y;
  GetVertical2Vector(z, x, y);
  const REAL dt = M_PI * 0.5;
  const REAL r0 = d.norm() * rad_ratio;
  const CV3 p1 = p0 + node_ratio * d;
  const CV3 p2 = p0 + d;
  //
  const size_t i0 = vtx_xyz.size() / 3;
  CVec3f(affine_trans.MultVec3_Homography(p0.data())).push_back_to_vector(vtx_xyz);
  CVec3f(affine_trans.MultVec3_Homography(p2.data())).push_back_to_vector(vtx_xyz);
  for (int idiv = 0; idiv < 4; idiv++) {
    const REAL s0 = (REAL) (r0 * std::sin(idiv * dt));
    const REAL c0 = (REAL) (r0 * std::cos(idiv * dt));
    const CV3 q0 = p1 + s0 * x + c0 * y;
    CVec3f(affine_trans.MultVec3_Homography(q0.data())).push_back_to_vector(vtx_xyz);
  }
  for (int idiv = 0; idiv < 4; idiv++) {
    tri_vtx.push_back(i0);
    tri_vtx.push_back(i0 + 2 + (1 + idiv) % 4);
    tri_vtx.push_back(i0 + 2 + (0 + idiv) % 4);
    //
    tri_vtx.push_back(i0 + 1);
    tri_vtx.push_back(i0 + 2 + (0 + idiv) % 4);
    tri_vtx.push_back(i0 + 2 + (1 + idiv) % 4);
  }
}

}

#endif //MSH_ARROW_H_
