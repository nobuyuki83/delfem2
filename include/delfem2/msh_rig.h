//
// Created by Nobuyuki Umetani on 2021/12/16.
//

#ifndef DFM2_MSH_RIG_H_
#define DFM2_MSH_RIG_H_

#include "delfem2/msh_arrow.h"

namespace delfem2 {

template<typename REAL>
void Mesh_RigBones_Octahedron(
  std::vector<REAL> &vtx_xyz,
  std::vector<int> &tri_vtx,
  const std::vector<delfem2::CRigBone> &aBone) {
  namespace dfm2 = delfem2;
  for (unsigned int ibone = 0; ibone < aBone.size(); ++ibone) {
    const CRigBone &bone = aBone[ibone];
    const int ibone_p = aBone[ibone].ibone_parent;
    if (ibone_p < 0 || ibone_p >= (int) aBone.size()) { continue; }
    const CRigBone &bone_p = aBone[ibone_p];
    const CVec3d p0(bone_p.invBindMat[3], bone_p.invBindMat[7], bone_p.invBindMat[11]);
    const CVec3d p1(bone.invBindMat[3], bone.invBindMat[7], bone.invBindMat[11]);
    dfm2::CMat4d At(bone_p.affmat3Global);
    Mesh_ArrowOcta<float>(
      vtx_xyz, tri_vtx,
      At.cast<float>(),
      dfm2::CVec3f(0.f, 0.f, 0.f),
      (p0 - p1).cast<float>(),
      0.1f,
      0.2f);
  }
}

}

#endif //DFM2_MSH_RIG_H_
