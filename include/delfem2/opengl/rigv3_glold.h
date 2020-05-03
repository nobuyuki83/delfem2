/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef RIGV3_GLOLD_H
#define RIGV3_GLOLD_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include "delfem2/dfm2_inline.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/opengl/gizmo_glold.h"

namespace delfem2 {
namespace opengl{

DFM2_INLINE void Draw_RigBone(
    int ibone,
    bool is_selected,
    int ielem_selected,
    const std::vector<delfem2::CRigBone>& aBone,
    double rad_bone_sphere,
    double rad_rot_hndlr);

DFM2_INLINE void DrawBone(
    const std::vector<delfem2::CRigBone>& aBone,
    int ibone_selected,
    int ielem_selected,
    double rad_bone_sphere,
    double rad_rot_hndlr);

DFM2_INLINE void DrawJoints(
    const std::vector<double>& aJntPos,
    const std::vector<int>& aIndBoneParent);

void Draw
 (CGizmo_Rig<float>& giz,
  const std::vector<CRigBone>& aBone)
{
  if( giz.mode_edit == CGizmo_Rig<float>::MODE_EDIT::TRNSL ){ // translation
    giz.gizmo_trnsl.pos[0] = aBone[0].transRelative[0];
    giz.gizmo_trnsl.pos[1] = aBone[0].transRelative[1];
    giz.gizmo_trnsl.pos[2] = aBone[0].transRelative[2];
    opengl::Draw(giz.gizmo_trnsl);
  }
  else if( giz.mode_edit == CGizmo_Rig<float>::MODE_EDIT::ROT ){ // translation
    if( giz.ipicked_bone != -1 ){
      assert( giz.ipicked_bone >= 0 && giz.ipicked_bone < (int)aBone.size() );
      giz.gizmo_rot.pos = aBone[giz.ipicked_bone].Pos().Float();
      { // set quaternion
        CMat3<double> m3;
        m3.SetMat4(aBone[giz.ipicked_bone].affmat3Global);
        CQuat<double> qj;
        m3.GetQuat_RotMatrix(qj.q);
        qj.CopyTo(giz.gizmo_rot.quat);
      }
      opengl::Draw(giz.gizmo_rot);
    }
  }
}

}
}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/opengl/rigv3_glold.cpp"
#endif


#endif /* rigmesh_hpp */
