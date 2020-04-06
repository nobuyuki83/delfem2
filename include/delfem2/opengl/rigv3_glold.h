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
#include "delfem2/vec3.h"
#include "delfem2/rig_v3q.h"

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

}
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/rigv3_glold.cpp"
#endif


#endif /* rigmesh_hpp */
