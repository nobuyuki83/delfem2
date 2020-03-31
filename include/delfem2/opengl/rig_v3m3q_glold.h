/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef GLOLD_RIG_V23Q_H
#define GLOLD_RIG_V23Q_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>

#include "delfem2/vec3.h"
#include "delfem2/rig_v3q.h"

namespace delfem2 {
namespace opengl{

void Draw_RigBone(
    int ibone,
    bool is_selected,
    int ielem_selected,
    const std::vector<delfem2::CRigBone>& aBone,
    double rad_bone_sphere,
    double rad_rot_hndlr);

void DrawBone(const std::vector<delfem2::CRigBone>& aBone,
    int ibone_selected,
    int ielem_selected,
    double rad_bone_sphere,
    double rad_rot_hndlr);

void DrawJoints(
    const std::vector<double>& aJntPos,
    const std::vector<int>& aIndBoneParent);

}
}


#endif /* rigmesh_hpp */
