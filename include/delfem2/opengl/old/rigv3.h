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

#include "delfem2/opengl/old/gizmo.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/dfm2_inline.h"

namespace delfem2::opengl {

DFM2_INLINE void DrawBoneReference_Line(
    const std::vector<delfem2::CRigBone> &aBone,
    double rad_bone_sphere );

DFM2_INLINE void DrawBoneCurrent_Line(
    const std::vector<delfem2::CRigBone>& aBone,
    int ibone_selected,
    double rad_bone_sphere );

DFM2_INLINE void DrawBoneCurrent_Octahedron(
    const std::vector<CRigBone>& aBone,
    unsigned int ibone_selected,
    unsigned int ielem_selected,
    double rad_bone_sphere,
    double rad_rot_hndlr);

DFM2_INLINE void DrawBoneReference_Octahedron(
    const std::vector<CRigBone> &aBone,
    unsigned int ibone_selected,
    [[maybe_unused]] unsigned int ielem_selected,
    [[maybe_unused]] double rad_bone_sphere,
    [[maybe_unused]] double rad_rot_hndlr);

DFM2_INLINE void DrawJoints(
    const std::vector<double>& aJntPos,
    const std::vector<int>& aIndBoneParent);

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/old/rigv3.cpp"
#endif


#endif /* RIGV3_GLOLD_H */
