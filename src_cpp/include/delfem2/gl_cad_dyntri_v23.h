/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef cad_dyntri_v23_gl_h
#define cad_dyntri_v23_gl_h


#include <stdio.h>
#include <vector>

#include "delfem2/dyntri.h"
#include "delfem2/vec2.h"
#include "delfem2/vec3.h"
#include "delfem2/cad2d.h"

/////////////////////////////////////////

void DrawMeshDynTri_FaceNorm(const std::vector<ETri>& aSTri,
                             const std::vector<CVector3>& aVec3);

void DrawMeshDynTri_Edge(const std::vector<ETri>& aSTri,
                         const std::vector<CVector3>& aVec3);


void DrawMeshDynTri_Edge(const std::vector<ETri>& aSTri,
                         const std::vector<CVector2>& aVec2);

void DrawMeshDynTri_FaceNorm(const std::vector<ETri>& aSTri,
                             const std::vector<CVector2>& aVec2);

void DrawMeshDynTri3D_Edge(const std::vector<double>& aXYZ,
                           const std::vector<ETri>& aSTri);

/////////////////////////////////////////

void Draw_CCad2D(const CCad2D& cad2d);


#endif /* cad_dyntri_v23_gl_hpp */
