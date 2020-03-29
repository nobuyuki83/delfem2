/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_CADDTRI_V3_GLOLD_H
#define DFM2_CADDTRI_V3_GLOLD_H

#include <stdio.h>
#include <vector>

#include "delfem2/dtri.h"
#include "delfem2/vec3.h"
#include "delfem2/cad2d_v2dtri.h"

// ---------------------------

namespace delfem2{
namespace opengl{


void DrawMeshDynTri_FaceNorm(const std::vector<CDynTri>& aSTri,
                             const std::vector<CVec3d>& aVec3);

void DrawMeshDynTri_FaceNorm(const std::vector<CDynTri>& aSTri,
                             const double* aXYZ);

void DrawMeshDynTri_Edge(const std::vector<CDynTri>& aSTri,
                         const std::vector<CVec3d>& aVec3);

void DrawMeshDynTri3D_Edge(const std::vector<double>& aXYZ,
                           const std::vector<CDynTri>& aSTri);

}
}

#endif /* cad_dyntri_v23_gl_hpp */
