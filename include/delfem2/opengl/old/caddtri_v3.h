/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_CADDTRI_V3_GLOLD_H
#define DFM2_CADDTRI_V3_GLOLD_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/cad2_dtri2.h"
#include "delfem2/vec3.h"
#include <stdio.h>
#include <vector>


// ---------------------------

namespace delfem2{
namespace opengl{


DFM2_INLINE void DrawMeshDynTri_FaceNorm(
    const std::vector<CDynTri>& aSTri,
    const std::vector<CVec3d>& aVec3);

DFM2_INLINE void DrawMeshDynTri_FaceNorm(
    const std::vector<CDynTri>& aSTri,
    const double* aXYZ);

DFM2_INLINE void DrawMeshDynTri_FaceNormTex(
    const std::vector<CDynTri>& aSTri,
    const double* aXYZ,
    const std::vector<CVec2d>& aVec2);

DFM2_INLINE void DrawMeshDynTri_Edge(
    const std::vector<CDynTri>& aSTri,
    const std::vector<CVec3d>& aVec3);

DFM2_INLINE void DrawMeshDynTri3D_Edge(
    const std::vector<double>& aXYZ,
    const std::vector<CDynTri>& aSTri);

}
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/old/caddtri_v3.cpp"
#endif

#endif /* cad_dyntri_v23_gl_hpp */
