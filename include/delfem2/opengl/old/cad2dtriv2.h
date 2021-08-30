/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_OPENGL_OLD_CAD2DTRIV2_H
#define DFM2_OPENGL_OLD_CAD2DTRIV2_H

#include <vector>

#include "delfem2/dfm2_inline.h"
#include "delfem2/cad2_dtri2.h"

// ---------------------------

namespace delfem2{
namespace opengl{

DFM2_INLINE void DrawMeshDynTri_FaceNorm(
    const std::vector<CDynTri>& aSTri,
    const double* aXYZ);

DFM2_INLINE void DrawMeshDynTri_Edge(
    const std::vector<CDynTri>& aSTri,
    const std::vector<CVec2d>& aVec2);

DFM2_INLINE void DrawMeshDynTri_FaceNorm(
    const std::vector<CDynTri>& aSTri,
    const std::vector<CVec2d>& aVec2);

// --------------------------------------------

DFM2_INLINE void Draw_CCad2DEdge(
    const delfem2::CCad2D_EdgeGeo& edge,
    bool is_selected,
    int ipicked_elem);

DFM2_INLINE void Draw_CCad2D(
    const delfem2::CCad2D& cad2d);
  
}
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/old/cad2dtriv2.cpp"
#endif

#endif
