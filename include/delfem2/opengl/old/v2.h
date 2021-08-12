/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_OPENGL_OLD_V2_H
#define DFM2_OPENGL_OLD_V2_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec2.h"
#include <vector>

namespace delfem2 {
namespace opengl {

// ------------------------------------------------------------------------------------
// vec2 starts here

DFM2_INLINE void myGlVertex(
    unsigned int i,
    const std::vector<CVec2d> &aP);

DFM2_INLINE void myGlVertex(
    const CVec2d &v);

DFM2_INLINE void drawPolyLine(
    const std::vector<CVec2d> &aP);

DFM2_INLINE void drawPolyLine2D(
    const std::vector<CVec2d> &aP);

DFM2_INLINE void Draw_MeshTri(
    const std::vector<CVec2d> &aP,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE void Draw_MeshTri_Edge(
    const std::vector<CVec2d> &aP,
    const std::vector<unsigned int> &aTri);

DFM2_INLINE void myGlVertex2(
    int i,
    const std::vector<double> &vec);

} // end namespace opengl
} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#include "delfem2/opengl/old/v2.cpp"
#endif

#endif
