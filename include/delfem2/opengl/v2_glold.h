/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @discussion It might be nice to separate dependency of v2 from v3 and quat
 * because 2D application does not use v3 and quat
 */



#ifndef DFM2_V2_GLOLD_H
#define DFM2_V2_GLOLD_H

#include <vector>
#include "delfem2/vec2.h"

namespace delfem2{
namespace opengl
{

// ------------------------------------------------------------------------------------
// vec2 starts here

void myGlVertex(unsigned int i,
                const std::vector<CVec2d>& aP);

void myGlVertex(const CVec2d& v);

void drawPolyLine(const std::vector<CVec2d>& aP);

void drawPolyLine2D(const std::vector<CVec2d>& aP);

void Draw_MeshTri(const std::vector<CVec2d>& aP,
                  const std::vector<unsigned int>& aTri);

void Draw_MeshTri_Edge(const std::vector<CVec2d>& aP,
                       const std::vector<unsigned int>& aTri);


// ------------------------------------------------------------------------------------
// vec3 starts here

void myGlVertex2(int i, const std::vector<double>& vec);

  
} // end namespace opengl
} // end namespace delfem2

#endif
