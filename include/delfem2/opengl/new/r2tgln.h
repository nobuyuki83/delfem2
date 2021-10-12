/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_OPENGL_NEW_R2TGLN_H
#define DFM2_OPENGL_NEW_R2TGLN_H

#include <stdio.h>
#include <vector>

#include "delfem2/opengl/new/mshcolor.h" // shader definitions
#include "delfem2/opengl/new/shdr_msh.h"
#include "delfem2/opengl/new/shdr_mshtex.h"
#include "delfem2/opengl/new/shdr_points.h"
#include "delfem2/opengl/r2t.h"
#include "delfem2/dfm2_inline.h"

namespace delfem2{
namespace opengl{

class CRender2Tex_DrawNewGL
{
public:
  CRender2Tex_DrawNewGL(){
    pointSize = 3;
    isDrawTex = true;
    draw_len_axis = 1.0;
  }
  // --------------
  void Draw(
      const ::delfem2::opengl::CRender2Tex& r2t,
      const float mat4_projection[16],
      const float mat4_modelview[16]) const;
  // ------------
  virtual void InitGL(); // override function
  virtual void SetDepth(const ::delfem2::opengl::CRender2Tex& r2t);
public:
  bool isDrawTex;
  double draw_len_axis;
  float pointSize;
  CShader_Mesh drawer_view_frustrum;
  Drawer_RectangleTex drawer_projected_image{-1,+1,-1,+1,+1};
  Drawer_Coords drawer_projected_points;
};
  
}
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/new/r2tgln.cpp"
#endif

#endif /* DFM2_OPENGL_NEW_R2TGLN_H */
