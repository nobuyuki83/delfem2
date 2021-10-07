/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_OPENGL_NEW_SHDR_POINTS_H
#define DFM2_OPENGL_NEW_SHDR_POINTS_H

#include "delfem2/opengl/new/funcs.h" // CGL4_VAO_Mesh
#include "delfem2/color.h"
#include "delfem2/dfm2_inline.h"
#include <cstdio>
#include <vector>

// -------------------------------------

namespace delfem2::opengl {

class CShader_Points{
public:
  void InitGL();

  template <typename REAL>
  void SetCoords(std::vector<REAL> &vtx_coords, unsigned int ndim);

  void Draw(
      GLenum gl_primitive_type,
      float mP[16],
      float mMV[16]) const;
public:
  VertexArrayObject vao; // gl4
  int shaderProgram = -1;
  int Loc_MatrixProjection = -1;
  int Loc_MatrixModelView = -1;
  int Loc_Color = -1;
  unsigned int num_vtx = 0;
  delfem2::CColor color_face = delfem2::CColor(0.0,0.0,0.0,0.0);
};

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/new/shdr_points.cpp"
#endif

#endif /* gl4_msh_hpp */
