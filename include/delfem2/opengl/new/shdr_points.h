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
#include <stdio.h>
#include <vector>

// -------------------------------------

namespace delfem2 {
namespace opengl {

class CShader_Points{
public:
  void Compile();

  template <typename REAL>
  void Initialize(
      std::vector<REAL>& aXYZd);

  template <typename REAL>
  void UpdateVertex(
      std::vector<REAL>& aXYZd);

  void Draw(
      float mP[16],
      float mMV[16]) const;
public:
  CGL4_VAO_Mesh vao; // gl4
  int shaderProgram;
  int Loc_MatrixProjection;
  int Loc_MatrixModelView;
  int Loc_Color;
  unsigned int nPoint = 0;
  delfem2::CColor color_face = delfem2::CColor(0.0,0.0,0.0,0.0);
};

}
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/new/shdr_points.cpp"
#endif

#endif /* gl4_msh_hpp */
