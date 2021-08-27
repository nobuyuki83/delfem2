/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_OPENGL_NEW_SHDR_MSH_H
#define DFM2_OPENGL_NEW_SHDR_MSH_H

#include <stdio.h>
#include <vector>

#include "delfem2/opengl/new/funcs.h" // CGL4_VAO_Mesh
#include "delfem2/color.h"
#include "delfem2/dfm2_inline.h"

// -------------------------------------

namespace delfem2 {
namespace opengl {

class CShader_Mesh {
 public:
  void Compile();

  template<typename REAL>
  void Initialize(
      std::vector<REAL> &aXYZd,
      unsigned int ndim,
      std::vector<unsigned int> &aLine,
      int gl_primitive_type);

  template<typename REAL>
  void UpdateVertex(
      std::vector<REAL> &aXYZd,
      unsigned int ndim,
      std::vector<unsigned int> &aLine);

  void Draw(float mP[16], float mMV[16]) const;

 public:
  CGL4_VAO_Mesh vao; // gl4
  float color[3] = {0, 0, 0};
  int shaderProgram;
  int Loc_MatrixProjection;
  int Loc_MatrixModelView;
  int Loc_Color;
};

}
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/new/shdr_msh.cpp"
#endif

#endif /* gl4_msh_hpp */
