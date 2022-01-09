/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_OPENGL_NEW_DRAWER_MSHDISP_H
#define DFM2_OPENGL_NEW_DRAWER_MSHDISP_H

#include "delfem2/opengl/new/funcs.h" // CGL4_VAO_Mesh
#include "delfem2/color.h"
#include "delfem2/dfm2_inline.h"
#include <stdio.h>
#include <vector>

// -------------------------------------

namespace delfem2::opengl {

class CShader_TriMesh_Disp{
public:
  CShader_TriMesh_Disp(){
  }

  template <typename REAL>
  void Initialize(
      std::vector<REAL>& aPosD,
      unsigned int ndim,
      std::vector<unsigned int>& aTri,
      std::vector<REAL>& aDispD);

  template <typename REAL>
  void UpdateVertex(
      std::vector<REAL>& aPosD,
      unsigned int ndim,
      std::vector<REAL>& aDispD);

  void Compile();
  void Draw(float mP[16], float mMV[16]);
  
public:
  VertexArrayObject vao; // gl4
  int shaderProgram;
  int Loc_MatrixProjection;
  int Loc_MatrixModelView;
  int Loc_Color0, Loc_Color1;
  int Loc_ValMin, Loc_ValMax;
};

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/new/drawer_mshdisp.cpp"
#endif

#endif /* DFM2_OPENGL_NEW_DRAWER_MSHDISP_H */
