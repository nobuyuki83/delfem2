/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @brief classes and functions that depend on cad,dtri,vec2,vec3 and OpenGL ver. 4
 */

#ifndef DFM2_OPENGL_NEW_DYNMESH_H
#define DFM2_OPENGL_NEW_DYNMESH_H

#include "delfem2/opengl/new/funcs.h" // for CGL4_VAO_Mesh
#include "delfem2/vec2.h"
#include "delfem2/dtri.h"
#include "delfem2/dfm2_inline.h"

namespace delfem2::opengl{

class CShader_MeshDTri2D
{
public:
  CShader_MeshDTri2D(){
  }
  void MakeBuffer(const std::vector<CVec2d>& aVec2,
                  const std::vector<CDynTri>& aETri);
  void Draw(const float mP[16],
            const float mMV[16]) const;
  void Compile();
public:
  int shdr0_program; // for face
  int shdr0_Loc_MatrixProjection;
  int shdr0_Loc_MatrixModelView;
  int shdr0_Loc_Color;
  
  VertexArrayObject vao;
};
  
} // namespace delfem2::opnegl


#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/new/drawer_dynmesh.cpp"
#endif

#endif /* DFM2_OPENGL_NEW_V23DTRICAD_H */
