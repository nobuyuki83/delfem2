#ifndef DFM2_OPENGL_NEW_SHDR_MSHTEX_H
#define DFM2_OPENGL_NEW_SHDR_MSHTEX_H

#include <cstdio>
#include <vector>

#include "delfem2/opengl/new/funcs.h" // CGL4_VAO_Mesh
#include "delfem2/dfm2_inline.h"

// -------------------------------------

namespace delfem2::opengl {

class CShader_MeshTex {
public:
  void SetElement(
      std::vector<unsigned int> &elem_vtx,
      int gl_primitive_type);

  template <typename REAL>
  void setCoords(
      std::vector<REAL> &aXYZd,
      unsigned int ndim);

  template <typename REAL>
  void setTexCoords(
      std::vector<REAL> &aTex);

  void InitGL();

  void Draw(const float mat4_projection[16],
            const float mat4_modelview[16]) const;

public:
  VertexArrayObject vao; // gl4
  int shaderProgram = -1;
  int Loc_MatrixProjection = -1;
  int Loc_MatrixModelView = -1;
  int Loc_Texture = -1;
};

}


#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/new/shdr_mshtex.cpp"
#endif



#endif
