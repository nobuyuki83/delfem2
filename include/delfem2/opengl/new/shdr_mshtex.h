#ifndef DFM2_OPENGL_NEW_SHDR_MSHTEX_H
#define DFM2_OPENGL_NEW_SHDR_MSHTEX_H

#include "delfem2/opengl/new/funcs.h" // CGL4_VAO_Mesh
#include "delfem2/dfm2_inline.h"
#include <stdio.h>
#include <vector>

// -------------------------------------

namespace delfem2 {
namespace opengl {

class CShader_MeshTex {
public:
  void setElement(
      std::vector<unsigned int> &aTri,
      int gl_primitive_type);

  template <typename REAL>
  void setCoords(
      std::vector<REAL> &aXYZd,
      unsigned int ndim);

  template <typename REAL>
  void setTexCoords(
      std::vector<REAL> &aTex);

  void Compile();

  void Draw(float mP[16], float mMV[16]) const;

public:
  CGL4_VAO_Mesh vao; // gl4
  int shaderProgram;
  int Loc_MatrixProjection;
  int Loc_MatrixModelView;
  int Loc_Texture;
};

}
}


#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/opengl/new/shdr_mshtex.cpp"
#endif



#endif
