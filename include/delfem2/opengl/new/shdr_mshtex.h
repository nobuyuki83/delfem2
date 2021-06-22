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
  void Initialize(
      std::vector<double> &aXYZd,
      std::vector<unsigned int> &aTri,
      int gl_primitive_type,
      std::vector<double> &aTex);

  void UpdateVertex(
      std::vector<double> &aXYZd,
      std::vector<unsigned int> &aTri,
      std::vector<double> &aTex);

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


#ifdef DFM2_HEADER_ONLY
#  include "delfem2/opengl/new/shdr_mshtex.cpp"
#endif



#endif
