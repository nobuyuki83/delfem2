//
// Created by Nobuyuki Umetani on 2021/11/29.
//

#ifndef DFM2_DRAWER_MESHUNINDEX_H_
#define DFM2_DRAWER_MESHUNINDEX_H_

#include "delfem2/opengl/funcs.h"
#include "delfem2/opengl/new/funcs.h"

namespace delfem2::opengl {

class Drawer_MeshUnIndexed{
 public:
  void InitGL(){

    const std::string glsl33vert_projection =
      "uniform mat4 matrixProjection;\n"
      "uniform mat4 matrixModelView;\n"
      "layout (location = 0) in vec3 posIn;\n"
      "layout (location = 1) in vec3 nrmIn;\n"
      "layout (location = 2) in vec3 rgbIn;\n"
      "out vec3 nrmPrj;\n"
      "out vec3 col;\n"
      "void main()\n"
      "{\n"
      "  gl_Position = matrixProjection * matrixModelView * vec4(posIn.x, posIn.y, posIn.z, 1.0);\n"
      "  vec4 v0 = matrixModelView * vec4(nrmIn.x, nrmIn.y, nrmIn.z, 0.0);\n"
      "  nrmPrj = v0.xyz;\n"
      "  if( length(nrmIn) < 1.e-30 ){ nrmPrj = vec3(0.f, 0.f, 1.f); }\n"
      "  col = rgbIn;\n"
      "  nrmPrj = normalize(nrmPrj);\n"
      "}";

    const std::string glsl33frag =
      "uniform vec3 color;\n"
      "in vec3 nrmPrj;\n"
      "in vec3 col;\n"
      "out vec4 FragColor;\n"
      "void main()\n"
      "{\n"
      "  FragColor = abs(nrmPrj.z)*vec4(col,1);\n"
      "}";

#ifdef EMSCRIPTEN
    shaderProgram = GL24_CompileShader(
      (std::string("#version 300 es\n")+
      glsl33vert_projection).c_str(),
      (std::string("#version 300 es\n")+
      std::string("precision highp float;\n")+
      glsl33frag).c_str());
#else
    shaderProgram = delfem2::opengl::GL24_CompileShader(
      (std::string("#version 330 core\n")+
        glsl33vert_projection).c_str(),
      (std::string("#version 330 core\n")+
        glsl33frag).c_str());
#endif

    if( !glIsProgram(shaderProgram) ){
      std::cout << "shader doesnot exist" << std::endl;
    }
    glUseProgram(shaderProgram);
    Loc_MatrixProjection = glGetUniformLocation(shaderProgram,  "matrixProjection");
    Loc_MatrixModelView  = glGetUniformLocation(shaderProgram,  "matrixModelView");
  }

  template <typename REAL>
  void Initialize(
    std::vector<REAL>& elem_xyz,
    std::vector<REAL>& elem_nrm,
    std::vector<REAL>& elem_rgb,
    unsigned int ndim) {

    np = elem_xyz.size() / ndim;

    if( !glIsVertexArray(vao.idx_vao) ){ glGenVertexArrays(1, &vao.idx_vao); }
    vao.Delete_EBOs();

    glBindVertexArray(vao.idx_vao); // opengl4
    vao.ADD_VBO(0, elem_xyz);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(
      0,
      ndim, convertToGlType<REAL>(), GL_FALSE, ndim*sizeof(REAL),
      (void*)nullptr);
    // ---
    vao.ADD_VBO(1, elem_nrm);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(
      1,
      ndim, convertToGlType<REAL>(), GL_FALSE, ndim*sizeof(REAL),
      (void*)nullptr);
    // ---
    vao.ADD_VBO(2, elem_rgb);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(
      2,
      ndim, convertToGlType<REAL>(), GL_FALSE, ndim*sizeof(REAL),
      (void*)nullptr);
  }


  void Draw(const float mat4_projection[16],
            const float mat4_modelview[16]) const
  {
    glUseProgram(shaderProgram);
    glUniformMatrix4fv(
      Loc_MatrixProjection, 1, GL_FALSE,
      TransposeMat4ForOpenGL(mat4_projection, true).data());
    glUniformMatrix4fv(
      Loc_MatrixModelView, 1, GL_FALSE,
      TransposeMat4ForOpenGL(mat4_modelview, false).data());
    vao.DrawArray(GL_TRIANGLES,np);
  }

  void UpdateTriRgb(
    const std::vector<double>& tri_rgb){
    vao.ADD_VBO(2,tri_rgb);
  }

 public:
  VertexArrayObject vao; // gl4
  int shaderProgram;
  int Loc_MatrixProjection;
  int Loc_MatrixModelView;
  unsigned int np;
};

}



#endif //DFM2_DRAWER_MESHUNINDEX_H_
