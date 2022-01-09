/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/new/drawer_mshcolormap.h"

#include <sstream>
#ifdef EMSCRIPTEN
#include <emscripten/emscripten.h>
#define GLFW_INCLUDE_ES3
#include <GLFW/glfw3.h>
#elif defined(USE_GLEW)
#include <GL/glew.h>
#else
#include <glad/glad.h>
#endif

#include "delfem2/opengl/funcs.h" // compile shader
#include "delfem2/opengl/new/funcs.h" // CGL4_VAO_Mesh
#include "delfem2/mshmisc.h"

namespace dfm2 = delfem2;

// ----------------------------------------------------------------

void delfem2::opengl::Drawer_MeshColormap::AddConnectivity(
    std::vector<unsigned int> &aTri,
    int gl_primitive_type) {
  if (!glIsVertexArray(vao.idx_vao)) { glGenVertexArrays(1, &vao.idx_vao); }
  vao.Add_EBO(aTri, gl_primitive_type);
  /*
  //
  std::vector<unsigned int> aLine;
  MeshLine_MeshElem(
      aLine,
      aTri.data(), aTri.size() / 3, dfm2::MESHELEM_TRI,
      aPosD.size() / ndim);
  vao.Add_EBO(aLine, GL_LINES);
  //
  this->UpdateVertex(aPosD, ndim, aValD);
   */
}

template<typename REAL>
void delfem2::opengl::Drawer_MeshColormap::SetCoordinates(
    std::vector<REAL> &aPosD,
    [[maybe_unused]] unsigned int ndim) {
  glBindVertexArray(vao.idx_vao); // opengl4
  vao.ADD_VBO(0, aPosD);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(
      0, ndim, convertToGlType<REAL>(),
      GL_FALSE, ndim * sizeof(REAL), (void *) 0); // gl24
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::Drawer_MeshColormap::SetCoordinates(
    std::vector<float>& aPosD,
    unsigned int ndim);
template void delfem2::opengl::Drawer_MeshColormap::SetCoordinates(
    std::vector<double>& aPosD,
    unsigned int ndim);
#endif

template<typename REAL>
void delfem2::opengl::Drawer_MeshColormap::SetValues(
    std::vector<REAL> &aValD) {
  glBindVertexArray(vao.idx_vao); // opengl4
  vao.ADD_VBO(1, aValD);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(
      1, 1, convertToGlType<REAL>(),
      GL_FALSE, 1 * sizeof(REAL), (void *) 0); // gl24
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::Drawer_MeshColormap::SetValues(
    std::vector<float>& aValD);
template void delfem2::opengl::Drawer_MeshColormap::SetValues(
    std::vector<double>& aValD);
#endif

void delfem2::opengl::Drawer_MeshColormap::InitGL() {
  const std::string glsl33vert_projection =
      "uniform mat4 matrixProjection;\n"
      "uniform mat4 matrixModelView;\n"
      "layout (location = 0) in vec3 posIn;\n"
      "layout (location = 1) in float valIn;\n"
      "out float val;\n"
      "void main()\n"
      "{\n"
      "  gl_Position = matrixProjection * matrixModelView * vec4(posIn.x, posIn.y, posIn.z, 1.0);\n"
      "  val = valIn;\n"
      "}\0";

  std::string str_colors;
  {
    std::stringstream ss;
    ss << "const int ncolor = " << colors.size() << ";\n";
    ss << "vec3[ncolor] colors = vec3[] (\n";
    for (unsigned int ic = 0; ic < colors.size(); ++ic) {
      const auto &c = colors[ic];
      ss << " vec3(" << std::get<0>(c) << "," << std::get<1>(c) << "," << std::get<2>(c);
      if (ic == colors.size() - 1) {
        ss << "));" << std::endl;
      } else {
        ss << ")," << std::endl;
      }
    }
    str_colors = ss.str();
  }

  const std::string glsl33frag =
      str_colors +
          "uniform vec3 uniform_color;\n"
          "uniform bool use_uniform_color;\n"
          "uniform float val_min;\n"
          "uniform float val_max;\n"
          "in float val;\n"
          "out vec4 FragColor;\n"
          "void main()\n"
          "{\n"
          "  if( !use_uniform_color ) {\n"
          "    float scaled_value = (val-val_min)/(val_max-val_min) * (ncolor-1);\n"
          "    int idx_color = int(scaled_value);\n"
//          "    if( scaled_value < 0 ){ idx_color -= 1; }\n"
          "    float r01 = scaled_value - float(idx_color);\n"
          "    if( idx_color < 0 ){ idx_color = 0; r01 = 0.;}\n"
          "    if( idx_color > ncolor-2 ){ idx_color = ncolor-2; r01 = 1.; }\n"
          "    vec3 clr01 = (1.f-r01)*colors[idx_color] + r01*colors[idx_color+1];\n"
          "    FragColor = vec4(clr01.x, clr01.y, clr01.z, 1.0f);\n"
          "  }\n"
          "  else {\n"
          "    FragColor = vec4(uniform_color, 1.0f);\n"
          "  }\n"
          "}\n\0";

  // std::cout << glsl33frag << std::endl;

#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader((std::string("#version 300 es\n")+
                                      glsl33vert_projection).c_str(),
                                     (std::string("#version 300 es\n")+
                                      std::string("precision highp float;\n")+
                                      glsl33frag).c_str());
#else
  shaderProgram = dfm2::opengl::GL24_CompileShader(
      (std::string("#version 330 core\n") +
          glsl33vert_projection).c_str(),
      (std::string("#version 330 core\n") +
          glsl33frag).c_str());
#endif

  if (!glIsProgram(shaderProgram)) {
    std::cout << "shader doesnot exist" << std::endl;
  }
  glUseProgram(shaderProgram);
  Loc_MatrixProjection = glGetUniformLocation(shaderProgram, "matrixProjection");
  Loc_MatrixModelView = glGetUniformLocation(shaderProgram, "matrixModelView");
  Loc_UniformColor = glGetUniformLocation(shaderProgram, "uniform_color");
  Loc_UseUniformColor = glGetUniformLocation(shaderProgram, "use_uniform_color");
  Loc_ValMin = glGetUniformLocation(shaderProgram, "val_min");
  Loc_ValMax = glGetUniformLocation(shaderProgram, "val_max");
}

void delfem2::opengl::Drawer_MeshColormap::Draw(
    const float mat4_projection[16],
    const float mat4_modelview[16],
    float val_min,
    float val_max) const {
  glUseProgram(shaderProgram);
  glUniformMatrix4fv(
      Loc_MatrixProjection, 1, GL_FALSE,
      TransposeMat4ForOpenGL(mat4_projection, true).data());
  glUniformMatrix4fv(
      Loc_MatrixModelView, 1, GL_FALSE,
      TransposeMat4ForOpenGL(mat4_modelview, false).data());
  // glUniformMatrix4fv(Loc_MatrixProjection, 1, GL_FALSE, mP);
  // glUniformMatrix4fv(Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform1i(Loc_UseUniformColor, false);
  glUniform1f(Loc_ValMin, val_min);
  glUniform1f(Loc_ValMax, val_max);
  vao.Draw(0); // draw face
  /*
  //
  glUniform3f(Loc_UniformColor, 0, 0, 0);
  glUniform1i(Loc_UseUniformColor, true);
  vao.Draw(1); // draw line
   */
}