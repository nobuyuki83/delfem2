#include "delfem2/opengl/new/drawer_mshtex.h"

#ifdef EMSCRIPTEN
  #include <GLFW/glfw3.h>
#else
  #include "glad/glad.h" // gl3.0+
#endif

#include "delfem2/opengl/funcs.h"

void delfem2::opengl::Drawer_MeshTex::SetElement(
    std::vector<unsigned int>& elem_vtx,
    int gl_primitive_type)
{
  if( !glIsVertexArray(vao.VAO) ){ glGenVertexArrays(1, &vao.VAO); }
  vao.Delete_EBOs();
  vao.Add_EBO(elem_vtx, gl_primitive_type);
}

template <typename REAL>
void delfem2::opengl::Drawer_MeshTex::SetCoords(
    std::vector<REAL>& vtx_coords,
    unsigned int ndim)
{
  if (!glIsVertexArray(vao.VAO)) { glGenVertexArrays(1, &vao.VAO); }
  vao.ADD_VBO(0, vtx_coords);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(
      0,
      ndim, convertToGlType<REAL>(), GL_FALSE,
      ndim * sizeof(REAL), (void *) nullptr); // gl24
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::Drawer_MeshTex::SetCoords(
    std::vector<float>& aXYZd,
    unsigned int ndim);
template void delfem2::opengl::Drawer_MeshTex::SetCoords(
    std::vector<double>& aXYZd,
    unsigned int ndim);
#endif


template <typename REAL>
void delfem2::opengl::Drawer_MeshTex::SetTexUV(
    std::vector<REAL>& aTex)
{
  if (!glIsVertexArray(vao.VAO)) { glGenVertexArrays(1, &vao.VAO); }
  vao.ADD_VBO(1,aTex);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(
      1,
      2, convertToGlType<REAL>(), GL_FALSE,
      2*sizeof(REAL), (void*)nullptr); // gl24
}
#ifdef DFM2_STATIC_LIBRARY
template void delfem2::opengl::Drawer_MeshTex::SetTexUV(
    std::vector<float>& aTex);
template void delfem2::opengl::Drawer_MeshTex::SetTexUV(
    std::vector<double>& aTex);
#endif

void delfem2::opengl::Drawer_MeshTex::InitGL()
{
  const std::string glsl33vert_projection =
      "uniform mat4 matrixProjection;\n"
      "uniform mat4 matrixModelView;\n"
      "layout (location = 0) in vec3 posIn;\n"
      "layout (location = 1) in vec2 texIn;\n"
      "out vec2 texPrj;\n"
      "void main()\n"
      "{\n"
      "  gl_Position = matrixProjection * matrixModelView * vec4(posIn.x, posIn.y, posIn.z, 1.0);\n"
      "  texPrj = texIn;\n"
      "}";

  const std::string glsl33frag =
      "in vec2 texPrj;\n"
      "out vec4 FragColor;\n"
      "uniform sampler2D myTextureSampler;\n"
      "void main()\n"
      "{\n"
      "  FragColor = texture(myTextureSampler,texPrj);\n"
//      "  FragColor = vec4(texPrj.x, texPrj.y, 1.0, 0.0);\n"
      "}";

#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader(
      (std::string("#version 300 es\n")+
      glsl33vert_projection).c_str(),
      (std::string("#version 300 es\n")+
      std::string("precision highp float;\n")+
      glsl33frag).c_str());
#else
  shaderProgram = ::delfem2::opengl::GL24_CompileShader(
      (std::string("#version 330 core\n")+
      glsl33vert_projection).c_str(),
      (std::string("#version 330 core\n")+
      glsl33frag).c_str());
#endif

  if( !glIsProgram(shaderProgram) ){
    std::cerr << "shader doesnot exist" << std::endl;
  }
  glUseProgram(shaderProgram);
  loc_projection_matrix = glGetUniformLocation(shaderProgram,  "matrixProjection");
  loc_modelview_matrix  = glGetUniformLocation(shaderProgram,  "matrixModelView");
  loc_texture = glGetUniformLocation(shaderProgram, "myTextureSampler");
  // std::cout << "Loc_Texture: " << loc_texture << std::endl;
}


void delfem2::opengl::Drawer_MeshTex::Draw(
    const float mat4_projection[16],
    const float mat4_modelview[16]) const
{
  glUseProgram(shaderProgram);
  glUniformMatrix4fv(
      loc_projection_matrix, 1, GL_FALSE,
      TransposeMat4ForOpenGL(mat4_projection, true).data());
  glUniformMatrix4fv(
      loc_modelview_matrix, 1, GL_FALSE,
      TransposeMat4ForOpenGL(mat4_modelview, false).data());
  glUniform1i(loc_texture, 0);
  vao.Draw(0); // draw face
}

// ====================================================

void delfem2::opengl::Drawer_MeshMixTwoTextures::InitGL() {

  const std::string glsl33vert_projection =
      "uniform mat4 matrixProjection;\n"
      "uniform mat4 matrixModelView;\n"
      "layout (location = 0) in vec3 posIn;\n"
      "layout (location = 1) in vec2 texIn;\n"
      "out vec2 texPrj;\n"
      "void main()\n"
      "{\n"
      "  gl_Position = matrixProjection * matrixModelView * vec4(posIn.x, posIn.y, posIn.z, 1.0);\n"
      "  texPrj = texIn;\n"
      "}";

  const std::string glsl33frag =
      "in vec2 texPrj;\n"
      "out vec4 FragColor;\n"
      "uniform float ratio;\n"
      "uniform sampler2D tex0;\n"
      "uniform sampler2D tex1;\n"
      "void main()\n"
      "{\n"
      "  FragColor = (1-ratio)*texture(tex0,texPrj) + ratio*texture(tex1,texPrj);\n"
      "}";
#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader(
      (std::string("#version 300 es\n")+
      glsl33vert_projection).c_str(),
      (std::string("#version 300 es\n")+
      std::string("precision highp float;\n")+
      glsl33frag).c_str());
#else
  shaderProgram = ::delfem2::opengl::GL24_CompileShader(
      (std::string("#version 330 core\n")+
          glsl33vert_projection).c_str(),
      (std::string("#version 330 core\n")+
          glsl33frag).c_str());
#endif

  if( !glIsProgram(shaderProgram) ){
    std::cout << "shader doesnot exist" << std::endl;
  }
  glUseProgram(shaderProgram);
  loc_projection_matrix = glGetUniformLocation(shaderProgram,  "matrixProjection");
  loc_modelview_matrix  = glGetUniformLocation(shaderProgram,  "matrixModelView");
  loc_ratio = glGetUniformLocation(shaderProgram, "ratio");
  this->loc_texture = glGetUniformLocation(shaderProgram, "tex0");
  this->loc_overlaytex = glGetUniformLocation(shaderProgram, "tex1");
  // std::cout << "Loc_Texture: " << loc_Texture << " " << loc_tex1 << std::endl;
}

void delfem2::opengl::Drawer_MeshMixTwoTextures::Draw(
    const float mat4_projection[16],
    const float mat4_modelview[16],
    float ratio) const
{
  glUseProgram(shaderProgram);
  glUniformMatrix4fv(
      loc_projection_matrix, 1, GL_FALSE,
      TransposeMat4ForOpenGL(mat4_projection, true).data());
  glUniformMatrix4fv(
      loc_modelview_matrix, 1, GL_FALSE,
      TransposeMat4ForOpenGL(mat4_modelview, false).data());
  glUniform1i(loc_texture, 0);
  glUniform1i(loc_overlaytex, 1);
  glUniform1f(loc_ratio, ratio);
  vao.Draw(0); // draw face
}
