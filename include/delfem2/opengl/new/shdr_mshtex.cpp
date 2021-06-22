#ifdef EMSCRIPTEN
  #include <GLFW/glfw3.h>
#else
  #include "glad/glad.h" // gl3.0+
#endif

#include "delfem2/opengl/new/shdr_mshtex.h"
#include "delfem2/opengl/funcs.h"

void delfem2::opengl::CShader_MeshTex::Initialize(
    std::vector<double>& aXYZd,
    std::vector<unsigned int>& aElem,
    int gl_primitive_type,
    std::vector<double>& aTex)
{
  if( !glIsVertexArray(vao.VAO) ){ glGenVertexArrays(1, &vao.VAO); }
  vao.Delete_EBOs();
  vao.Add_EBO(aElem,gl_primitive_type);
  this->UpdateVertex(aXYZd, aElem, aTex);
}

void delfem2::opengl::CShader_MeshTex::UpdateVertex(
    std::vector<double>& aXYZd,
    std::vector<unsigned int>& aTri,
    std::vector<double>& aTex)
{
  glBindVertexArray(vao.VAO); // opengl4

  vao.ADD_VBO(0,aXYZd);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0); // gl24

  vao.ADD_VBO(1,aTex);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0); // gl24
}


void delfem2::opengl::CShader_MeshTex::Compile()
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
      "}\0";

  const std::string glsl33frag =
      "in vec2 texPrj;\n"
      "out vec4 FragColor;\n"
      "uniform sampler2D myTextureSampler;\n"
      "void main()\n"
      "{\n"
      "  FragColor = texture(myTextureSampler,texPrj);\n"
//      "  FragColor = vec4(texPrj.x, texPrj.y, 1.0, 0.0);\n"
      "}\0";

#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader((std::string("#version 300 es\n")+
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
  Loc_MatrixProjection = glGetUniformLocation(shaderProgram,  "matrixProjection");
  Loc_MatrixModelView  = glGetUniformLocation(shaderProgram,  "matrixModelView");
  Loc_Texture = glGetUniformLocation(shaderProgram, "myTextureSampler");
  std::cout << "Loc_Texture: " << Loc_Texture << std::endl;
}


void delfem2::opengl::CShader_MeshTex::Draw(
    float mP[16],
    float mMV[16]) const
{
  glUseProgram(shaderProgram);
  glUniformMatrix4fv(Loc_MatrixProjection, 1, GL_FALSE, mP);
  glUniformMatrix4fv(Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform1i(Loc_Texture, 0);
  vao.Draw(0); // draw face
}
