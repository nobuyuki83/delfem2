#include <iostream>
#include <math.h>

#if defined(_MSC_VER)
  #include <windows.h>
#endif

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#endif

#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/primitive.h"

#include "delfem2/gl24_funcs.h"
#include "delfem2/gl4_funcs.h"
#include "../glfw_funcs.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

CNav3D_GLFW nav;
int shaderProgram;
int Loc_MatrixProjection;
int Loc_MatrixModelView;
int Loc_Color;
CGL4_VAO_Mesh vao_face;
CGL4_VAO_Mesh vao_edge;

void draw(GLFWwindow* window)
{
  float asp;
  {
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);
    asp = width / (float) height;
//    std::cout << width << " " << height << " " << asp << std::endl;
  }
  
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );

  glUseProgram(shaderProgram);
  float mP[16]; nav.camera.Affine4f_Projection(mP, asp, 10);
  glUniformMatrix4fv(Loc_MatrixProjection, 1, GL_FALSE, mP);
  float mMV[16]; nav.camera.Affine4f_ModelView(mMV);
  glUniformMatrix4fv(Loc_MatrixModelView, 1, GL_FALSE, mMV);
  glUniform3f(Loc_Color, 1,0,0);
  vao_face.Draw();
  glUniform3f(Loc_Color, 0,0,0);
  vao_edge.Draw();
  
  
  
  glfwSwapBuffers(window);
  glfwPollEvents();
}

void callback_key(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS){
    glfwSetWindowShouldClose(window, GL_TRUE);
  }
}

void callback_resize(GLFWwindow* window, int width, int height)
{
  glViewport(0, 0, width, height);
}

void callback_mouse_button(GLFWwindow* window, int button, int action, int mods)
{
  nav.Mouse(window,button,action,mods);

}

void callback_cursor_position(GLFWwindow* window, double xpos, double ypos)
{
  nav.Motion(window,xpos,ypos);
}

void callback_scroll(GLFWwindow* window, double xoffset, double yoffset)
{
  nav.camera.scale *= pow(1.01,yoffset);
}


int main(void)
{
  GLFWwindow* window = myGLFW_OpenWindow(800,600);
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, callback_resize);
  glfwSetKeyCallback(            window, callback_key);
  glfwSetMouseButtonCallback(    window, callback_mouse_button);
  glfwSetCursorPosCallback(      window, callback_cursor_position);
  glfwSetScrollCallback(         window, callback_scroll);
  
  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  {
    std::vector<double> aXYZd;
    std::vector<unsigned int> aTri;
    MeshTri3D_Torus(aXYZd, aTri,
                    1.0, 0.2);
    std::vector<unsigned int> aLine;
    MeshLine_MeshElem(aLine,
                      aTri.data(), aTri.size()/3, MESHELEM_TRI,
                      aXYZd.size()/3);
    std::cout << "nXYZ: " << aXYZd.size()/3 << "  nTri: "  << aTri.size()/3 << std::endl;
    Normalize(aXYZd);
    std::vector<double> aNrmd(aXYZd.size());
    Normal_MeshTri3D(aNrmd.data(),
                     aXYZd.data(), aXYZd.size()/3,
                     aTri.data(), aTri.size()/3);
    std::vector<float> aXYZf(aXYZd.begin(),aXYZd.end());
    std::vector<float> aNrmf(aNrmd.begin(),aNrmd.end());
    ///
    {
      int VAO = GL4_VAO_MeshTri3D_FaceNormal(aXYZf.data(),aXYZf.size()/3,3,
                                             aTri.data(),aTri.size()/3,
                                             aNrmf.data());
      int nTri = aTri.size()/3;
      std::cout << "VAO: " << VAO << std::endl;
      vao_face.VAO = VAO;
      vao_face.nElem = nTri;
      vao_face.nNoel = 3;
      vao_face.GL_ELEM_TYPE = GL_TRIANGLES;
    }
    {
      int VAO = GL4_VAO_MeshTri3D(aXYZf.data(),aXYZf.size()/3,3,
                                  aLine.data(),aLine.size()/2);
      std::cout << "VAO: " << VAO << std::endl;
      vao_edge.VAO = VAO;
      vao_edge.nElem = aLine.size()/2;
      vao_edge.nNoel = 2;
      vao_edge.GL_ELEM_TYPE = GL_LINES;
    }
  }

#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader((std::string("#version 300 es\n")+
                                      glsl33vert_projection).c_str(),
                                     (std::string("#version 300 es\n")+
                                      std::string("precision highp float;\n")+
                                      glsl33frag).c_str());
#else
  shaderProgram = GL24_CompileShader((std::string("#version 330 core\n")+
                                      glsl33vert_projection).c_str(),
                                     (std::string("#version 330 core\n")+
                                      glsl33frag).c_str());
#endif
  
  
  {
    if( !glIsProgram(shaderProgram) ){
      std::cout << "shader doesnot exist" << std::endl;
    }
    glUseProgram(shaderProgram);
    Loc_MatrixProjection = glGetUniformLocation(shaderProgram,  "matrixProjection");
    Loc_MatrixModelView  = glGetUniformLocation(shaderProgram,  "matrixModelView");
    Loc_Color            = glGetUniformLocation(shaderProgram,  "color");
    std::cout << "projectionMatrixLoc: " << Loc_MatrixProjection << "   shaderProgram: " << shaderProgram << "  LocColor: " << Loc_Color << std::endl;
  }
  
  nav.camera.view_height = 1.5;
  nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, window, 60, 1);
#else
  while (!glfwWindowShouldClose(window)) { draw(window); }
#endif
  
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

