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
#include "delfem2/mshio.h"
#include "delfem2/primitive.h"

#include "delfem2/gl24_funcs.h"
#include "delfem2/gl4_funcs.h"
#include "../glfw_funcs.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

CNav3D_GLFW nav;
int shaderProgram;
int projectionMatrixLoc;
CGL4_VAO_Mesh vao_mesh;

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
  glClear(GL_COLOR_BUFFER_BIT);

  glUseProgram(shaderProgram);
  float mP[16]; nav.camera.Affine4f_Projection(mP, asp, 10);
  glUniformMatrix4fv(projectionMatrixLoc, 1, GL_FALSE, mP);
  vao_mesh.Draw();
  
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


int main(void)
{
  GLFWwindow* window = myGLFW_OpenWindow(800,600);
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, callback_resize);
  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, callback_key);
  
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
    MeshTri3D_Torus(aXYZd, aTri, 1.0, 0.2);
    std::cout << "nXYZ: " << aXYZd.size()/3 << "  nTri: "  << aTri.size()/3 << std::endl;
//    Read_Obj(std::string(PATH_INPUT_DIR)+"/bunny_1k.obj", aXYZ, aTri);
    Normalize(aXYZd);
    std::vector<float> aXYZf(aXYZd.begin(),aXYZd.end());
    int VAO = GL4_VAO_MeshTri3D(aXYZf.data(),aXYZf.size()/3,3,
                                aTri.data(),aTri.size()/3);
    int nTri = aTri.size()/3;
    std::cout << "VAO: " << VAO << std::endl;
    vao_mesh.VAO = VAO;
    vao_mesh.nTri = nTri;
  }

#ifdef EMSCRIPTEN
  shaderProgram = GL24_CompileShader(glsles3vert_projection.c_str(),
                                     glsles3frag.c_str());
#else
  shaderProgram = GL24_CompileShader(glsl33vert_projection.c_str(),
                                     glsl33frag.c_str());
#endif
  
  
  {
    if( !glIsProgram(shaderProgram) ){
      std::cout << "shader doesnot exist" << std::endl;
    }
    glUseProgram(shaderProgram);
    projectionMatrixLoc = glGetUniformLocation(shaderProgram,  "projectionMatrix");
    std::cout << "projectionMatrixLoc: " << projectionMatrixLoc << "   shaderProgram: " << shaderProgram << std::endl;
  }
  
  nav.camera.view_height = 1.5;
  
  
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, window, 60, 1);
#else
  while (!glfwWindowShouldClose(window)) { draw(window); }
#endif
  
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

