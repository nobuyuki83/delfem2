#include <iostream>
#include <math.h>
#include "delfem2/primitive.h"

// ---

#if defined(_MSC_VER)
  #include <windows.h>
#endif

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#endif

#include "delfem2/opengl/gl24_funcs.h"
#include "delfem2/opengl/gl4_funcs.h"
#include "delfem2/opengl/gl4_mshcolor.h"
#include "delfem2/opengl/glfw_cam.h"
#include "delfem2/opengl/glfw_viewer.hpp"

// ---------------------------
CShader_TriMesh shdr;
delfem2::opengl::CViewer_GLFW viewer;
// ---------------------------


void draw(GLFWwindow* window)
{
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  

  float mP[16], mMV[16];
  viewer.nav.Matrix_MVP(mMV, mP, window);
  shdr.Draw(mP, mMV);

  glfwSwapBuffers(window);
  glfwPollEvents();
}

int main(void)
{
  viewer.Init_GLnew();
  
  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  shdr.Compile();
  {
    std::vector<double> aXYZd;
    std::vector<unsigned int> aTri;
    delfem2::MeshTri3D_Torus(aXYZd, aTri,
                             1.0, 0.2);
    shdr.Initialize(aXYZd, aTri);
  }
 
  viewer.nav.camera.view_height = 2.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

