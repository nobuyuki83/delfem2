#include <iostream>
#include <cmath>
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

#include "delfem2/opengl/gl_funcs.h"
#include "delfem2/opengl/glnew_funcs.h"
#include "delfem2/opengl/glnew_mshcolor.h"
//
#include "delfem2/opengl/glfw/cam_glfw.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// ---------------------------
dfm2::opengl::CShader_TriMesh shdr;
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
  
  viewer.DrawEnd_oldGL();
}

int main()
{
  viewer.Init_newGL();
  
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
    delfem2::MeshTri3_Torus(aXYZd, aTri,
                            1.0, 0.2,
                            32,18);
    shdr.Initialize(aXYZd, aTri);
  }
 
  viewer.nav.camera.view_height = 2.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

