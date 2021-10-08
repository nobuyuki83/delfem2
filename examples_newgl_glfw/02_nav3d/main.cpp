
#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#else
  #include <glad/glad.h>
#endif

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/new/shdr_mshtri.h"
#include "delfem2/mshprimitive.h"

#if defined(_MSC_VER)
  #include <windows.h>
#endif
#include <GLFW/glfw3.h>

#include <iostream>
#include <cmath>

namespace dfm2 = delfem2;

// ---------------------------
// global variables
dfm2::opengl::CShader_TriMesh shdr;
delfem2::glfw::CViewer3 viewer(2.0);

// ---------------------------

void draw()
{
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  shdr.Draw(viewer.GetProjectionMatrix().data(),
            viewer.GetModelViewMatrix().data());
  viewer.SwapBuffers();
  glfwPollEvents();
}

int main()
{
  dfm2::glfw::InitGLNew();
  viewer.InitGL();
  
  // glad: load all OpenGL function pointers
  // ---------------------------------------
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif

  shdr.Compile();
  {
    std::vector<double> aXYZd;
    std::vector<unsigned int> aTri;
    delfem2::MeshTri3_Torus(aXYZd, aTri,
                            1.0, 0.2,
                            32,18);
    shdr.Initialize(aXYZd, 3, aTri);
  }
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(); }
#endif
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

