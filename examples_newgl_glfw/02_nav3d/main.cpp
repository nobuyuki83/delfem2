
#include <iostream>
#include <cmath>

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#else
  #include <glad/glad.h>
#endif

#if defined(_MSC_VER)
  #include <windows.h>
#endif
#include <GLFW/glfw3.h>

#include "delfem2/msh_primitive.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/new/drawer_mshtri.h"

namespace dfm2 = delfem2;

// ---------------------------

struct MyData {
  dfm2::opengl::CShader_TriMesh shdr;
  delfem2::glfw::CViewer3 viewer;
};

// ---------------------------

void draw(MyData* data)
{
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  data->shdr.Draw(
      data->viewer.GetProjectionMatrix().data(),
      data->viewer.GetModelViewMatrix().data());
  data->viewer.SwapBuffers();
  glfwPollEvents();
}

int main()
{
  MyData data;
  data.viewer.projection = std::make_unique<delfem2::Projection_LookOriginFromZplus>(2, false);

  dfm2::glfw::InitGLNew();
  data.viewer.OpenWindow();
  
  // glad: load all OpenGL function pointers
  // ---------------------------------------
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif

  data.shdr.InitGL();
  {
    std::vector<float> aXYZd;
    std::vector<unsigned int> aTri;
    delfem2::MeshTri3_Torus(
        aXYZd, aTri,
        1.f, 0.2f, 32,18);
    data.shdr.Initialize(aXYZd, 3, aTri);
  }
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, &data, 0, 1);
#else
  while (!glfwWindowShouldClose(data.viewer.window)) { draw(&data); }
#endif
  
  glfwDestroyWindow(data.viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

