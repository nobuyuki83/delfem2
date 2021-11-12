
#if defined(_MSC_VER)
#include <windows.h>
#endif

#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#else
#  include <glad/glad.h>
#endif

#include "delfem2/cad2_dtri2.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/new/funcs.h"
#include "delfem2/opengl/new/drawer_cad2.h"
#include <GLFW/glfw3.h>
#include <iostream>

// end of header
// -----------------------------------------------------

class CCAD2D_Viewer : public delfem2::glfw::CViewer2 {
 public:
  CCAD2D_Viewer() {
    std::vector<double> aXY = {-1, -1, +1, -1, +1, +1, -1, +1};
    cad.AddPolygon(aXY);
    cad.aEdge[0].SetCubicBezierCurve({-0.8,-0.7}, {+0.8,-0.7});
    {
      mesher.edge_length = -1;
      mesher.Meshing(dmesh, cad);
    }
  }
  void mouse_press(const float src[2]) override {
    cad.Pick(src[0], src[1], 1.0);
  }
  void mouse_drag(const float src0[2], const float src1[2]) override {
    cad.DragPicked(src1[0], src1[1], src0[0], src0[1]);
    {
      mesher.edge_length = -1;
      mesher.Meshing(dmesh, cad);
    }
    shdr_cad.MakeBuffer(cad,dmesh);
  }
 public:
  delfem2::CCad2D cad;
  delfem2::opengl::CShader_Cad2D shdr_cad;
  delfem2::CMeshDynTri2D dmesh;
  delfem2::CMesher_Cad2D mesher;
};

CCAD2D_Viewer viewer;

// -----------------------------------------------------

void draw(GLFWwindow *window) {
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);

  delfem2::CMat4f mMV = viewer.GetModelViewMatrix();
  delfem2::CMat4f mP = viewer.GetProjectionMatrix();
  const delfem2::CMat4f mZ = delfem2::CMat4f::ScaleXYZ(1, 1, -1);
  viewer.shdr_cad.Draw(
      (mZ * mP).transpose().data(),
      mMV.transpose().data(),
      viewer.cad);
  glfwSwapBuffers(window);
  glfwPollEvents();
}

int main() {
  delfem2::glfw::InitGLNew();
  viewer.view_height = 2.0;
  viewer.OpenWindow();
  
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif
  viewer.shdr_cad.Compile();
  viewer.shdr_cad.MakeBuffer(viewer.cad, viewer.dmesh);
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

