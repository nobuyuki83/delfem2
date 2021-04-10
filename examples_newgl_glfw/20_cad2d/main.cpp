
#if defined(_MSC_VER)
  #include <windows.h>
#endif

#ifdef EMSCRIPTEN
  #include <emscripten/emscripten.h>
  #define GLFW_INCLUDE_ES3
#else
  #include <glad/glad.h>
#endif

#include "delfem2/glfw/viewer2.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/new/funcs.h"
#include "delfem2/opengl/new/v23dtricad.h"
#include "delfem2/cad2_dtri2.h"
#include <GLFW/glfw3.h>
#include <iostream>

// end of header
// -----------------------------------------------------

class CCAD2D_Viewer : public delfem2::glfw::CViewer2
{
public:
  CCAD2D_Viewer(){
    std::vector<double> aXY = {-1,-1, +1,-1, +1,+1, -1,+1};
    cad.AddPolygon(aXY);
    {
      double param[4] = {0.2, 0.3, 0.8, +0.3};
      std::vector<double> vparam(param,param+4);
      cad.SetEdgeType( 0, delfem2::CCad2D_EdgeGeo::BEZIER_CUBIC, vparam );
    }
  }
  void mouse_press(const float src[2]) override {
    cad.Pick(src[0], src[1], 1.0);
  }
  void mouse_drag(const float src0[2], const float src1[2]) override {
    cad.DragPicked(src1[0],src1[1], src0[0],src0[1]);
    shdr_cad.MakeBuffer(cad);
  }
public:
  delfem2::CCad2D cad;
  delfem2::opengl::CShader_Cad2D shdr_cad;
};

CCAD2D_Viewer viewer;

// -----------------------------------------------------

void draw(GLFWwindow* window)
{
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  int nw, nh; glfwGetFramebufferSize(window, &nw, &nh);
  const float asp = (float)nw/nh;
  float mMV[16], mP[16]; viewer.Mat4_MVP_OpenGL(mMV, mP, asp);
  viewer.shdr_cad.Draw(mP, mMV, viewer.cad);
  
  glfwSwapBuffers(window);
  glfwPollEvents();
}

int main()
{
  delfem2::glfw::InitGLNew();
  viewer.view_height = 2.0;
  viewer.InitGL();
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif
  viewer.shdr_cad.Compile();
  viewer.shdr_cad.MakeBuffer(viewer.cad);
  
#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

