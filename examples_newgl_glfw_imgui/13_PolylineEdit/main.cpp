#include <iostream>
#include <vector>
#if defined(_MSC_VER)
#include <windows.h>
#endif
//
#ifdef EMSCRIPTEN
#include <emscripten/emscripten.h>
#define GLFW_INCLUDE_ES3
#define GL_GLEXT_PROTOTYPES
#define EGL_EGLEXT_PROTOTYPES
#else
#include <glad/glad.h>
#endif
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "delfem2/polyline_elastic_edit2.h"
#include "delfem2/geo_polyline2.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/new/drawer_polyline.h"

namespace dfm2 = delfem2;

// ---------------------------------------------------------------

class MyViewer : public delfem2::glfw::CViewer2 {
 public:
  void mouse_press(const float src[2]) override {
    if( mode == MODE_SKETCH) {
      vtx_xy.clear();
    }
    else if(mode == MODE_DRAG){
      const auto [ivtx,ratio] = delfem2::FindNearestPointInPolyline(
          vtx_xy,
          delfem2::CVec2f(src[0],src[1]).cast<double>());
      drag.Initialize(vtx_xy, ivtx, src[0], src[1]);
    }
  }

  void mouse_drag(
      const float src0[2],
      const float src1[2]) override {
    if (mode == MODE_SKETCH) {
      vtx_xy.emplace_back(src1[0], src1[1]);
    } else if (mode == MODE_SMOOTHING) {
      const auto[ivtx, ratio] = delfem2::FindNearestPointInPolyline(
          vtx_xy,
          delfem2::CVec2f(src0[0], src0[1]).cast<double>());
      unsigned int nvtx = vtx_xy.size();
      for (int jvtx = (int) ivtx - 2; jvtx < (int) ivtx + 2; ++jvtx) {
        if (jvtx <= 0 || jvtx >= (int) nvtx - 1) { continue; }
        unsigned int jvtx0 = jvtx - 1;
        unsigned int jvtx1 = jvtx;
        unsigned int jvtx2 = jvtx + 1;
        vtx_xy[jvtx1] = (vtx_xy[jvtx0] + vtx_xy[jvtx2]) * 0.5;
      }
    } else if (mode == MODE_DRAG) {
      drag.Drag(vtx_xy, src1[0], src1[1]);
    }
  }

  void mouse_release() override {
    if( mode == MODE_SKETCH) {
      vtx_xy = delfem2::Polyline_Resample_Polyline(vtx_xy, 0.05);
    }
  }

  void InitGL(){
    drawer_polyline.InitGL();
    drawer_polyline.radius_cylinder = 0.01;
    drawer_polyline.radius_sphere = 0.02;
    drawer_polyline.sphere.color = {1,0,0,1};
    drawer_polyline.cylinder.color = {0,0,0,1};
  }

  void Draw()  {
    ::glEnable(GL_DEPTH_TEST);
    ::glEnable(GL_TEXTURE_2D);
    drawer_polyline.Draw(
        vtx_xy,2,
        GetProjectionMatrix().data(),
        GetModelViewMatrix().data());
  }
 public:
  std::vector<dfm2::CVec2d> vtx_xy;
  delfem2::opengl::Drawer_Polyline drawer_polyline;
  int mode = 0;
  enum {
    MODE_SKETCH,
    MODE_SMOOTHING,
    MODE_DRAG,
  };
  delfem2::PolylineElasticEdit2 drag;
} viewer;

void draw(GLFWwindow *window) {
  glfwPollEvents();
  ::glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);

  // feed inputs to dear imgui, start new frame
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  { // render your GUI
    ImGui::Begin("Triangle Mesh");
    ImGui::RadioButton("sketch", &(viewer.mode), MyViewer::MODE_SKETCH);
    ImGui::RadioButton("smoothing", &(viewer.mode), MyViewer::MODE_SMOOTHING);
    ImGui::RadioButton("drag", &(viewer.mode), MyViewer::MODE_DRAG);
    ImGui::End();
  }

  viewer.Draw();

  // Render dear imgui into screen
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);

  glfwSwapBuffers(window);
}

int main(int, char **) {

  delfem2::glfw::InitGLNew();
  viewer.OpenWindow();
  
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif
  viewer.InitGL();

  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui_ImplGlfw_InitForOpenGL(viewer.window, true);// Setup Platform/Renderer bindings
  ImGui::StyleColorsClassic(); // Setup Dear ImGui style

#ifdef EMSCRIPTEN
  ImGui_ImplOpenGL3_Init("#version 100");
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, window, 60, 1);
#else
  ImGui_ImplOpenGL3_Init("#version 150");
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif

  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  glfwDestroyWindow(viewer.window);
  glfwTerminate();

  return 0;
}
