/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif

#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#else
#  include <glad/glad.h>
#endif

#include "delfem2/polyline_elastic_rod2.h"
#include "delfem2/vec2.h"
#include "delfem2/mshuni.h"
#include "delfem2/geo_polyline.h"
#include "delfem2/colormap.h"
#include "delfem2/glfw/util.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/opengl/new/drawer_polyline.h"
#include "delfem2/opengl/new/drawer_polylinecolormap.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

class MyViewer : public delfem2::glfw::CViewer2 {
 public:
  void mouse_press(const float src[2]) override {
    if (this->nav.imodifier == GLFW_MOD_SHIFT) { return; }
    if (this->nav.ibutton == GLFW_MOUSE_BUTTON_LEFT) {  // sketch
      vtx_xy.clear();
      vtx_xy.emplace_back(src[0], src[1]);
    } else if (this->nav.ibutton == GLFW_MOUSE_BUTTON_MIDDLE) {  // smooth
    } else if (this->nav.ibutton == GLFW_MOUSE_BUTTON_RIGHT) {  // drag
      const delfem2::CVec2f s0(src[0], src[1]);
      double t = delfem2::Nearest_Polyline(vtx_xy, s0);
      dfm2::CVec2f s1 = delfem2::Sample_Polyline(vtx_xy, t);
      if ((s0 - s1).norm() > 0.1) {
        ivtx = UINT_MAX;
        return;
      }
      ivtx = int(t);
      ls_.Initialize(vtx_xy.size());
    }
  }

  void mouse_drag(const float src0[2], const float src1[2]) override {
    if (this->nav.ibutton == GLFW_MOUSE_BUTTON_LEFT) {
      vtx_xy.emplace_back(src1[0], src1[1]);
    } else if (this->nav.ibutton == GLFW_MOUSE_BUTTON_MIDDLE) {  // smooth
      const delfem2::CVec2f s0(src1[0], src1[1]);
      double t = delfem2::Nearest_Polyline(vtx_xy, s0);
      dfm2::CVec2f s1 = delfem2::Sample_Polyline(vtx_xy, t);
      if ((s0 - s1).norm() > 0.1) {
        ivtx = UINT_MAX;
        return;
      }
      ivtx = int(t);
      delfem2::Smooth_Polyline(vtx_xy, ivtx, 0.8);
    } else if (this->nav.ibutton == GLFW_MOUSE_BUTTON_RIGHT) {  // drag
      if (ivtx >= vtx_xy.size()) { return; }
      DragPolylineElastic_Rod2(
          vtx_xy, ls_,
          ivtx, {src1[0], src1[1]},
          edge_length, 0.001);
    }
  }

  void mouse_release() override {
    if (this->nav.ibutton == GLFW_MOUSE_BUTTON_LEFT) {
      vtx_xy = delfem2::Polyline_Resample_Polyline(
          vtx_xy,
          static_cast<float>(edge_length));
      if (vtx_xy.size() < 4) { vtx_xy.clear(); }
    } else if (this->nav.ibutton == GLFW_MOUSE_BUTTON_MIDDLE) {  // smooth
    } else if (this->nav.ibutton == GLFW_MOUSE_BUTTON_RIGHT) {  // drag
    }
  }

  void InitGL() {
    drawer_polyline_.InitGL();
    //
    //constexpr auto& colormap = delfem2::colormap_magma<double>;
    //constexpr auto& colormap = delfem2::colormap_parula<double>;
    //constexpr auto& colormap = delfem2::colormap_turbo<double>;
    //constexpr auto& colormap = delfem2::colormap_inferno<double>;
    //constexpr auto& colormap = delfem2::colormap_plasma<double>;
    constexpr auto& colormap = delfem2::colormap_viridis<double>;
    //constexpr auto& colormap = delfem2::colormap_cividis<double>;
    //constexpr auto& colormap = delfem2::colormap_github<double>;
    //constexpr auto& colormap = delfem2::colormap_hot<double>;
    unsigned int ncolor = sizeof(colormap)/sizeof(colormap[0]);
    drawer_polylinecolormap_.sphere.colors.assign(colormap, colormap+ncolor);
    drawer_polylinecolormap_.cylinder.colors.assign(colormap, colormap+ncolor);
    drawer_polylinecolormap_.InitGL();
  }

  void Draw() {
    drawer_polyline_.radius_cylinder = 0.01;
    drawer_polyline_.radius_sphere = 0.01;
    drawer_polyline_.Draw(
        vtx_xy, 2,
        this->GetProjectionMatrix().data(),
        this->GetModelViewMatrix().data());
  }

  void DrawColormap() {
    const unsigned int nvtx = vtx_xy.size();
    std::vector<float> vtx_val(nvtx);
    for(unsigned int iv=0;iv<nvtx;++iv){
      vtx_val[iv] = static_cast<float>(iv);
    }

    drawer_polylinecolormap_.radius_cylinder = 0.01;
    drawer_polylinecolormap_.radius_sphere = 0.01;
    drawer_polylinecolormap_.Draw(
        vtx_xy, 2,
        this->GetProjectionMatrix().data(),
        this->GetModelViewMatrix().data(),
        vtx_val, 0, nvtx);
  }

 public:
  delfem2::opengl::Drawer_Polyline drawer_polyline_;
  delfem2::opengl::Drawer_PolylineColormap drawer_polylinecolormap_;

  const double edge_length = 0.05;
  std::vector<delfem2::CVec2f> vtx_xy;
  unsigned int ivtx;
  delfem2::LinearSystemSolver_BlockPentaDiagonal<2> ls_;
};

int main() {

  delfem2::glfw::InitGLNew();
  MyViewer viewer;
  viewer.view_height = 0.5;
  viewer.trans[0] = -0.5;
  viewer.OpenWindow();
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif
  viewer.InitGL();

  while (true) {
    ::glClearColor(0.8, 1.0, 1.0, 1.0);
    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //viewer.Draw();
    viewer.DrawColormap();
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  viewer.ExitIfClosed();
  return 0;
}


