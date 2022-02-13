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
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/polyline_elastic_rod2.h"
#include "delfem2/vec2.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/geo_polyline.h"
#include "delfem2/glfw/util.h"
#include "delfem2/glfw/viewer2.h"

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
      ls.Initialize(vtx_xy.size());
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
          vtx_xy, ls,
          ivtx, {src1[0], src1[1]},
          edge_length, 0.001);
    }
  }

  void mouse_release() override {
    if (this->nav.ibutton == GLFW_MOUSE_BUTTON_LEFT) {
      vtx_xy = delfem2::Polyline_Resample_Polyline(vtx_xy, edge_length);
      if (vtx_xy.size() < 4) { vtx_xy.clear(); }
    } else if (this->nav.ibutton == GLFW_MOUSE_BUTTON_MIDDLE) {  // smooth
    } else if (this->nav.ibutton == GLFW_MOUSE_BUTTON_RIGHT) {  // drag
    }
  }

  void Draw() {
    ::glColor3d(1, 0, 0);
    ::glLineWidth(1);
    ::glBegin(GL_LINE_STRIP);
    for (auto &ivtx: vtx_xy) {
      ::glVertex2fv(ivtx.data());
    }
    ::glEnd();
    //
    glPointSize(5);
    ::glBegin(GL_POINTS);
    for (auto &ivtx: vtx_xy) {
      ::glVertex2fv(ivtx.data());
    }
    ::glEnd();
  }

 public:
  const double edge_length = 0.05;
  std::vector<delfem2::CVec2f> vtx_xy;
  unsigned int ivtx;
  delfem2::LinearSystemSolver_BlockPentaDiagonal<2> ls;
};

int main() {
  dfm2::glfw::InitGLOld();
  MyViewer viewer;
  viewer.view_height = 0.5;
  viewer.trans[0] = -0.5;
  viewer.OpenWindow();

  while (true) {
    viewer.DrawBegin_oldGL();
    viewer.Draw();
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  viewer.ExitIfClosed();
  return 0;
}


