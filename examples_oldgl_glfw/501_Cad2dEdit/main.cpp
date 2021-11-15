/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <stack>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cad2.h"
#include "delfem2/cad2_mesher.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/cad2dtriv2.h"

#ifndef M_PI
#  define M_PI 3.141592653589793
#endif
// -------------------------------------

class Cad2Viewer : public delfem2::glfw::CViewer3 {
 public:
  explicit Cad2Viewer(delfem2::CCad2D& cad) : CViewer3(1.5), cad(cad) {
    delfem2::CMesher_Cad2D mesher;
    mesher.edge_length = -1;
    mesher.Meshing(dmesh, cad);
  }
  //
  void mouse_press(
      const float src[3],
      [[maybe_unused]] const float dir[3]) override {
    cad_ui.Pick(src[0], src[1], 1.5, cad);
  }

  void mouse_drag(
      const float src0[3],
      const float src1[3],
      [[maybe_unused]] const float dir[3]) override {
    cad_ui.DragPicked(cad, src1[0], src1[1], src0[0], src0[1]);
    delfem2::CMesher_Cad2D mesher;
    mesher.edge_length = -1;
    mesher.Meshing(dmesh, cad);
  }
  //
  void Draw() const {
    this->DrawBegin_oldGL();
    delfem2::opengl::DrawCad2Vtxs(cad, cad_ui.ivtx_picked);
    delfem2::opengl::DrawCad2Edges(cad, cad_ui.iedge_picked);
    delfem2::opengl::DrawMeshDynTri_Edge(dmesh.aETri,dmesh.aVec2);
    ::glColor3d(0.8, 0.8, 0.8);
    delfem2::opengl::DrawMeshDynTri_FaceNorm(dmesh.aETri, dmesh.aVec2);
    this->SwapBuffers();
  }
 public:
  delfem2::CCad2D &cad;
  delfem2::Cad2_Ui cad_ui;
  delfem2::CMeshDynTri2D dmesh;
};

int main() {
  delfem2::glfw::InitGLOld();
  {
    delfem2::CCad2D cad;
    cad.AddPolygon(std::vector<double>{-1, -1, +1, -1, +1, +1, -1, +1});
    Cad2Viewer viewer(cad);
    viewer.OpenWindow();
    while (!glfwWindowShouldClose(viewer.window)) {
      viewer.Draw();
      glfwPollEvents();
    }
    glfwDestroyWindow(viewer.window);
  }
  {
    delfem2::CCad2D cad;
    cad.AddPolygon(std::vector<double>{-1, -1, +1, -1, +1, +1, -1, +1});
    cad.aEdge[0].SetCubicBezierCurve(
        {-0.8, -1.2},
        {+0.8,-1.2});
    cad.aEdge[2].SetQuadraticBezierCurve(
        {0.0, +1.2});
    Cad2Viewer viewer(cad);
    viewer.OpenWindow();
    while (!glfwWindowShouldClose(viewer.window)) {
      viewer.Draw();
      glfwPollEvents();
    }
    glfwDestroyWindow(viewer.window);
  }
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

