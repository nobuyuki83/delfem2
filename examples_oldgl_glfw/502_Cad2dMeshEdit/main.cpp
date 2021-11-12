/*
 * Copyright (c) 2019 Nobuyuki Umetani
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

#include "delfem2/cad2_mesh_deformation.h"
#include "delfem2/cagedef.h"
#include "delfem2/cad2.h"
#include "delfem2/cad2_mesher.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/cad2dtriv2.h"

#ifndef M_PI
#  define M_PI 3.141592653589793
#endif

namespace dfm2 = delfem2;

// ----------------------------------------------

class CCadMesh2DVeiwer : public delfem2::glfw::CViewer2 {
 public:
  explicit CCadMesh2DVeiwer(delfem2::CCad2D &cad)
      : CViewer2(), cad(cad) {
    {
      delfem2::CBoundingBox2<double> bb = cad.BB();
      this->view_height = static_cast<float>(bb.LengthDiagonal());
    }
    cad.is_draw_face = false;
    {
      delfem2::CMeshDynTri2D dmsh;
      mesher.edge_length = 0.1;
      mesher.Meshing(dmsh, cad);
      dmsh.Export_StlVectors(vtx_xy, tri_vtx);
    }
    vtx_w.resize(vtx_xy.size() / 2, 0.0);
  }
  void mouse_press(
      const float src[2]) override {
    picked_pos = {src[0], src[1]};
    vtx_xy_when_picked = vtx_xy;
    cad.Pick(src[0], src[1], view_height);
    if (cad.ivtx_picked != UINT_MAX) { // vertex is picked
      picked_pos = {
          static_cast<float>(cad.aVtx[cad.ivtx_picked].pos.x),
          static_cast<float>(cad.aVtx[cad.ivtx_picked].pos.y)};
    }
    SetCadMeshDeformationWeight(
        vtx_w,
        cad, mesher, vtx_xy);
  }
  void mouse_drag(
      const float src0[2],
      const float src1[2]) override {
    cad.DragPicked(src1[0], src1[1], src0[0], src0[1]);
    for (unsigned int ip = 0; ip < vtx_xy.size() / 2; ++ip) {
      vtx_xy[ip * 2 + 0] = vtx_xy_when_picked[ip * 2 + 0] + vtx_w[ip] * (src1[0] - picked_pos[0]);
      vtx_xy[ip * 2 + 1] = vtx_xy_when_picked[ip * 2 + 1] + vtx_w[ip] * (src1[1] - picked_pos[1]);
    }
  }
  void Draw() {
    this->DrawBegin_oldGL();
    delfem2::opengl::Draw_CCad2D(cad);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0.8, 0.8, 0.8);
    delfem2::opengl::DrawMeshTri2D_Face(tri_vtx, vtx_xy);
    ::glLineWidth(1);
    delfem2::opengl::DrawMeshTri2D_Edge(tri_vtx, vtx_xy);
    this->SwapBuffers();
  }
 public:
  delfem2::CCad2D &cad;
  delfem2::CMesher_Cad2D mesher;
  std::vector<double> vtx_xy;
  std::vector<unsigned int> tri_vtx;
  std::vector<double> vtx_xy_when_picked;
  std::vector<double> vtx_w;
  std::array<float, 2> picked_pos{0.f, 0.f};
};

int main() {
  delfem2::glfw::InitGLOld();
  // ----------------------------------
  { // case 0
    delfem2::CCad2D cad;
    cad.AddPolygon({-1, -1, +1, -1, +1, +1, -1, +1});
    CCadMesh2DVeiwer viewer(cad);
    viewer.OpenWindow();
    while (!glfwWindowShouldClose(viewer.window)) {
      viewer.Draw();
      glfwPollEvents();
    }
    glfwDestroyWindow(viewer.window);
  }
  { // case 1
    delfem2::CCad2D cad;
    cad.AddPolygon({
      -1, -1,
      +1, -1,
      +1, +1,
      -1, +1});
    cad.aEdge[0].SetQuadraticBezierCurve({0,-0.5});
    cad.aEdge[2].SetCubicBezierCurve({0.5, 1.5}, {-0.5, 1.5});
    CCadMesh2DVeiwer viewer(cad);
    viewer.OpenWindow();
    while (!glfwWindowShouldClose(viewer.window)) {
      viewer.Draw();
      glfwPollEvents();
    }
    glfwDestroyWindow(viewer.window);
  }
  { // case 2
    delfem2::CCad2D cad;
    cad.AddPolygon({
      -0.84, +1.20,
      -1.00, +0.25,
      -1.00, -1.40,
      +1.00, -1.40,
      +1.00, +0.25,
      +0.84, +1.20,
      +0.30, +1.40,
      -0.30, +1.40});
    cad.aEdge[0].SetCubicBezierCurve({-0.74,+1.0}, {-0.7, +0.45});
    cad.aEdge[4].SetCubicBezierCurve({+0.7,+0.45},{0.74,+1.0});
    cad.aEdge[6].SetCubicBezierCurve({+0.1,+1.2},{-0.1,+1.2});
    CCadMesh2DVeiwer viewer(cad);
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


