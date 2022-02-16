/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <utility>
#include <vector>
#include <cstdlib>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/vec2.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/msh_primitive.h"
#include "delfem2/slice.h"
#include "delfem2/color.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/eigen/implicit_rbf_approximation.h"

// ----------------------

void Draw(
    const std::vector<double> &vtx_xy,
    const std::vector<unsigned int> &tri_vtx,
    const std::vector<double> &vtx_val,
    const std::vector<delfem2::CSegInfo> &aSeg,
    const std::vector<delfem2::CVec2d> &stroke) {
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  delfem2::opengl::DrawMeshTri2D_Edge(
      vtx_xy.data(), vtx_xy.size() / 2,
      tri_vtx.data(), tri_vtx.size() / 3);
  std::vector<std::pair<double, delfem2::CColor> > colorMap;
  delfem2::ColorMap_BlueCyanGreenYellowRed(
      colorMap,
      -1, +1);
  delfem2::opengl::DrawMeshTri2D_ScalarP1(
      vtx_xy.data(), vtx_xy.size() / 2,
      tri_vtx.data(), tri_vtx.size() / 3,
      vtx_val.data(), 1,
      colorMap);
  ::glLineWidth(5);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (auto &iseg: aSeg) {
    double pA[2], pB[2];
    iseg.Pos2D(pA, pB,
               vtx_xy.data(), tri_vtx.data());
    ::glVertex2dv(pA);
    ::glVertex2dv(pB);
  }
  ::glEnd();
  // ---
  ::glColor3d(1, 0, 0);
  ::glPointSize(5);
  ::glBegin(GL_POINTS);
  for (auto &ip: stroke) {
    ::glVertex3d(ip.x, ip.y, 0.1);
  }
  ::glEnd();
}

void Initialize(
    std::vector<double> &vtx_xy,
    std::vector<unsigned int> &tri_vtx,
    int ndiv) {
  std::vector<unsigned int> aQuad;
  delfem2::MeshQuad2D_Grid(
      vtx_xy, aQuad,
      ndiv, ndiv);
  delfem2::convert2Tri_Quad(
      tri_vtx,
      aQuad);
  delfem2::Translate_Points2(
      vtx_xy,
      -ndiv * 0.5, -ndiv * 0.5);
  delfem2::Scale_PointsX(
      vtx_xy,
      1.0 / ndiv);
}

int main() {
  std::vector<double> vtx_xy;
  std::vector<unsigned int> tri_vtx;
  Initialize(
      vtx_xy, tri_vtx,
      32);
  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  // ---------------------
  double time0 = glfwGetTime();
  glfwSetWindowTitle(viewer.window, "rbf without linear");
  while (!glfwWindowShouldClose(viewer.window)) {
    std::vector<delfem2::CVec2d> stroke;
    ImplicitRbfApproximation rbf(
        [](double r) { return std::exp(-r); },
        false);
    auto ndiv = static_cast<unsigned int>(28. * sin(glfwGetTime()) + 30.0);
    for (unsigned int i = 0; i < ndiv; ++i) {
      stroke.emplace_back(0.3 * sin(i * 0.1), 0.4 * cos(i * 0.1));
    }
    rbf.SetPolyline2(stroke, 1.0e-3);
    std::vector<double> vtx_val(vtx_xy.size() / 2, 0.0);
    for (unsigned int ip = 0; ip < vtx_xy.size() / 2; ++ip) {
      vtx_val[ip] = rbf.Evaluate2(vtx_xy[ip * 2 + 0], vtx_xy[ip * 2 + 1]);
    }
    std::vector<delfem2::CSegInfo> segments;
    delfem2::AddContour(
        segments,
        0.0,
        tri_vtx.data(), tri_vtx.size() / 3,
        vtx_val.data());
    // ------------------
    viewer.DrawBegin_oldGL();
    Draw(vtx_xy, tri_vtx, vtx_val, segments, stroke);
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwGetTime() - time0 > 6 ){ break; }
  }
  // -----------------------
  time0 = glfwGetTime();
  glfwSetWindowTitle(viewer.window, "rbf with linear");
  while (!glfwWindowShouldClose(viewer.window)) {
    std::vector<delfem2::CVec2d> stroke;
    ImplicitRbfApproximation rbf(
        [](double r) { return std::exp(-r); },
        true);
    auto ndiv = static_cast<unsigned int>(28. * sin(glfwGetTime()) + 30.0);
    for (unsigned int i = 0; i < ndiv; ++i) {
      stroke.emplace_back(0.3 * sin(i * 0.1), 0.4 * cos(i * 0.1));
    }
    rbf.SetPolyline2(stroke, 1.0e-3);
    std::vector<double> vtx_val(vtx_xy.size() / 2, 0.0);
    for (unsigned int ip = 0; ip < vtx_xy.size() / 2; ++ip) {
      vtx_val[ip] = rbf.Evaluate2(vtx_xy[ip * 2 + 0], vtx_xy[ip * 2 + 1]);
    }
    std::vector<delfem2::CSegInfo> segments;
    delfem2::AddContour(
        segments,
        0.0,
        tri_vtx.data(), tri_vtx.size() / 3,
        vtx_val.data());
    // ------------------
    viewer.DrawBegin_oldGL();
    Draw(vtx_xy, tri_vtx, vtx_val, segments, stroke);
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwGetTime() - time0 > 6 ){ break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
