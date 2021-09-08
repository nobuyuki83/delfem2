/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cstdlib>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/points.h"
#include "delfem2/mshuni.h"
#include "delfem2/slice.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/color.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"

// -------------------------

std::vector<double> aXY;
std::vector<unsigned int> aTri;
std::vector<double> aVal;
std::vector<delfem2::CSegInfo> aSeg;

// ---------------------------

void myGlutDisplay() {
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  delfem2::opengl::DrawMeshTri2D_Edge(
      aXY.data(), aXY.size() / 2,
      aTri.data(), aTri.size() / 3);

  std::vector<std::pair<double, delfem2::CColor> > colorMap;
  delfem2::ColorMap_BlueCyanGreenYellowRed(
      colorMap,
      -1, 1);
  delfem2::opengl::DrawMeshTri2D_ScalarP1(
      aXY.data(), aXY.size() / 2,
      aTri.data(), aTri.size() / 3,
      aVal.data(), 1,
      colorMap);
  ::glLineWidth(5);
  ::glBegin(GL_LINES);
  ::glColor3d(0, 0, 0);
  for (auto &iseg : aSeg) {
    double pA[2], pB[2];
    iseg.Pos2D(pA, pB,
               aXY.data(), aTri.data());
    ::glVertex2dv(pA);
    ::glVertex2dv(pB);
  }
  ::glEnd();

}

void Hoge() {
  { // make mesh
    const int ndiv = 16;
    std::vector<unsigned int> aQuad;
    delfem2::MeshQuad2D_Grid(
        aXY, aQuad,
        ndiv, ndiv);
    delfem2::convert2Tri_Quad(
        aTri,
        aQuad);
    delfem2::Translate_Points2(
        aXY,
        -ndiv * 0.5, -ndiv * 0.5);
    delfem2::Scale_PointsX(
        aXY,
        1.0 / ndiv);
  }

  { // make value
    const size_t np = aXY.size() / 2;
    aVal.resize(np);
    for (size_t ip = 0; ip < np; ++ip) {
      double x0 = aXY[ip * 2 + 0];
      double y0 = aXY[ip * 2 + 1];
      aVal[ip] = sqrt(x0 * x0 + y0 * y0) * 4.0 - 1.5;
    }
  }

}

int main() {
  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();

  viewer.camera.view_height = 1.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;

  Hoge();

  while (!glfwWindowShouldClose(viewer.window)) {
    {
      static int iframe = 0;
      double thres = 0.9 * sin(iframe * 0.001);
      aSeg.clear();
      delfem2::AddContour(aSeg,
                          thres,
                          aTri.data(), aTri.size() / 3,
                          aVal.data());
      iframe += 1;
    }
    // ------------------
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
