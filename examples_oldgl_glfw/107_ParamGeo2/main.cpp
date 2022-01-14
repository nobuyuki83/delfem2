/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/vec2.h"
#include "delfem2/geo_curve_ndegree.h"
#include "delfem2/geo_curve_cubic.h"
// #include "delfem2/geo_curve_quadratic.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/v2.h"

namespace dfm2 = delfem2;

// -----------------------------------------------

void SetExample(
    std::vector<double> &aKnotFlat,
    unsigned int ndegree,
    unsigned int ncp) {
  const int ndiv = ncp - ndegree;
  //
  std::vector<double> aKnot;
  aKnot.assign(ndiv + 1, 0);
  for (int idiv = 0; idiv < ndiv + 1; ++idiv) {
    aKnot[idiv] = (double) idiv / ndiv;
  }
  //
  std::vector<int> aKnotMulti;
  aKnotMulti.assign(ndiv + 1, 1);
  aKnotMulti[0] = ndegree + 1;
  aKnotMulti[ndiv] = ndegree + 1;
  dfm2::FlatKnot(aKnotFlat, aKnotMulti, aKnot);
  for (unsigned int ik = 0; ik < aKnotFlat.size(); ++ik) {
    std::cout << "knot" << ik << " " << aKnotFlat[ik] << std::endl;
  }
}

void myGlutDisplay(
    const std::vector<dfm2::CVec2d> &aCtrlPoint,
    const std::vector<dfm2::CVec2d> &polyline0) {
  ::glDisable(GL_LIGHTING);

  ::glLineWidth(2);
  ::glPointSize(5);
  ::glColor3d(1, 0, 0);

  ::glColor3d(0, 0.5, 0);
  delfem2::opengl::drawPolyLine2D(aCtrlPoint);

  ::glPointSize(2);
  ::glColor3d(0, 0, 0);
  delfem2::opengl::drawPolyLine2D(polyline0);
}

int main() {
  const unsigned int ncp = 10;
  std::vector<dfm2::CVec2d> aCtrlPoint;
  const int nsmpl = 100;
  std::vector<dfm2::CVec2d> polyline0; // current test
  //
  delfem2::glfw::CViewer3 viewer(5);
  delfem2::glfw::InitGLOld();
  for (unsigned int ndegree = 2; ndegree < 5; ++ndegree) {
    viewer.OpenWindow();
    std::vector<double> aKnotFlat;
    SetExample(
        aKnotFlat,
        ndegree, ncp);
    aCtrlPoint.resize(ncp);
    while (!glfwWindowShouldClose(viewer.window)) {
      double t0 = glfwGetTime();
      for (unsigned int icp = 0; icp < ncp; ++icp) {
        aCtrlPoint[icp] = {
            icp + 0.3 * cos(icp * t0) - ncp * 0.5,
            aCtrlPoint[icp].y = sin(icp * t0)};
      }
      dfm2::SampleBSpline(
          polyline0,
          nsmpl, ndegree, aKnotFlat, aCtrlPoint);
      viewer.DrawBegin_oldGL();
      myGlutDisplay(aCtrlPoint, polyline0);
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
    }
    glfwDestroyWindow(viewer.window);
  }

  viewer.OpenWindow();
  aCtrlPoint.resize(ncp);
  while (!glfwWindowShouldClose(viewer.window)) {
    double t0 = glfwGetTime();
    for (unsigned int icp = 0; icp < ncp; ++icp) {
      aCtrlPoint[icp] = {
          icp + 0.3 * cos(icp * t0) - ncp * 0.5,
          aCtrlPoint[icp].y = sin(icp * t0)};
    }
    for (unsigned int ismpl = 0; ismpl < nsmpl + 1; ++ismpl) {
      polyline0[ismpl] = Sample_CatmullRomSplineCurve(
          static_cast<double>(ismpl) / (nsmpl) * (aCtrlPoint.size() - 1),
          aCtrlPoint);
    }
    viewer.DrawBegin_oldGL();
    myGlutDisplay(aCtrlPoint, polyline0);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);

  glfwTerminate();
  exit(EXIT_SUCCESS);
}


