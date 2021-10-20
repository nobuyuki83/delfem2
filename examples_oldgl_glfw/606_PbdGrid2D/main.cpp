/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/pbd_geo3.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/jagarray.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// -------------------------------------------------

void stepTime(
    std::vector<double> &aXY1,
    std::vector<double> &aUV1,
    std::vector<double> &aTmp,
    double dt,
    int nitr,
    const std::vector<unsigned int> &clstr_ind,
    const std::vector<unsigned int> &clstr,
    const std::vector<int> &aBC,
    [[maybe_unused]] const std::vector<unsigned int> &aQuad,
    const std::vector<double> &aXY0) {
  const size_t ndof = aXY0.size();
  for (unsigned int idof = 0; idof < ndof; idof++) {
    aTmp[idof] = aXY1[idof] + aUV1[idof] * dt;
  }
  for (size_t ip = 0; ip < aXY0.size() / 2; ++ip) {
    if (aBC[ip] == 0) { continue; }
    aTmp[ip * 2 + 0] = aXY1[ip * 2 + 0];
    aTmp[ip * 2 + 1] = aXY1[ip * 2 + 1];
  }
  // deform
  for (int itr = 0; itr < nitr; itr++) {
    dfm2::PBD_ConstProj_Rigid2D(
        aTmp.data(),
        0.5,
        clstr_ind.data(), clstr_ind.size(),
        clstr.data(), clstr.size(),
        aXY0.data(), aXY0.size());
  }
  for (size_t ip = 0; ip < aXY0.size() / 2; ++ip) {
    if (aBC[ip] == 0) { continue; }
    aTmp[ip * 2 + 0] = aXY1[ip * 2 + 0];
    aTmp[ip * 2 + 1] = aXY1[ip * 2 + 1];
  }
  for (unsigned int idof = 0; idof < ndof; ++idof) {
    aUV1[idof] = (aTmp[idof] - aXY1[idof]) * (1.0 / dt);
  }
  for (unsigned int idof = 0; idof < ndof; idof++) {
    aXY1[idof] = aTmp[idof];
  }
}

int main() {
  // --------------------------
  std::vector<double> aXY0;
  std::vector<double> aXY1;
  std::vector<double> aUV1;
  std::vector<double> aXYt;
  std::vector<unsigned int> aQuad;
  std::vector<unsigned int> clstr_ind, clstr;
  std::vector<int> aBC;
  const int nX = 8;
  const int nY = 8;
  delfem2::MeshQuad2D_Grid(aXY0, aQuad, nX, nY);
  aBC.assign(aXY0.size() / 2, 0);
  for (int ix = 0; ix < nX + 1; ++ix) { aBC[ix] = 1; }
  aXY1 = aXY0;
  aXYt = aXY0;
  aUV1.resize(aXY0.size());
  {
    std::vector<unsigned int> psup_ind, psup;
    dfm2::JArray_PSuP_MeshElem(
        psup_ind, psup,
        aQuad.data(), aQuad.size() / 4, 4,
        aXY0.size() / 2);
    dfm2::JArray_AddDiagonal(
        clstr_ind, clstr,
        psup_ind.data(), psup_ind.size(),
		psup.data(), psup.size());
  }

  dfm2::glfw::CViewer2 viewer;
  {
    viewer.view_height = 10.0;
    viewer.trans[0] = -5.0;
    viewer.trans[1] = -5.0;
  }
  delfem2::glfw::InitGLOld();
  viewer.CreateWindow();

  const double dt = 1.0 / 60.0; // frame-rate is fixed to 60FPS.
  double time_last_update = 0.0;
  while (!glfwWindowShouldClose(viewer.window)) {
    // control of the frame rate
    const double time_now = glfwGetTime();
    if (time_now - time_last_update < dt) {
      glfwPollEvents();
      continue;
    }
    time_last_update = time_now;
    {
      for (int ix = 0; ix < nX + 1; ++ix) {
        aXY1[ix * 2 + 0] = ix + 2 * sin(time_now * 10);
        aXY1[ix * 2 + 1] = 0;
      }
      stepTime(aXY1, aUV1, aXYt,
               dt, 1,
               clstr_ind, clstr,
               aBC,
               aQuad, aXY0);
    }
    //
    viewer.DrawBegin_oldGL();
    delfem2::opengl::DrawMeshQuad2D_Edge(
        aXY1.data(), aXY1.size() / 2,
        aQuad.data(), aQuad.size() / 4);
    viewer.SwapBuffers();
    glfwPollEvents();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

