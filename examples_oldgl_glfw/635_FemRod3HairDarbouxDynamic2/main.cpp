/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/hair_darboux.h"
#include "delfem2/lsmats.h"
#include "delfem2/hair_sparse.h"
#include "delfem2/mshuni.h"
#include "delfem2/ls_pentadiagonal.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/funcs.h"

namespace dfm2 = delfem2;

// -------------------------------------t

void myGlutDisplay(
    const std::vector<dfm2::CVec3d> &aP,
    const std::vector<dfm2::CVec3d> &aS,
    const std::vector<unsigned int> &aIP_HairRoot) {
  assert(!aIP_HairRoot.empty());
  const unsigned int nhair = static_cast<unsigned int>(aIP_HairRoot.size()) - 1;
  for (unsigned int ihair = 0; ihair < nhair; ++ihair) {
    const unsigned int ips = aIP_HairRoot[ihair];
    const unsigned int ipe = aIP_HairRoot[ihair + 1];
    assert(aP.size() == aS.size());
    ::glDisable(GL_LIGHTING);
    ::glColor3d(1, 0, 0);
    ::glPointSize(3);
    ::glBegin(GL_POINTS);
    for (unsigned int ip = ips; ip < ipe; ++ip) {
      ::glVertex3d(aP[ip].x, aP[ip].y, aP[ip].z);
    }
    ::glEnd();
    // ------------
    ::glColor3d(0, 0, 0);
    ::glLineWidth(3);
    ::glBegin(GL_LINES);
    unsigned int ns = ipe - ips - 1;
    for (unsigned int is = 0; is < ns; ++is) {
      const unsigned int ip0 = ips + is + 0;
      assert(ip0 < aP.size());
      const unsigned int ip1 = ips + is + 1;
      assert(ip1 < aP.size());
      ::glVertex3d(aP[ip0].x, aP[ip0].y, aP[ip0].z);
      ::glVertex3d(aP[ip1].x, aP[ip1].y, aP[ip1].z);
    }
    ::glEnd();
    // --------------
    ::glBegin(GL_LINES);
    for (unsigned int is = 0; is < ns; ++is) {
      const unsigned int ip0 = ips + is + 0;
      assert(ip0 < aP.size());
      const unsigned int ip1 = ips + is + 1;
      assert(ip1 < aP.size());
      dfm2::CVec3d p01 = 0.5 * (aP[ip0] + aP[ip1]);
      double l01 = (aP[ip0] - aP[ip1]).norm();
      dfm2::opengl::myGlVertex(p01);
      dfm2::opengl::myGlVertex(p01 + (l01 * 0.5) * aS[is]);
    }
    ::glEnd();
  }
}

int main() {
  std::mt19937 reng(std::random_device{}());
  std::uniform_real_distribution<double> dist01(0.0, 1.0);
  dfm2::glfw::CViewer3 viewer(1.5);
  // -----------------------
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  // -------
  while (true) {
    std::vector<dfm2::CVec3d> aP0; // initial position
    std::vector<dfm2::CVec3d> aS0; // initial director vector
    { // make the un-deformed shape of hair
      const delfem2::CHairShape hs{
        10, 0.1, dist01(reng), dist01(reng),
        {-1., 0., 0.} };
      std::vector<unsigned int> aIP;
      MakeProblemSetting_Spiral(
          aP0, aS0, aIP,
          {hs});
      assert(aS0.size() == aP0.size());
    }
    for (int itr = 0; itr < 10; ++itr) { // relax director vectors
      ParallelTransport_RodHair(
          aP0, aS0,
          {0,static_cast<unsigned int>(aP0.size())});
    }
    std::vector<int> aBCFlag; // boundary condition
    dfm2::MakeBCFlag_RodHair( // set fixed boundary condition
        aBCFlag,
        {0,static_cast<unsigned int>(aP0.size())});
    assert(aBCFlag.size() == aP0.size() * 4);

    dfm2::CMatrixSparse<double> mats; // sparse matrix
    {
      std::vector<unsigned int> psup_ind, psup;
      delfem2::JArray_PSuP_Hair(
          psup_ind,psup,
          {0,static_cast<unsigned int>(aP0.size())});
      mats.Initialize(aP0.size(),4,true);
      mats.SetPattern(
          psup_ind.data(), psup_ind.size(),
          psup.data(), psup.size());
    }
    dfm2::BlockPentaDiagonalMatrix<4> pdiamat;
    std::cout << aP0.size() << std::endl;
    pdiamat.Initialize(aP0.size());
    // -----------------
    std::vector<dfm2::CVec3d> aP = aP0, aS = aS0;
    std::vector<dfm2::CVec3d> aPV(aP0.size(), dfm2::CVec3d(0, 0, 0)); // velocity
    std::vector<dfm2::CVec3d> aPt = aP; // temporally positions
    double dt = 0.01;
    double mass = 1.0e-2;
    dfm2::CVec3d gravity(0, -10, 0);
    const double stiff_stretch = 10000 * (dist01(reng) + 1.);
    const double stiff_bendtwist[3] = {
        1000 * (dist01(reng) + 1.),
        1000 * (dist01(reng) + 1.),
        1000 * (dist01(reng) + 1.)};
    for (int iframe = 0; iframe < 100; ++iframe) {
      for (unsigned int ip = 0; ip < aP.size(); ++ip) {
        if (aBCFlag[ip * 4 + 0] == 0) {
          aPt[ip] = aP[ip] + dt * aPV[ip] + (dt * dt / mass) * gravity;
        }
      }
      dfm2::MakeDirectorOrthogonal_RodHair(aS, aPt);
      Solve_RodHair(
          aPt, aS, mats,
          stiff_stretch, stiff_bendtwist, mass / (dt * dt),
          aP0, aS0, aBCFlag, {0,static_cast<unsigned int>(aP.size())});
      for (unsigned int ip = 0; ip < aP.size(); ++ip) {
        if (aBCFlag[ip * 4 + 0] != 0) { continue; }
        aPV[ip] = (aPt[ip] - aP[ip]) / dt;
        aP[ip] = aPt[ip];
      }
      // -------------
      viewer.DrawBegin_oldGL();
      myGlutDisplay(
          aP, aS,
          {0,static_cast<unsigned int>(aP.size())});
      viewer.SwapBuffers();
      glfwPollEvents();
      if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
    }
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
