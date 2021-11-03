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

#include "delfem2/femrod.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/lsmats.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/funcs.h"

namespace dfm2 = delfem2;

// -------------------------------------

void myGlutDisplay(
    const std::vector<dfm2::CVec3d> &aP,
    const std::vector<dfm2::CVec3d> &aS,
    std::vector<unsigned int> &aIP_HairRoot) {
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
    // ------------
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

class CHairShape {
 public:
  //! number of points
  unsigned int np;

  //! axis position increment for the helix
  double pitch;

  //! radius of helix
  double rad0;

  //! angle increment for the helix
  double dangle;

  //! root shape
  double p0[3];
};

void MakeProblemSetting_Spiral(
    std::vector<dfm2::CVec3d> &aP0,
    std::vector<dfm2::CVec3d> &aS0,
    std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<CHairShape> &aHairShape) {
  aIP_HairRoot.assign(1, 0);
  aP0.clear();
  aS0.clear();
  for (unsigned int ihair = 0; ihair < aHairShape.size(); ++ihair) {
    const unsigned int np = aHairShape[ihair].np;
    const double pitch = aHairShape[ihair].pitch;
    const double dangle = aHairShape[ihair].dangle;
    const double rad0 = aHairShape[ihair].rad0;
    const double *p0 = aHairShape[ihair].p0;
    for (unsigned int ip = 0; ip < np; ++ip) {
      dfm2::CVec3d p = dfm2::CVec3d(
          p0[0] + ip * pitch,
          p0[1] + rad0 * cos(dangle * ip),
          p0[2] + rad0 * sin(dangle * ip));
      aP0.push_back(p);
    }
    const unsigned int np0 = aIP_HairRoot[ihair];
    for (unsigned int is = 0; is < np - 1; ++is) {
      const dfm2::CVec3d v = (aP0[np0 + is + 1] - aP0[np0 + is + 0]).normalized();
      dfm2::CVec3d s(1.3, 1.5, 1.7);
      s = (s - (s.dot(v)) * v).normalized();
      aS0.push_back(s);
    }
    aS0.emplace_back(1, 0, 0);
    aIP_HairRoot.push_back(static_cast<unsigned int>(aP0.size()));
  }
}

int main() {
  dfm2::glfw::CViewer3 viewer(1.5);
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  // ---------------
  std::random_device rd;
  std::mt19937 reng(rd());
  std::uniform_real_distribution<double> dist01(0.0, 1.0);
  // --------------
  while (true) {
    std::vector<dfm2::CVec3d> aP0, aS0;
    std::vector<unsigned int> aIP_HairRoot;
    {
      std::vector<CHairShape> aHairShape;
      double rad0 = dist01(reng);
      double dangle = dist01(reng);
      for (int ihair = 0; ihair < 10; ++ihair) {
        CHairShape hs{30, 0.1, rad0, dangle,
                      {-1.0,
                       (dist01(reng) - 0.5) * 2.0,
                       (dist01(reng) - 0.5) * 2.0}};
        aHairShape.push_back(hs);
      }
      MakeProblemSetting_Spiral(aP0, aS0, aIP_HairRoot,
                                aHairShape); // dangle
      assert(aS0.size() == aP0.size());
    }
    for (int itr = 0; itr < 10; ++itr) {
      dfm2::ParallelTransport_RodHair(aP0, aS0, aIP_HairRoot);
    }
    std::vector<int> aBCFlag;
    dfm2::MakeBCFlag_RodHair(
        aBCFlag,
        aIP_HairRoot);
    dfm2::CMatrixSparse<double> mats;
    dfm2::MakeSparseMatrix_RodHair(
        mats,
        aIP_HairRoot, 4);
    // -----------------
    std::vector<dfm2::CVec3d> aS = aS0, aP = aP0;
    // apply random perturbation
    for (unsigned int ip = 0; ip < aP.size(); ++ip) {
      aP[ip] = aP0[ip];
      auto rnd = dfm2::CVec3d::Random(dist01, reng) * 0.1;
      if (aBCFlag[ip * 4 + 0] == 0) { aP[ip].p[0] += rnd.x; }
      if (aBCFlag[ip * 4 + 1] == 0) { aP[ip].p[1] += rnd.y; }
      if (aBCFlag[ip * 4 + 2] == 0) { aP[ip].p[2] += rnd.z; }
      if (aBCFlag[ip * 4 + 3] == 0) {
        assert(ip != aP.size() - 1);
        aS[ip] += dfm2::CVec3d::Random(dist01, reng) * 0.1;
      }
    }
    dfm2::MakeDirectorOrthogonal_RodHair(aS, aP);
    const double stiff_stretch = dist01(reng) + 1.0;
    const double stiff_bendtwist[3] = {
        dist01(reng) + 1.0,
        dist01(reng) + 1.0,
        dist01(reng) + 1.0};
    // --------------
    for (int iframe = 0; iframe < 100; ++iframe) {
      // static minimization of the rod deformation
      dfm2::MakeDirectorOrthogonal_RodHair(aS, aP);
      Solve_RodHair(
          aP, aS, mats,
          stiff_stretch, stiff_bendtwist, 0.0,
          aP0, aS0, aBCFlag, aIP_HairRoot);
      // ----------
      viewer.DrawBegin_oldGL();
      myGlutDisplay(aP, aS, aIP_HairRoot);
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
