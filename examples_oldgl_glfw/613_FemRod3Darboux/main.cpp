/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/fem_rod3_darboux.h"
#include "delfem2/fem_distance3.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/lsmats.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/view_vectorx.h"
#include "delfem2/mshuni.h"
#include "delfem2/jagarray.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/funcs.h"

namespace dfm2 = delfem2;

// -------------------------------------

void myGlutDisplay(
    const std::vector<dfm2::CVec3d> &aP,
    const std::vector<dfm2::CVec3d> &aS,
    const std::vector<unsigned int> &aElemSeg) {
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1, 0, 0);
  ::glPointSize(10);
  ::glBegin(GL_POINTS);
  for (const auto &p : aP) {
    ::glVertex3d(p.x, p.y, p.z);
  }
  ::glEnd();
  // ------------
  ::glColor3d(0, 0, 0);
  ::glLineWidth(3);
  ::glBegin(GL_LINES);
  for (unsigned int iseg = 0; iseg < aElemSeg.size() / 2; ++iseg) {
    unsigned int i0 = aElemSeg[iseg * 2 + 0];
    assert(i0 < aP.size());
    unsigned int i1 = aElemSeg[iseg * 2 + 1];
    assert(i1 < aP.size());
    ::glVertex3d(aP[i0].x, aP[i0].y, aP[i0].z);
    ::glVertex3d(aP[i1].x, aP[i1].y, aP[i1].z);
  }
  ::glEnd();
  // --------------
  ::glBegin(GL_LINES);
  for (unsigned int iseg = 0; iseg < aElemSeg.size() / 2; ++iseg) {
    unsigned int i0 = aElemSeg[iseg * 2 + 0];
    assert(i0 < aP.size());
    unsigned int i1 = aElemSeg[iseg * 2 + 1];
    assert(i1 < aP.size());
    dfm2::CVec3d p01 = 0.5 * (aP[i0] + aP[i1]);
    double l01 = (aP[i0] - aP[i1]).norm();
    dfm2::opengl::myGlVertex(p01);
    dfm2::opengl::myGlVertex(p01 + (l01 * 0.5) * aS[iseg]);
  }
  ::glEnd();
}

void MakeProblemSetting_Spiral(
    std::vector<dfm2::CVec3d> &aP0,
    std::vector<dfm2::CVec3d> &aS0,
    std::vector<unsigned int> &aElemSeg,
    std::vector<unsigned int> &aElemRod,
    unsigned int np,
    double pitch,
    double rad0,
    double dangle) {
  aP0.resize(np);
  for (unsigned int ip = 0; ip < np; ++ip) {
    aP0[ip] = dfm2::CVec3d(
        -1.0 + ip * pitch,
        rad0 * cos(dangle * ip),
        rad0 * sin(dangle * ip));
  };
  // -------------------------
  // below: par segment data
  const unsigned int ns = np - 1;
  aElemSeg.resize(ns * 2);
  for (unsigned int is = 0; is < ns; ++is) {
    aElemSeg[is * 2 + 0] = is + 0;
    aElemSeg[is * 2 + 1] = is + 1;
  }
  { // initial director vector
    aS0.resize(ns, dfm2::CVec3d(1, 0, 0));
    for (unsigned int is = 0; is < ns; ++is) {
      unsigned int ip0 = aElemSeg[is * 2 + 0];
      unsigned int ip1 = aElemSeg[is * 2 + 1];
      const dfm2::CVec3d v = (aP0[ip1] - aP0[ip0]).normalized();
      aS0[is] = (aS0[is] - (aS0[is].dot(v)) * v).normalized();
    }
  }
  // --------------------------
  // below: par rod element data
  const unsigned int nr = ns - 1;
  aElemRod.resize(nr * 5);
  for (unsigned int ir = 0; ir < nr; ++ir) {
    aElemRod[ir * 5 + 0] = ir + 0;
    aElemRod[ir * 5 + 1] = ir + 1;
    aElemRod[ir * 5 + 2] = ir + 2;
    aElemRod[ir * 5 + 3] = np + ir + 0;
    aElemRod[ir * 5 + 4] = np + ir + 1;
  };
  // smoothing
  for (int itr = 0; itr < 10; ++itr) {
    for (unsigned int ir = 0; ir < nr; ++ir) {
      const unsigned int ip0 = aElemRod[ir * 5 + 0];
      const unsigned int ip1 = aElemRod[ir * 5 + 1];
      const unsigned int ip2 = aElemRod[ir * 5 + 2];
      const unsigned int is0 = aElemRod[ir * 5 + 3] - np;
      assert(is0 < ns);
      const unsigned int is1 = aElemRod[ir * 5 + 4] - np;
      assert(is1 < ns);
      const dfm2::CMat3d CMat3 = dfm2::Mat3_MinimumRotation(aP0[ip1] - aP0[ip0], aP0[ip2] - aP0[ip1]);
      dfm2::CVec3d s1 = CMat3 * aS0[is0] + aS0[is1];
      const dfm2::CVec3d v = (aP0[ip2] - aP0[ip1]).normalized();
      aS0[is1] = (s1 - (s1.dot(v)) * v).normalized();
    }
  }
}

DFM2_INLINE void Solve_DispRotSeparate(
    std::vector<delfem2::CVec3d> &aP,
    std::vector<delfem2::CVec3d> &aS,
    delfem2::CMatrixSparse<double> &mats,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    const std::vector<delfem2::CVec3d> &aP0,
    const std::vector<delfem2::CVec3d> &aS0,
    const std::vector<unsigned int> &aElemSeg,
    const std::vector<unsigned int> &aElemRod,
    const std::vector<int> &aBCFlag) {
  using namespace delfem2;
  assert(mats.nrowdim_ == 3);
  assert(mats.ncoldim_ == 3);
  const size_t nNode = aBCFlag.size() / 3;
  assert(aP.size() + aS.size() == nNode);
  mats.setZero();
  std::vector<double> vec_r;
  vec_r.assign(nNode * 3, 0.0);
  std::vector<unsigned int> tmp_buffer;
  double W = 0;
  for (unsigned int iseg = 0; iseg < aElemSeg.size() / 2; ++iseg) {
    const unsigned int i0 = aElemSeg[iseg * 2 + 0];
    const unsigned int i1 = aElemSeg[iseg * 2 + 1];
    const unsigned int *aINoel = aElemSeg.data() + iseg * 2;
    const double L0 = (aP0[i0] - aP0[i1]).norm();
    const CVec3d aPE[2] = {aP[i0], aP[i1]};
    // --------------
    CVec3d dW_dP[2];
    CMat3d ddW_ddP[2][2];
    W += WdWddW_SquareLengthLineseg3D(dW_dP, ddW_ddP,
                                      stiff_stretch, aPE, L0);
    {
      double eM[2][2][3][3];
      for (int in = 0; in < 2; ++in) {
        for (int jn = 0; jn < 2; ++jn) {
          ddW_ddP[in][jn].CopyTo(&eM[in][jn][0][0]);
        }
      }
//      mats.Mearge(2, aINoel, 2, aINoel, 9, eM, tmp_buffer);
      Merge<2, 2, 3, 3, double>(mats, aINoel, aINoel, eM, tmp_buffer);
    }
    {
      for (int inoel = 0; inoel < 2; inoel++) {
        const unsigned int ip = aINoel[inoel];
        vec_r[ip * 3 + 0] -= dW_dP[inoel].x;
        vec_r[ip * 3 + 1] -= dW_dP[inoel].y;
        vec_r[ip * 3 + 2] -= dW_dP[inoel].z;
      }
    }
  }
  for (unsigned int irod = 0; irod < aElemRod.size() / 5; ++irod) {
    const unsigned int *aINoel = aElemRod.data() + irod * 5;
    const size_t nP = aP.size();
    const CVec3d aPE[3] = {aP[aINoel[0]], aP[aINoel[1]], aP[aINoel[2]]};
    const CVec3d aSE[2] = {aS[aINoel[3] - nP], aS[aINoel[4] - nP]};
    CVec3d Darboux0;
    {
      const CVec3d aPE0[3] = {aP0[aINoel[0]], aP0[aINoel[1]], aP0[aINoel[2]]};
      const CVec3d aSE0[2] = {aS0[aINoel[3] - nP], aS0[aINoel[4] - nP]};
      Darboux0 = Darboux_Rod(aPE0, aSE0);
    }
    // ------
    CVec3d dW_dP[3];
    double dW_dt[2];
    CMat3d ddW_ddP[3][3];
    CVec3d ddW_dtdP[2][3];
    double ddW_ddt[2][2];
    W += WdWddW_Rod3Approx(
        dW_dP, dW_dt, ddW_ddP, ddW_dtdP, ddW_ddt,
        stiff_bendtwist,
        aPE, aSE, Darboux0);
    {
      double eM[5][5][3][3];
      for (int i = 0; i < 5 * 5 * 3 * 3; ++i) { (&eM[0][0][0][0])[i] = 0.0; }
      for (int in = 0; in < 3; ++in) {
        for (int jn = 0; jn < 3; ++jn) {
          ddW_ddP[in][jn].CopyTo(&eM[in][jn][0][0]);
        }
      }
      for (int in = 0; in < 3; ++in) {
        for (int jn = 0; jn < 2; ++jn) {
          eM[3 + jn][in][0][0] = eM[in][jn + 3][0][0] = ddW_dtdP[jn][in].x;
          eM[3 + jn][in][0][1] = eM[in][jn + 3][1][0] = ddW_dtdP[jn][in].y;
          eM[3 + jn][in][0][2] = eM[in][jn + 3][2][0] = ddW_dtdP[jn][in].z;
        }
      }
      for (int in = 0; in < 2; ++in) {
        for (int jn = 0; jn < 2; ++jn) {
          eM[in + 3][jn + 3][0][0] = ddW_ddt[in][jn];
        }
      }
      Merge<5, 5, 3, 3, double>(mats, aINoel, aINoel, eM, tmp_buffer);
//      mats.Mearge(5, aINoel, 5, aINoel, 9, &eM[0][0][0][0], tmp_buffer);
    }
    {
      for (int inoel = 0; inoel < 3; inoel++) {
        const unsigned int ip = aINoel[inoel];
        vec_r[ip * 3 + 0] -= dW_dP[inoel].x;
        vec_r[ip * 3 + 1] -= dW_dP[inoel].y;
        vec_r[ip * 3 + 2] -= dW_dP[inoel].z;
      }
      for (int inoel = 0; inoel < 2; inoel++) {
        const unsigned int in0 = aINoel[3 + inoel];
        vec_r[in0 * 3 + 0] -= dW_dt[inoel];
      }
    }
  }
  //  std::cout << CheckSymmetry(mats) << std::endl;
  //  mats.AddDia(0.00001);
  std::cout << "energy:" << W << std::endl;
  //    std::cout << "sym: " << CheckSymmetry(mats) << std::endl;
  mats.SetFixedBC(aBCFlag.data());
  setRHS_Zero(vec_r, aBCFlag, 0);
  std::vector<double> vec_x;
  vec_x.assign(nNode * 3, 0.0);
  {
    const std::size_t n = vec_r.size();
    std::vector<double> tmp0(n), tmp1(n);
    auto vr = delfem2::ViewAsVectorXd(vec_r);
    auto vu = delfem2::ViewAsVectorXd(vec_x);
    auto vs = delfem2::ViewAsVectorXd(tmp0);
    auto vt = delfem2::ViewAsVectorXd(tmp1);
    auto aConvHist = delfem2::Solve_CG(
        vr, vu, vs, vt,
        1.0e-4, 300, mats);
    if (!aConvHist.empty()) {
      std::cout << "            conv: " << aConvHist.size();
      std::cout << " " << aConvHist[0];
      std::cout << " " << aConvHist[aConvHist.size() - 1] << std::endl;
    }
  }
  /*
   {
   auto aConvHist = Solve_BiCGStab(vec_r,vec_x,
   1.0e-4, 300, mats);
   if( aConvHist.size() > 0 ){
   std::cout << "            conv: " << aConvHist.size() << " " << aConvHist[0] << " " << aConvHist[aConvHist.size()-1] << std::endl;
   }
   }
   */
  assert(aS.size() == aElemSeg.size() / 2);
  for (unsigned int is = 0; is < aS.size(); ++is) {
    const unsigned int i0 = aElemSeg[is * 2 + 0];
    const unsigned int i1 = aElemSeg[is * 2 + 1];
    CVec3d V01 = aP[i1] - aP[i0];
    CVec3d du(vec_x[i1 * 3 + 0] - vec_x[i0 * 3 + 0],
              vec_x[i1 * 3 + 1] - vec_x[i0 * 3 + 1],
              vec_x[i1 * 3 + 2] - vec_x[i0 * 3 + 2]);
    const size_t np = aP.size();
    const double dtheta = vec_x[np * 3 + is * 3];
    CVec3d frm[3];
    RodFrameTrans(frm,
                  aS[is], V01, du, dtheta);
    aS[is] = frm[0];
  }
  for (unsigned int ip = 0; ip < aP.size(); ++ip) {
    aP[ip].p[0] += vec_x[ip * 3 + 0];
    aP[ip].p[1] += vec_x[ip * 3 + 1];
    aP[ip].p[2] += vec_x[ip * 3 + 2];
  }
  for (unsigned int iseg = 0; iseg < aElemSeg.size() / 2; ++iseg) {
    const unsigned int i0 = aElemSeg[iseg * 2 + 0];
    const unsigned int i1 = aElemSeg[iseg * 2 + 1];
    const CVec3d &p0 = aP[i0];
    const CVec3d &p1 = aP[i1];
    const CVec3d e01 = (p1 - p0).normalized();
    aS[iseg] -= (aS[iseg].dot(e01)) * e01;
    aS[iseg].normalize();
  }
}

int main() {
  dfm2::glfw::CViewer3 viewer(1.5);
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  //
  std::random_device rd;
  std::mt19937 reng(rd());
  std::uniform_real_distribution<double> dist01(0.0, 1.0);
  std::uniform_real_distribution<double> dist03(0.0, 3.0);
  // ------
  while (true) {
    std::vector<dfm2::CVec3d> aP0, aS0;
    std::vector<unsigned int> aElemSeg, aElemRod;
    MakeProblemSetting_Spiral(aP0, aS0,
                              aElemSeg, aElemRod,
                              30,
                              0.1,  // np
                              dist01(reng),  // rad0
                              dist01(reng));  // dangle
    std::vector<int> aBCFlag;
    {
      const auto np = static_cast<unsigned int>(aP0.size());
      const auto ns = static_cast<unsigned int>(aS0.size());
      const unsigned int nNode = np + ns;
      aBCFlag.assign(nNode * 3, 0);
      {
        aBCFlag[0 * 3 + 0] = 1;
        aBCFlag[0 * 3 + 1] = 1;
        aBCFlag[0 * 3 + 2] = 1;
        aBCFlag[1 * 3 + 0] = 1;
        aBCFlag[1 * 3 + 1] = 1;
        aBCFlag[1 * 3 + 2] = 1;
        aBCFlag[(np + 0) * 3 + 0] = 1;  //
        for (unsigned int is = 0; is < ns; ++is) {
          aBCFlag[(np + is) * 3 + 1] = 1;  // fix the unused dof
          aBCFlag[(np + is) * 3 + 2] = 1;  // fix the unused dof
        }
      }
    }
    dfm2::CMatrixSparse<double> mats;
    {
      size_t nNode = aElemSeg.size() / 2 + aP0.size();
      std::vector<unsigned int> psup_ind, psup;
      dfm2::JArray_PSuP_MeshElem(
          psup_ind, psup,
          aElemRod.data(), aElemRod.size() / 5, 5, nNode);
      dfm2::JArray_Sort(psup_ind, psup);
      mats.Initialize(nNode, 3, true);
      mats.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
    }
    // -----------------
    std::vector<dfm2::CVec3d> aS = aS0, aP = aP0;
    assert(aS.size() == aElemSeg.size() / 2);
    // apply random deviation
    for (unsigned int ip = 0; ip < aP.size(); ++ip) {
      aP[ip] = aP0[ip];
      if (aBCFlag[ip * 3 + 0] == 0) { aP[ip].p[0] += dist03(reng); }
      if (aBCFlag[ip * 3 + 1] == 0) { aP[ip].p[1] += dist03(reng); }
      if (aBCFlag[ip * 3 + 2] == 0) { aP[ip].p[2] += dist03(reng); }
    }
    const auto ns = static_cast<unsigned int>(aS.size());
    for (unsigned int is = 0; is < ns; ++is) {
      aS[is] = aS0[is];
      const auto np = static_cast<unsigned int>(aP.size());
      if (aBCFlag[(np + is) * 3 + 0] == 0) { aS[is].p[0] += dist03(reng); }
      if (aBCFlag[(np + is) * 3 + 1] == 0) { aS[is].p[1] += dist03(reng); }
      if (aBCFlag[(np + is) * 3 + 2] == 0) { aS[is].p[2] += dist03(reng); }
    }
    for (unsigned int iseg = 0; iseg < aElemSeg.size() / 2; ++iseg) {
      const unsigned int i0 = aElemSeg[iseg * 2 + 0];
      const unsigned int i1 = aElemSeg[iseg * 2 + 1];
      const dfm2::CVec3d &p0 = aP[i0];
      const dfm2::CVec3d &p1 = aP[i1];
      const dfm2::CVec3d e01 = (p1 - p0).normalized();
      assert(iseg < aS.size());
      aS[iseg] -= (aS[iseg].dot(e01)) * e01;
      aS[iseg].normalize();
    }
    const double stiff_stretch = dist01(reng) + 1.;
    const double stiff_bendtwist[3] = {
        dist01(reng) + 1.,
        dist01(reng) + 1.,
        dist01(reng) + 1.};
    for (int iframe = 0; iframe < 70; ++iframe) {
      Solve_DispRotSeparate(
          aP, aS, mats,
          stiff_stretch, stiff_bendtwist,
          aP0, aS0, aElemSeg, aElemRod, aBCFlag);
      viewer.DrawBegin_oldGL();
      myGlutDisplay(aP, aS, aElemSeg);
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
