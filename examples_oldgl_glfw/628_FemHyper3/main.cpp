/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <random>
#include <vector>
#include <cstdlib>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/mshprimitive.h"
#include "delfem2/femsolidhyper.h"
#include "delfem2/lsilu_mats.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/lsmats.h"
#include "delfem2/lsvecx.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/jagarray.h"
#include "delfem2/mshuni.h"
#include "delfem2/femutil.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// -------------------------------------------------------

void InitializeMatrix(
    dfm2::CMatrixSparse<double> &smat,
    delfem2::CPreconditionerILU<double> &ilu,
    const std::vector<unsigned int> &aHex,
    const std::vector<double> &aXYZ) {
  const size_t np = aXYZ.size() / 3;
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind, psup,
      aHex.data(), aHex.size() / 8, 8,
      aXYZ.size() / 3);
  dfm2::JArray_Sort(psup_ind, psup);
  smat.Initialize(np, 3, true);
  smat.SetPattern(
      psup_ind.data(), psup_ind.size(),
      psup.data(), psup.size());
  ilu.SetPattern0(smat);
}

void Simulation(
    std::vector<double> &aDisp,
    std::vector<double> &aVelo,
    dfm2::CMatrixSparse<double> &smat,
    delfem2::CPreconditionerILU<double> &silu,
    //
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aHex,
    const std::vector<int> &aBCFlag,
    //
    double dt,
    double mass,
    double c1,
    double c2,
    double stiff_comp,
    const double gravity[3]) {
  const size_t np = aXYZ0.size() / 3;
  const size_t nDoF = np * 3;
  for (unsigned int ip = 0; ip < np; ++ip) {
    aDisp[ip * 3 + 0] += dt * aVelo[ip * 3 + 0];
    aDisp[ip * 3 + 1] += dt * aVelo[ip * 3 + 1];
    aDisp[ip * 3 + 2] += dt * aVelo[ip * 3 + 2];
  }
  std::vector<unsigned int> tmp_buffer;
  smat.setZero();
  std::vector<double> vecb(nDoF, 0.0);
  for (unsigned int ih = 0; ih < aHex.size() / 8; ++ih) {
    double aP0[8][3];
    delfem2::FetchData<8, 3>(aP0, aHex.data() + ih * 8, aXYZ0.data());
    double aU[8][3];
    delfem2::FetchData<8, 3>(aU, aHex.data() + ih * 8, aDisp.data());
    //
    double W = 0.0;
    double dW[8][3];
    for (int i = 0; i < 8 * 3; ++i) { (&dW[0][0])[i] = 0.0; }
    double ddW[8][8][3][3];
    for (int i = 0; i < 8 * 8 * 3 * 3; ++i) { (&ddW[0][0][0][0])[i] = 0.0; }
    double vol = 0.0;
    delfem2::AddWdWddW_Solid3HyperMooneyrivlin2Reduced_Hex(
        W, dW, ddW, vol,
        c1, c2, aP0, aU, 1);
    {
      double vol1 = 0.0;
      delfem2::AddWdWddW_Solid3Compression_Hex(
          W, dW, ddW, vol1,
          stiff_comp, aP0, aU, 0);
    }
    const unsigned int *aIP = aHex.data() + ih * 8;
    for (unsigned int ino = 0; ino < 8; ++ino) {
      unsigned int ip0 = aIP[ino];
      vecb[ip0 * 3 + 0] -= dW[ino][0];
      vecb[ip0 * 3 + 1] -= dW[ino][1];
      vecb[ip0 * 3 + 2] -= dW[ino][2];
    }
    delfem2::Merge<8, 8, 3, 3, double>(smat, aIP, aIP, ddW, tmp_buffer);
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    smat.val_dia_[ip * 9 + 0] += mass / (dt * dt);
    smat.val_dia_[ip * 9 + 4] += mass / (dt * dt);
    smat.val_dia_[ip * 9 + 8] += mass / (dt * dt);
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    vecb[ip * 3 + 0] += mass * gravity[0];
    vecb[ip * 3 + 1] += mass * gravity[1];
    vecb[ip * 3 + 2] += mass * gravity[2];
  }
  smat.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vecb, aBCFlag, 0);
  // --------------------------------
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  silu.CopyValue(smat);
  silu.Decompose();
  std::vector<double> vecx(vecb.size());
  {
    const std::size_t n = vecb.size();
    std::vector<double> tmp0(n), tmp1(n);
    std::vector<double> aHist = Solve_PCG(
        dfm2::CVecXd(vecb),
        dfm2::CVecXd(vecx),
        dfm2::CVecXd(tmp0),
        dfm2::CVecXd(tmp1),
        conv_ratio, iteration,
        smat, silu);
//    std::cout << "nconv:" << aHist.size() << std::endl;
  }
  // -----------------------------
  dfm2::XPlusAY(
      aDisp,
      nDoF, aBCFlag, 1.0, vecx);
  for (unsigned int ip = 0; ip < np; ++ip) {
    aVelo[ip * 3 + 0] += vecx[ip * 3 + 0] / dt;
    aVelo[ip * 3 + 1] += vecx[ip * 3 + 1] / dt;
    aVelo[ip * 3 + 2] += vecx[ip * 3 + 2] / dt;
  }
}

int main() {
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aHex;
  dfm2::MeshHex3_Grid(
      aXYZ0, aHex,
      10, 10, 1, 0.2);
  std::vector<double> aMass(aXYZ0.size() / 3);

  dfm2::CMatrixSparse<double> smat;
  delfem2::CPreconditionerILU<double> silu;
  InitializeMatrix(
      smat, silu,
      aHex, aXYZ0);

  std::vector<double> aDisp(aXYZ0.size(), 0.0);
  std::vector<int> aBCFlag(aXYZ0.size(), 0); // 0: free, 1: fix BC
  {
    std::mt19937 rnddev(std::random_device{}());
    std::uniform_real_distribution<double> dist_m1p1(-1, 1);
    for (unsigned int ip = 0; ip < aXYZ0.size() / 3; ++ip) {
      double x0 = aXYZ0[ip * 3 + 0];
      if (x0 > 1.0e-10) { continue; }
      aBCFlag[ip * 3 + 0] = 1;
      aBCFlag[ip * 3 + 1] = 1;
      aBCFlag[ip * 3 + 2] = 1;
    }
  }
  std::vector<double> aVelo(aXYZ0.size(), 0.0);
  const double dt = 0.01;
  const double mass = 0.5;
  const double stiff_c1 = 1000.0;
  const double stiff_c2 = 1000.0;
  const double stiff_comp = 10000.0;
  const double gravity[3] = {0, 0, -10};

  // ----------------------
  delfem2::glfw::CViewer3 viewer;
  viewer.camera.view_height = 1.5;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::ZTOP;
  viewer.camera.theta = 0.1;
  viewer.camera.psi = 0.1;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();

  delfem2::opengl::setSomeLighting();
  while (!glfwWindowShouldClose(viewer.window)) {
    Simulation(
        aDisp, aVelo, smat, silu,
        aXYZ0, aHex, aBCFlag, dt, mass, stiff_c1, stiff_c2, stiff_comp, gravity);
    // -----
    viewer.DrawBegin_oldGL();
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0, 0, 0);
    delfem2::opengl::DrawMeshHex3D_EdgeDisp(
        aXYZ0.data(), aXYZ0.size() / 3,
        aHex.data(), aHex.size() / 8,
        aDisp.data());
    //
    ::glEnable(GL_LIGHTING);
//    dfm2::opengl::DrawMeshHex3D_FaceNorm(aXYZ0.data(), aHex.data(), aHex.size() / 8);
    delfem2::opengl::DrawHex3D_FaceNormDisp(aXYZ0, aHex, aDisp);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
