/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/vec2.h"
#include "delfem2/femutil.h"
#include "delfem2/femrod.h"
#include "delfem2/lsvecx.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/glfw/util.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

void DrawPolyline(const std::vector<double> &aXY) {
  ::glBegin(GL_LINE_STRIP);
  for (unsigned int ixy = 0; ixy < aXY.size() / 2; ++ixy) {
    ::glVertex2dv(aXY.data() + ixy * 2);
  }
  ::glEnd();
}

void Solve(
    std::vector<double> &axy,
    std::vector<double> &auv,
    const std::vector<double> &aXY,
    const std::vector<int> &aBCFlag,
    double dt,
    double stiff_stretch,
    double stiff_bend,
    double mass_point,
    const double gravity[2],
    dfm2::CMatrixSparse<double> &mats) {
  unsigned int np = static_cast<unsigned int>(aXY.size() / 2);
  assert(np >= 3);
  assert(mats.ncolblk_ == np && mats.nrowblk_ == np);
  assert(mats.ncoldim_ == 2 && mats.nrowdim_ == 2);
  for (unsigned int ip = 0; ip < np; ++ip) {
    axy[ip * 2 + 0] += dt * auv[ip * 2 + 0];
    axy[ip * 2 + 1] += dt * auv[ip * 2 + 1];
  }
  unsigned int nhinge = np - 2;
  std::vector<unsigned int> merge_buffer;
  mats.setZero();
  std::vector<double> vec_r(np * 2, 0.0);
  double W = 0.0;
  for (unsigned int ihinge = 0; ihinge < nhinge; ++ihinge) {
    const unsigned int aIP[3] = {ihinge, ihinge + 1, ihinge + 2};
    double aP[3][2];
    dfm2::FetchData<3, 2>(aP, aIP, aXY.data());
    double ap[3][2];
    dfm2::FetchData<3, 2>(ap, aIP, axy.data());
    const double aL[2] = {
        dfm2::Distance2(aP[0], aP[1]),
        dfm2::Distance2(aP[1], aP[2])};
    double We, dWe[3][2], ddWe[3][3][2][2];
    dfm2::WdWddW_Rod2(
        We, dWe, ddWe,
        ap, aL, stiff_stretch, stiff_stretch, stiff_bend);
    W += We;
    for (int ino = 0; ino < 3; ++ino) {
      vec_r[aIP[ino] * 2 + 0] -= dWe[ino][0];
      vec_r[aIP[ino] * 2 + 1] -= dWe[ino][1];
    }
    dfm2::Merge<3, 3, 2, 2>(mats, aIP, aIP, ddWe, merge_buffer);
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    mats.val_dia_[ip * 4 + 0] += mass_point / (dt * dt);
    mats.val_dia_[ip * 4 + 3] += mass_point / (dt * dt);
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    vec_r[ip * 2 + 0] += mass_point * gravity[0];
    vec_r[ip * 2 + 1] += mass_point * gravity[1];
  }
  assert(aBCFlag.size() == np * 2);
  mats.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_r, aBCFlag, 0);
  std::vector<double> vec_x;
  vec_x.assign(np * 2, 0.0);
  {
    const std::size_t n = vec_r.size();
    std::vector<double> tmp0(n), tmp1(n);
    auto aConvHist = Solve_CG(
        dfm2::CVecXd(vec_r),
        dfm2::CVecXd(vec_x),
        dfm2::CVecXd(tmp0),
        dfm2::CVecXd(tmp1),
        1.0e-4, 300, mats);
    if (!aConvHist.empty()) {
      std::cout << "            conv: " << aConvHist.size();
      std::cout << " " << aConvHist[0];
      std::cout << " " << aConvHist[aConvHist.size() - 1] << std::endl;
    }
  }
  std::cout << W << std::endl;
  for (unsigned int ip = 0; ip < np; ++ip) {
    axy[ip * 2 + 0] += vec_x[ip * 2 + 0];
    axy[ip * 2 + 1] += vec_x[ip * 2 + 1];
    auv[ip * 2 + 0] += vec_x[ip * 2 + 0] / dt;
    auv[ip * 2 + 1] += vec_x[ip * 2 + 1] / dt;
  }
}

int main(
    [[maybe_unused]] int argc,
    [[maybe_unused]] char *argv[]) {
  const unsigned int N = 11;
  const double dt = 1.0 / 60.0;
  const double stiff_stretch = 10.0;
  const double stiff_bend = 0.001;
  const double mass_point = 1.0e-5;
  const double gravity[2] = {0, -10};
  std::vector<double> aXY;  // vertices of polyline in the initial config
  for (unsigned int i = 0; i < N; ++i) {
    aXY.push_back(i * 0.1);
    aXY.push_back(0.0);
  }
  std::vector<double> aUV(aXY.size(), 0.);
  std::vector<int> aBCFlag(aXY.size(), 0);
  const std::vector<double> aXY0 = aXY;
  std::mt19937 rndeng(std::random_device{}());
  dfm2::CMatrixSparse<double> mats;
  dfm2::MakeSparseMatrix_RodHair(
      mats, {0, (unsigned int) (aXY.size() / 2)}, 2);

  // opengl set up
  dfm2::glfw::InitGLOld();
  dfm2::glfw::CViewer2 viewer;
  {
    viewer.view_height = 0.7f;
    viewer.trans[0] = -0.2f;
    viewer.trans[1] = +0.3f;
    viewer.title = "dynamic 2D rod demo";
  }
  viewer.InitGL();  // opengl start here

  int iframe = 0;
  bool isNormalFix = false;
  while (true) {
    if (iframe == 0) {
      isNormalFix = !isNormalFix;
      // set boundary condition flag
      aBCFlag.assign(aXY.size(), 0);
      aBCFlag[0] = 1;
      aBCFlag[1] = 1;
      if (isNormalFix) {
        aBCFlag[2] = 1;
        aBCFlag[3] = 1;
      }
      // set initial coordinates
      aXY = aXY0;
      std::uniform_real_distribution<double> distm1p1(-1., +1.);
      double mag = 0.1;
      for (unsigned int i = 0; i < aXY.size(); ++i) {
        if (aBCFlag[i] != 0) { continue; }
        aXY[i] += mag * distm1p1(rndeng);
      }
      aUV.assign(aXY.size(), 0.0);
    }
    Solve(
        aXY, aUV,
        aXY0, aBCFlag, dt, stiff_stretch, stiff_bend, mass_point, gravity,
        mats);
    iframe = (iframe + 1) % 100;
    //
    viewer.DrawBegin_oldGL();
    ::glLineWidth(1);
    ::glColor3d(0, 0, 0);
    DrawPolyline(aXY0);
    ::glPointSize(5);
    delfem2::opengl::DrawPoints2d_Points(aXY0);
    //
    ::glColor3d(1, 0, 0);
    DrawPolyline(aXY);
    delfem2::opengl::DrawPoints2d_Points(aXY);
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  viewer.ExitIfClosed();

  return 0;
}


