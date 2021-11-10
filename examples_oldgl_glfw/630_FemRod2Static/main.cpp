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
#include "delfem2/mshuni.h"
#include "delfem2/femutil.h"
#include "delfem2/fem_rod2.h"
#include "delfem2/lsmats.h"
#include "delfem2/view_vectorx.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/glfw/util.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

void DrawPolyline(const std::vector<double> &vtx_xy) {
  ::glBegin(GL_LINE_STRIP);
  for (unsigned int ixy = 0; ixy < vtx_xy.size() / 2; ++ixy) {
    ::glVertex2dv(vtx_xy.data() + ixy * 2);
  }
  ::glEnd();
}

void Solve(
    std::vector<double> &vtx_xy_deformed,
    const std::vector<double> &vtx_xy_initial,
    dfm2::CMatrixSparse<double> &sparse) {
  const double stiff_stretch = 1.0;
  const double stiff_bend = 0.01;
  const size_t np = vtx_xy_initial.size() / 2;
  assert(np >= 3);
  assert(sparse.ncolblk_ == np && sparse.nrowblk_ == np);
  assert(sparse.ncoldim_ == 2 && sparse.nrowdim_ == 2);
  unsigned int nhinge = np - 2;
  std::vector<unsigned int> merge_buffer;
  sparse.setZero();
  std::vector<double> vec_r(np * 2, 0.0);
  double W = 0.0;
  for (unsigned int ihinge = 0; ihinge < nhinge; ++ihinge) {
    const unsigned int aIP[3] = {ihinge, ihinge + 1, ihinge + 2};
    double aP[3][2];
    dfm2::FetchData<3, 2>(aP, aIP, vtx_xy_initial.data());
    double ap[3][2];
    dfm2::FetchData<3, 2>(ap, aIP, vtx_xy_deformed.data());
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
    dfm2::Merge<3, 3, 2, 2>(sparse, aIP, aIP, ddWe, merge_buffer);
  }
  std::vector<double> vec_x;
  vec_x.assign(np * 2, 0.0);
  {
    const std::size_t n = vec_r.size();
    std::vector<double> tmp0(n), tmp1(n);
    auto aConvHist = dfm2::Solve_CG(
        dfm2::ViewAsVectorXd(vec_r),
        dfm2::ViewAsVectorXd(vec_x),
        dfm2::ViewAsVectorXd(tmp0),
        dfm2::ViewAsVectorXd(tmp1),
        1.0e-4, 300, sparse);
    if (!aConvHist.empty()) {
      std::cout << "            conv: " << aConvHist.size();
      std::cout << " " << aConvHist[0];
      std::cout << " " << aConvHist[aConvHist.size() - 1] << std::endl;
    }
  }
  std::cout << W << std::endl;
  for (unsigned int ip = 0; ip < np; ++ip) {
    vtx_xy_deformed[ip * 2 + 0] += vec_x[ip * 2 + 0];
    vtx_xy_deformed[ip * 2 + 1] += vec_x[ip * 2 + 1];
  }
}

int main() {
  const unsigned int N = 11;
  std::vector<double> vtx_xy;  // vertices of polyline in the initial config
  for (unsigned int i = 0; i < N; ++i) {
    vtx_xy.push_back(i * 0.1);
    vtx_xy.push_back(0.0);
  }
  const std::vector<double> vtx_xy_initial = vtx_xy;
  std::mt19937 rndeng(std::random_device{}());
  dfm2::CMatrixSparse<double> sparse;
  {
    std::vector<unsigned int> psup_ind, psup;
    delfem2::JArray_PSuP_Hair(
        psup_ind,psup,
        {0,(unsigned int) (vtx_xy.size() / 2)});
    sparse.Initialize(vtx_xy.size() / 2, 2, true);
    sparse.SetPattern(
        psup_ind.data(), psup_ind.size(),
        psup.data(), psup.size());
  }

  dfm2::glfw::InitGLOld();
  dfm2::glfw::CViewer2 viewer;
  {
    viewer.view_height = 0.5;
    viewer.trans[0] = -0.5;
  }
  viewer.OpenWindow();

  int iframe = 0;
  while (true) {
    if (iframe == 0) {
      vtx_xy = vtx_xy_initial;
      std::uniform_real_distribution<double> distm1p1(-1., +1.);
      double mag = 1;
      for (unsigned int ixy = 0; ixy < vtx_xy.size() / 2; ++ixy) {
        vtx_xy[ixy * 2 + 0] += mag * distm1p1(rndeng);
        vtx_xy[ixy * 2 + 1] += mag * distm1p1(rndeng);
      }
    }
    Solve(
        vtx_xy,
        vtx_xy_initial, sparse);
    iframe = (iframe + 1) % 100;
    //
    viewer.DrawBegin_oldGL();
    ::glLineWidth(1);
    ::glColor3d(0, 0, 0);
    DrawPolyline(vtx_xy_initial);
    ::glPointSize(5);
    delfem2::opengl::DrawPoints2d_Points(vtx_xy_initial);
    //
    ::glColor3d(1, 0, 0);
    DrawPolyline(vtx_xy);
    delfem2::opengl::DrawPoints2d_Points(vtx_xy);
    viewer.SwapBuffers();
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { break; }
  }
  viewer.ExitIfClosed();

  return 0;
}


