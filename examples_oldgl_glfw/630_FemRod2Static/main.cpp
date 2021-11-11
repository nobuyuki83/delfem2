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

#include "delfem2/ls_solver_block_sparse.h"
#include "delfem2/ls_pentadiagonal.h"
#include "delfem2/vec2.h"
#include "delfem2/mshuni.h"
#include "delfem2/femutil.h"
#include "delfem2/fem_rod2.h"
#include "delfem2/lsmats.h"
#include "delfem2/view_vectorx.h"
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

void Draw(
    const std::vector<double> &vtx_xy_initial,
    const std::vector<double> &vtx_xy){
  ::glLineWidth(1);
  ::glColor3d(0, 0, 0);
  DrawPolyline(vtx_xy_initial);
  ::glPointSize(5);
  delfem2::opengl::DrawPoints2d_Points(vtx_xy_initial);
  ::glColor3d(1, 0, 0);
  DrawPolyline(vtx_xy);
  delfem2::opengl::DrawPoints2d_Points(vtx_xy);
}

template <class BLOCK_LINEAR_SOLVER>
void StepTime(
    std::vector<double> &vtx_xy_deformed,
    const std::vector<double> &vtx_xy_initial,
    BLOCK_LINEAR_SOLVER &sparse) {
  const double stiff_stretch = 1.0;
  const double stiff_bend = 0.01;
  const size_t np = vtx_xy_initial.size() / 2;
  assert(np >= 3);
  assert(sparse.nblk() == np && sparse.ndim() == 2 );
  sparse.BeginMerge();
  double W = 0.0;
  for (unsigned int ihinge = 0; ihinge < np - 2; ++ihinge) {
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
      sparse.vec_r[aIP[ino] * 2 + 0] -= dWe[ino][0];
      sparse.vec_r[aIP[ino] * 2 + 1] -= dWe[ino][1];
    }
    sparse.template Merge<3, 3, 2, 2>(aIP, aIP, ddWe);
  }
  sparse.Solve();
  std::cout << W << std::endl;
  for (unsigned int ip = 0; ip < np; ++ip) {
    vtx_xy_deformed[ip * 2 + 0] += sparse.vec_x[ip * 2 + 0];
    vtx_xy_deformed[ip * 2 + 1] += sparse.vec_x[ip * 2 + 1];
  }
}

int main() {
  std::vector<double> vtx_xy_initial;
  for (unsigned int i = 0; i < 11; ++i) {
    vtx_xy_initial.push_back(i * 0.1);
    vtx_xy_initial.push_back(0.0);
  }
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> distm1p1(-1., +1.);

  dfm2::glfw::InitGLOld();
  dfm2::glfw::CViewer2 viewer;
  viewer.view_height = 0.5;
  viewer.trans[0] = -0.5;
  viewer.OpenWindow();

  while (true) {
    {
      dfm2::LinearSystemSolver_BlockSparse sparse;
      {
        std::vector<unsigned int> psup_ind, psup;
        delfem2::JArray_PSuP_Hair(
            psup_ind,psup,
            {0,(unsigned int) (vtx_xy_initial.size() / 2)});
        sparse.Initialize(
            vtx_xy_initial.size() / 2, 2,
            psup_ind, psup);
      }
      std::vector<double> vtx_xy = vtx_xy_initial;
      for (unsigned int ixy = 0; ixy < vtx_xy.size() / 2; ++ixy) {
        vtx_xy[ixy * 2 + 0] += distm1p1(rndeng);
        vtx_xy[ixy * 2 + 1] += distm1p1(rndeng);
      }
      for(unsigned int iframe=0;iframe<100;++iframe) {
        StepTime(
            vtx_xy,
            vtx_xy_initial, sparse);
        viewer.DrawBegin_oldGL();
        Draw(vtx_xy_initial, vtx_xy);
        viewer.SwapBuffers();
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) { break; }
      }
    }
    if (glfwWindowShouldClose(viewer.window)) { break; }
    {
      dfm2::LinearSystemSolver_BlockPentaDiagonal<2> sparse;
      sparse.Initialize(vtx_xy_initial.size() / 2);
      std::vector<double> vtx_xy = vtx_xy_initial;
      for (unsigned int ixy = 0; ixy < vtx_xy.size() / 2; ++ixy) {
        vtx_xy[ixy * 2 + 0] += distm1p1(rndeng);
        vtx_xy[ixy * 2 + 1] += distm1p1(rndeng);
      }
      for(unsigned int iframe=0;iframe<100;++iframe) {
        StepTime(
            vtx_xy,
            vtx_xy_initial, sparse);
        viewer.DrawBegin_oldGL();
        Draw(vtx_xy_initial, vtx_xy);
        viewer.SwapBuffers();
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) { break; }
      }
    }
  }


  viewer.ExitIfClosed();

  return 0;
}


