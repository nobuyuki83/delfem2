/*
 * Copyright (c) 2020-2021 Nobuyuki Umetani
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

#include "delfem2/ls_solver_block_sparse.h"
#include "delfem2/ls_pentadiagonal.h"
#include "delfem2/mshuni.h"
#include "delfem2/fem_distance3.h"
#include "delfem2/fem_rod3_straight.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

// -------------------------------------

template <class BLOCK_LINEAR_SOLVER>
void OptimizeRod(
    std::vector<delfem2::CVec3d> &vtx_pos_deformed,
    double stiff_stretch,
    const std::vector<delfem2::CVec3d> &vtx_pos_initial,
    const std::vector<unsigned int> &hair_rootidx,
    BLOCK_LINEAR_SOLVER &mats) {
  namespace dfm2 = delfem2;
  mats.BeginMerge();
  double W = 0;
  for (unsigned int ihair = 0; ihair < hair_rootidx.size() - 1; ++ihair) {
    const unsigned int ips = hair_rootidx[ihair];  // index of point at the root
    const unsigned int ns = hair_rootidx[ihair + 1] - hair_rootidx[ihair] - 1;
    for (unsigned int is = 0; is < ns; ++is) {
      const unsigned int ip0 = ips + is + 0;
      const unsigned int ip1 = ips + is + 1;
      const double L0 = (vtx_pos_initial[ip0] - vtx_pos_initial[ip1]).norm();
      const dfm2::CVec3d aPE[2] = {vtx_pos_deformed[ip0], vtx_pos_deformed[ip1]};
      // --------------
      dfm2::CVec3d dW_dP[2];
      dfm2::CMat3d ddW_ddP[2][2];
      W += WdWddW_SquareLengthLineseg3D(
          dW_dP, ddW_ddP,
          stiff_stretch, aPE, L0);
      double ddW_ddP0[2][2][3][3];
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          ddW_ddP[i][j].CopyTo(&ddW_ddP0[i][j][0][0]);
        }
      }
      const unsigned int aINoel[2] = {ip0, ip1};
      mats.template Merge<2, 2, 3, 3>(aINoel, aINoel, ddW_ddP0);
      for (int in = 0; in < 2; in++) {
        const unsigned int ip = aINoel[in];
        mats.vec_r[ip * 3 + 0] -= dW_dP[in].x;
        mats.vec_r[ip * 3 + 1] -= dW_dP[in].y;
        mats.vec_r[ip * 3 + 2] -= dW_dP[in].z;
      }
    }
    for (unsigned int is = 0; is < ns - 1; ++is) {
      const unsigned int ip0 = ips + is + 0;
      const unsigned int ip1 = ips + is + 1;
      const unsigned int ip2 = ips + is + 2;
      const dfm2::CVec3d aPE[3] = {vtx_pos_deformed[ip0], vtx_pos_deformed[ip1], vtx_pos_deformed[ip2]};
      // --------------
      dfm2::CVec3d dW_dP[3];
      dfm2::CMat3d ddW_ddP[3][3];
      W += dfm2::WdWddW_Rod3BendStraight(
          dW_dP, ddW_ddP,
          aPE);
      double ddW_ddP0[3][3][3][3];
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          ddW_ddP[i][j].CopyTo(&ddW_ddP0[i][j][0][0]);
        }
      }
      const unsigned int aINoel[3] = {ip0, ip1, ip2};
      mats.template Merge<3, 3, 3, 3>(aINoel, aINoel, ddW_ddP0);
      for (int in = 0; in < 3; in++) {
        const unsigned int ip = aINoel[in];
        mats.vec_r[ip * 3 + 0] -= dW_dP[in].x;
        mats.vec_r[ip * 3 + 1] -= dW_dP[in].y;
        mats.vec_r[ip * 3 + 2] -= dW_dP[in].z;
      }
    }
  }
  std::cout << W << std::endl;
  if( W < 1.0e-10 ){ return; }
  mats.Solve();
  for (unsigned int ip = 0; ip < vtx_pos_deformed.size(); ++ip) {
    vtx_pos_deformed[ip].x += mats.vec_x[ip * 3 + 0];
    vtx_pos_deformed[ip].y += mats.vec_x[ip * 3 + 1];
    vtx_pos_deformed[ip].z += mats.vec_x[ip * 3 + 2];
  }
}

// -------------------------------

namespace dfm2 = delfem2;

void Draw(const std::vector<dfm2::CVec3d> &vec_pos) {
  ::glColor3d(0, 0, 0);
  //
  ::glBegin(GL_LINE_STRIP);
  for (const auto &p : vec_pos) {
    ::glVertex3dv(p.p);
  }
  ::glEnd();
  //
  ::glPointSize(3);
  ::glBegin(GL_POINTS);
  for (const auto &p : vec_pos) {
    ::glVertex3dv(p.p);
  }
  ::glEnd();
}

int main() {
  std::mt19937 reng(std::random_device{}());
  std::uniform_real_distribution<double> dist01(0.0, 1.0);
  //
  std::vector<dfm2::CVec3d> vtx_pos_ini;
  vtx_pos_ini.reserve(10);
  for (int ip = 0; ip < 10; ++ip) {
    vtx_pos_ini.emplace_back(ip * 0.1, 0., 0.);
  }
  const std::vector<unsigned int> vec_point_index_root = {
      0,
      static_cast<unsigned int>(vtx_pos_ini.size())};
  const double stiff_stretch = 10.0;
  // -----
  dfm2::glfw::CViewer3 viewer(1.5);
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  //
  while (!glfwWindowShouldClose(viewer.window)) {
    {  // solve using sparse linear solver
      dfm2::LinearSystemSolver_BlockSparse ls_solver;
      {
        std::vector<unsigned int> psup_ind, psup;
        delfem2::JArray_PSuP_Hair(
            psup_ind,psup,
            vec_point_index_root);
        ls_solver.Initialize(
            vtx_pos_ini.size(),3,
            psup_ind,psup);
      }
      std::vector<dfm2::CVec3d> vtx_pos = vtx_pos_ini;
      for (auto & vtx_po : vtx_pos) {
        vtx_po.x += dist01(reng) * 0.01;
        vtx_po.y += dist01(reng) * 0.01;
        vtx_po.z += dist01(reng) * 0.01;
      }
      for(unsigned int iframe=0;iframe<100;++iframe) {
        OptimizeRod(
            vtx_pos,
            stiff_stretch,
            vtx_pos_ini, vec_point_index_root,
            ls_solver);
        viewer.DrawBegin_oldGL();
        Draw(vtx_pos);
        viewer.SwapBuffers();
        glfwPollEvents();
      }
    }
    {  // solve using penta-diagonal linear solver
      dfm2::LinearSystemSolver_BlockPentaDiagonal<3> ls_solver;
      ls_solver.Initialize(vtx_pos_ini.size());
      std::vector<dfm2::CVec3d> vtx_pos = vtx_pos_ini;
      for (auto & vtx_po : vtx_pos) {
        vtx_po.x += dist01(reng) * 0.01;
        vtx_po.y += dist01(reng) * 0.01;
        vtx_po.z += dist01(reng) * 0.01;
      }
      for(unsigned int iframe=0;iframe<100;++iframe) {
        OptimizeRod(
            vtx_pos,
            stiff_stretch,
            vtx_pos_ini, vec_point_index_root,
            ls_solver);
        viewer.DrawBegin_oldGL();
        Draw(vtx_pos);
        viewer.SwapBuffers();
        glfwPollEvents();
      }
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
