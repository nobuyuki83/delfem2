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

#include "delfem2/lsitrsol.h"
#include "delfem2/lsmats.h"
#include "delfem2/lsvecx.h"
#include "delfem2/mshuni.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/fem_distance3.h"
#include "delfem2/fem_rod3_straight.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

// -------------------------------------

void OptimizeRod(
    std::vector<delfem2::CVec3d> &vec_pos,
    delfem2::CMatrixSparse<double> &mats,
    std::vector<double> &vec_r,
    std::vector<unsigned int> &tmp_buffer,
    double stiff_stretch,
    const std::vector<delfem2::CVec3d> &vec_pos_ini,
    const std::vector<unsigned int> &vec_point_index_root,
    const std::vector<int> &vec_flag_dof) {
  namespace dfm2 = delfem2;
  mats.setZero();
  vec_r.assign(vec_pos.size() * 3, 0.);
  double W = 0;
  for (unsigned int ihair = 0; ihair < vec_point_index_root.size() - 1; ++ihair) {
    const unsigned int ips = vec_point_index_root[ihair];  // index of point at the root
    // number of line segments
    const unsigned int ns = vec_point_index_root[ihair + 1] - vec_point_index_root[ihair] - 1;
    for (unsigned int is = 0; is < ns; ++is) {
      const unsigned int ip0 = ips + is + 0;
      const unsigned int ip1 = ips + is + 1;
      const double L0 = (vec_pos_ini[ip0] - vec_pos_ini[ip1]).norm();
      const dfm2::CVec3d aPE[2] = {vec_pos[ip0], vec_pos[ip1]};
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
      dfm2::Merge<2, 2, 3, 3, double>(
          mats,
          aINoel, aINoel, ddW_ddP0, tmp_buffer);
      for (int in = 0; in < 2; in++) {
        const unsigned int ip = aINoel[in];
        vec_r[ip * 3 + 0] -= dW_dP[in].x;
        vec_r[ip * 3 + 1] -= dW_dP[in].y;
        vec_r[ip * 3 + 2] -= dW_dP[in].z;
      }
    }
    for (unsigned int is = 0; is < ns - 1; ++is) {
      const unsigned int ip0 = ips + is + 0;
      const unsigned int ip1 = ips + is + 1;
      const unsigned int ip2 = ips + is + 2;
      const dfm2::CVec3d aPE[3] = {vec_pos[ip0], vec_pos[ip1], vec_pos[ip2]};
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
      dfm2::Merge<3, 3, 3, 3, double>(
          mats,
          aINoel, aINoel, ddW_ddP0, tmp_buffer);
      for (int in = 0; in < 3; in++) {
        const unsigned int ip = aINoel[in];
        vec_r[ip * 3 + 0] -= dW_dP[in].x;
        vec_r[ip * 3 + 1] -= dW_dP[in].y;
        vec_r[ip * 3 + 2] -= dW_dP[in].z;
      }
    }
  }
  mats.SetFixedBC(vec_flag_dof.data());
  dfm2::setRHS_Zero(vec_r, vec_flag_dof, 0);
  std::vector<double> vec_x;
  vec_x.assign(vec_pos.size() * 3, 0.0);
  {
    const std::size_t n = vec_r.size();
    std::vector<double> tmp0(n), tmp1(n);
    auto aConvHist = dfm2::Solve_CG(
        dfm2::CVecXd(vec_r),
        dfm2::CVecXd(vec_x),
        dfm2::CVecXd(tmp0),
        dfm2::CVecXd(tmp1),
        1.0e-4, 300, mats);
    if (!aConvHist.empty()) {
      std::cout << "            conv: ";
      std::cout << aConvHist.size() << " ";
      std::cout << aConvHist[0] << " ";
      std::cout << aConvHist[aConvHist.size() - 1] << std::endl;
    }
  }
  for (unsigned int ip = 0; ip < vec_pos.size(); ++ip) {
    vec_pos[ip].x += vec_x[ip * 3 + 0];
    vec_pos[ip].y += vec_x[ip * 3 + 1];
    vec_pos[ip].z += vec_x[ip * 3 + 2];
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
  std::vector<dfm2::CVec3d> vec_pos;
  vec_pos.reserve(10);
  for (int ip = 0; ip < 10; ++ip) {
    vec_pos.emplace_back(ip * 0.1, 0., 0.);
  }
  const std::vector<dfm2::CVec3d> vec_pos_ini = vec_pos;
  std::vector<int> vec_flag_dof(vec_pos.size() * 3, 0);
  vec_flag_dof[0] = 1;
  vec_flag_dof[1] = 1;
  vec_flag_dof[2] = 1;
  for (unsigned int ip = 0; ip < vec_pos.size(); ++ip) {
    if (vec_flag_dof[ip * 3 + 0] != 0) continue;
    vec_pos[ip].x += dist01(reng) * 0.01;
    vec_pos[ip].y += dist01(reng) * 0.01;
    vec_pos[ip].z += dist01(reng) * 0.01;
  }
  const std::vector<unsigned int> vec_point_index_root = {
      0,
      static_cast<unsigned int>(vec_pos_ini.size())};
  const double stiff_stretch = 10.0;
  //
  dfm2::CMatrixSparse<double> mats;
  {
    std::vector<unsigned int> psup_ind, psup;
    delfem2::JArray_PSuP_Hair(
        psup_ind,psup,
        vec_point_index_root);
    mats.Initialize(vec_pos_ini.size(),3,true);
    mats.SetPattern(
        psup_ind.data(), psup_ind.size(),
        psup.data(), psup.size());
  }
  std::vector<double> vec_r(vec_pos.size() * 3);
  std::vector<unsigned int> tmp_buffer;
  // -----
  dfm2::glfw::CViewer3 viewer(1.5);
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  //
  while (!glfwWindowShouldClose(viewer.window)) {
    {
      OptimizeRod(
          vec_pos, mats, vec_r, tmp_buffer,
          stiff_stretch,
          vec_pos_ini, vec_point_index_root, vec_flag_dof);
      static int iframe = 0;
      iframe += 1;
      if( iframe >= 100 ){
        for (unsigned int ip = 0; ip < vec_pos.size(); ++ip) {
          if (vec_flag_dof[ip * 3 + 0] != 0) continue;
          vec_pos[ip].x += dist01(reng) * 0.01;
          vec_pos[ip].y += dist01(reng) * 0.01;
          vec_pos[ip].z += dist01(reng) * 0.01;
        }
        iframe = 0;
      }
    }
    // --------
    viewer.DrawBegin_oldGL();
    Draw(vec_pos);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
