/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cstdlib>
#include <Eigen/Core>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should put before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/eigen/ls_dense.h"
#include "delfem2/eigen/ls_sparse.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/msh_primitive.h"
#include "delfem2/femsolidlinear.h"
#include "delfem2/femutil.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// --------------------------------------------------------------

void Simulation_Mat3(
    std::vector<double> &vtx_dispxyz,
    delfem2::CMatrixSparseBlock<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> &mA,
    //
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &hex_vtx,
    const std::vector<int> &vtx_bcflagxyz,
    //
    double mass,
    double myu,
    double lambda,
    const double gravity[3]) {
  const unsigned int np = vtx_xyz.size() / 3;
  const unsigned int nDoF = np * 3;
  mA.setZero();
  {
    double ddW[8][8][3][3];
    {
      const double gravity_zero[3] = {0, 0, 0};
      double aP0[8][3], aU[8][3];
      delfem2::FetchData<8, 3>(aP0, hex_vtx.data(), vtx_xyz.data());
      delfem2::FetchData<8, 3>(aU, hex_vtx.data(), vtx_dispxyz.data());
      double dW[8][3];
      std::fill_n(&dW[0][0], 8 * 3, 0.0);
      std::fill_n(&ddW[0][0][0][0], 8 * 8 * 3 * 3, 0.0);
      delfem2::elemMatRes_LinearSolidGravity3_Static_Q1(
          myu, lambda,
          0, gravity_zero,
          aP0, aU, ddW, dW);
    }
    std::vector<unsigned int> tmp_buffer;
    for (unsigned int ih = 0; ih < hex_vtx.size() / 8; ++ih) {
      const unsigned int *aIP = hex_vtx.data() + ih * 8;
      delfem2::Merge<8, 8, 3, 3, double>(mA, aIP, aIP, ddW, tmp_buffer);
    }
  }
  Eigen::VectorXd vec_b(nDoF);
  vec_b.setZero();
  for (unsigned int ip = 0; ip < np; ++ip) {
    vec_b[ip * 3 + 0] += mass * gravity[0];
    vec_b[ip * 3 + 1] += mass * gravity[1];
    vec_b[ip * 3 + 2] += mass * gravity[2];
  }
  { // comput rhs vectors
    const Eigen::VectorXd &vd = Eigen::Map<const Eigen::VectorXd>(vtx_dispxyz.data(), nDoF);
    AddMatVec(vec_b, 1.0, -1.0, mA, vd);
    std::cout << "energy" << vec_b.dot(vd) << std::endl;
  }
  SetFixedBC_Dia(mA, vtx_bcflagxyz.data(), 1.f);
  SetFixedBC_Col(mA, vtx_bcflagxyz.data());
  SetFixedBC_Row(mA, vtx_bcflagxyz.data());
  delfem2::setZero_Flag(vec_b, vtx_bcflagxyz, 0);
  // --------------------------------
  Eigen::VectorXd vec_x(vec_b.size());
  {
    double conv_ratio = 1.0e-6;
    int iteration = 1000;
    const std::size_t n = vec_b.size();
    Eigen::VectorXd tmp0(n), tmp1(n);
    std::vector<double> aConv = delfem2::Solve_CG(
        vec_b, vec_x, tmp0, tmp1,
        conv_ratio, iteration, mA);
    std::cout << aConv.size() << std::endl;
  }
  // ------------------------------
  dfm2::XPlusAY(
      vtx_dispxyz,
      vtx_bcflagxyz, 1.0, vec_x);
}

void Simulation_Mat4(
    std::vector<double> &vtx_dispxyz,
    delfem2::CMatrixSparseBlock<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>, 3> &mA,
    //
    const std::vector<double> &vtx_xyz,
    const std::vector<unsigned int> &hex_vtx,
    const std::vector<int> &vtx_bcflagxyz,
    //
    double mass,
    double myu,
    double lambda,
    const double gravity[3]) {
  const unsigned int np = vtx_xyz.size() / 3;
  mA.setZero();
  {
    double ddW[8][8][3][3];
    {
      const double gravity_zero[3] = {0, 0, 0};
      double aP0[8][3], aU[8][3];
      delfem2::FetchData<8, 3>(aP0, hex_vtx.data(), vtx_xyz.data());
      delfem2::FetchData<8, 3>(aU, hex_vtx.data(), vtx_dispxyz.data());
      double dW[8][3];
      std::fill_n(&dW[0][0], 8 * 3, 0.0);
      std::fill_n(&ddW[0][0][0][0], 8 * 8 * 3 * 3, 0.0);
      delfem2::elemMatRes_LinearSolidGravity3_Static_Q1(
          myu, lambda,
          0, gravity_zero,
          aP0, aU, ddW, dW);
    }
    std::vector<unsigned int> tmp_buffer;
    for (unsigned int ih = 0; ih < hex_vtx.size() / 8; ++ih) {
      const unsigned int *aIP = hex_vtx.data() + ih * 8;
      delfem2::Merge<8, 8, 3, 3, double>(mA, aIP, aIP, ddW, tmp_buffer);
    }
  }
  Eigen::Matrix<double, -1, 4, Eigen::RowMajor> vec_b(np, 4);
  vec_b.setZero();
  for (unsigned int ip = 0; ip < np; ++ip) {
    vec_b(ip, 0) += mass * gravity[0];
    vec_b(ip, 1) += mass * gravity[1];
    vec_b(ip, 2) += mass * gravity[2];
  }
  { // comput rhs vectors
    Eigen::Matrix<double, -1, 4, Eigen::RowMajor> vd(np, 4);
    for (unsigned int ip = 0; ip < np; ++ip) {
      vd(ip, 0) = vtx_dispxyz[ip * 3 + 0];
      vd(ip, 1) = vtx_dispxyz[ip * 3 + 1];
      vd(ip, 2) = vtx_dispxyz[ip * 3 + 2];
      vd(ip, 3) = 0.0;
    }
    AddMatVec(vec_b, 1.0, -1.0, mA, vd);
    std::cout << "energy" << delfem2::Dot(vec_b, vd) << std::endl;
  }
  SetFixedBC_Dia(mA, vtx_bcflagxyz.data(), 1.f);
  SetFixedBC_Col(mA, vtx_bcflagxyz.data());
  SetFixedBC_Row(mA, vtx_bcflagxyz.data());
  delfem2::setZero_Flag(vec_b, np, vtx_bcflagxyz, 0);
  // --------------------------------
  Eigen::Matrix<double, -1, 4, Eigen::RowMajor> vec_x(np, 4);
  {
    double conv_ratio = 1.0e-6;
    int iteration = 1000;
    Eigen::Matrix<double, -1, 4, Eigen::RowMajor> tmp0(np, 4), tmp1(np, 4);
    std::vector<double> aConv = delfem2::Solve_CG(
        vec_b, vec_x, tmp0, tmp1,
        conv_ratio, iteration, mA);
    std::cout << aConv.size() << std::endl;
  }
  // ------------------------------
  dfm2::XPlusAY(
      vtx_dispxyz,
      np, vtx_bcflagxyz, 1.0, vec_x);
}

int main() {
  std::vector<double> vtx_xyz;
  std::vector<unsigned int> hex_vtx;
  dfm2::MeshHex3_Grid(
      vtx_xyz, hex_vtx,
      20, 10, 10, 0.1);
  std::vector<double> aMass(vtx_xyz.size() / 3);

  delfem2::CMatrixSparseBlock<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> A3;
  delfem2::CMatrixSparseBlock<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>, 3> A4;
  {
    const unsigned int np = vtx_xyz.size() / 3;
    std::vector<unsigned int> psup_ind, psup;
    dfm2::JArray_PSuP_MeshElem(
        psup_ind, psup,
        hex_vtx.data(), hex_vtx.size() / 8, 8,
        vtx_xyz.size() / 3);
    A3.Initialize(np);
    A3.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
    A4.Initialize(np);
    A4.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
  }

  std::vector<double> vtx_dispxyz(vtx_xyz.size(), 0.0);
  std::vector<int> vtx_bcflagxyz(vtx_xyz.size(), 0.0); // 0: free, 1: fix BC
  {
    for (unsigned int ip = 0; ip < vtx_xyz.size() / 3; ++ip) {
      double x0 = vtx_xyz[ip * 3 + 0];
      if (x0 > 1.0e-10) { continue; }
      vtx_bcflagxyz[ip * 3 + 0] = 1;
      vtx_bcflagxyz[ip * 3 + 1] = 1;
      vtx_bcflagxyz[ip * 3 + 2] = 1;
    }
  }
  const double mass = 0.5;
  const double gravity[3] = {0, 0, -10};
  vtx_dispxyz.assign(vtx_xyz.size(), 0.0);
  for (unsigned int i = 0; i < vtx_dispxyz.size(); ++i) {
    if (vtx_bcflagxyz[i] != 0) { continue; }
    vtx_dispxyz[i] = (i % 10) * 1.0e-4;
  }
  Simulation_Mat3(
      vtx_dispxyz, A3,
      vtx_xyz, hex_vtx, vtx_bcflagxyz,
      mass, 1.0e+5, 1.e+5, gravity);

  vtx_dispxyz.assign(vtx_xyz.size(), 0.0);
  for (unsigned int i = 0; i < vtx_dispxyz.size(); ++i) {
    if (vtx_bcflagxyz[i] != 0) { continue; }
    vtx_dispxyz[i] = (i % 10) * 1.0e-4;
  }
  Simulation_Mat4(
      vtx_dispxyz, A4,
      vtx_xyz, hex_vtx, vtx_bcflagxyz,
      mass, 1.0e+5, 1.e+5, gravity);

  // ----------------------
  delfem2::glfw::CViewer3 viewer(1.5);
//  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::ZTOP;
//  viewer.camera.theta = 0.1;
//  viewer.camera.psi = 0.1;
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();

  delfem2::opengl::setSomeLighting();
  while (!glfwWindowShouldClose(viewer.window)) {
    // -----
    viewer.DrawBegin_oldGL();
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0, 0, 0);
    delfem2::opengl::DrawMeshHex3D_EdgeDisp(
        vtx_xyz.data(), vtx_xyz.size() / 3,
        hex_vtx.data(), hex_vtx.size() / 8,
        vtx_dispxyz.data());
    //
    ::glEnable(GL_LIGHTING);
//    dfm2::opengl::DrawMeshHex3D_FaceNorm(aXYZ0.data(), aHex.data(), aHex.size() / 8);
    delfem2::opengl::DrawHex3D_FaceNormDisp(vtx_xyz, hex_vtx, vtx_dispxyz);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
