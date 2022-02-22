/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
#include <Eigen/Core>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should put before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/msh_topology_uniform.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/dtri_topology.h"
#include "delfem2/eigen/ls_dense.h"
#include "delfem2/eigen/ls_sparse.h"
#include "delfem2/eigen/ls_ilu_sparse.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/fem_solidlinear.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

void MakeMesh(
    std::vector<double> &aXY1,
    std::vector<unsigned int> &aTri1,
    std::vector<int> &aBCFlag,
    unsigned int ndim) {
  std::vector<std::vector<double> > aaXY;
  const double len = 1.0;
  {
    aaXY.resize(1);
    aaXY[0].push_back(-len);
    aaXY[0].push_back(-len);
    aaXY[0].push_back(-len);
    aaXY[0].push_back(+len);
    aaXY[0].push_back(+len);
    aaXY[0].push_back(+len);
    aaXY[0].push_back(+len);
    aaXY[0].push_back(-len);
  }
  std::vector<delfem2::CDynPntSur> aPo2D;
  std::vector<delfem2::CDynTri> aETri;
  std::vector<delfem2::CVec2d> aVec2;
  delfem2::GenMesh(aPo2D, aETri, aVec2,
                   aaXY, 0.05, 0.05);
  MeshTri2D_Export(
      aXY1, aTri1,
      aVec2, aETri);
  const unsigned int np = aXY1.size() / 2;
  aBCFlag.assign(np * ndim, 0);
  for (unsigned int ip = 0; ip < np; ++ip) {
//    const double px = aXY1[ip*2+0];
    const double py = aXY1[ip * 2 + 1];
    if (fabs(py - len) > 0.0001) { continue; }
    for (unsigned int idim = 0; idim < ndim; ++idim) {
      aBCFlag[ip * 2 + idim] = 1;
    }
  }
  std::cout << "  ntri;" << aTri1.size() / 3 << "  nXY:" << aXY1.size() / 2 << std::endl;
}

void Solve1(
    std::vector<double> &aVal,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &aBCFlag) {
  const unsigned int np = aXY1.size() / 2;
  const unsigned int nDoF = np * 2;
  // -----------
  std::vector<unsigned int> psup_ind0, psup0;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind0, psup0,
      aTri1.data(), aTri1.size() / 3, 3,
      aXY1.size() / 2);
  // -------------
  delfem2::CMatrixSparseBlock<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d>> mA;
  mA.Initialize(np);
  mA.SetPattern(psup_ind0.data(), psup_ind0.size(), psup0.data(), psup0.size());
  // ----------------------
  double myu = 10.0;
  double lambda = 10.0;
  double rho = 1.0;
  double g_x = 0.0;
  double g_y = -3.0;
  mA.setZero();
  Eigen::VectorXd vec_b(nDoF);
  vec_b.setZero();
  dfm2::MergeLinSys_SolidLinear_Static_MeshTri2D(
      mA, vec_b.data(),
      myu, lambda, rho, g_x, g_y,
      aXY1.data(), aXY1.size() / 2,
      aTri1.data(), aTri1.size() / 3,
      aVal.data());
  SetFixedBC_Dia(mA, aBCFlag.data(), 1.f);
  SetFixedBC_Col(mA, aBCFlag.data());
  SetFixedBC_Row(mA, aBCFlag.data());
  delfem2::setZero_Flag(vec_b, aBCFlag, 0);
  // ---------------
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
//  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  // --------------
  {
    delfem2::CILU_SparseBlock<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d>> ilu;
    delfem2::ILU_SetPattern0(ilu, mA);
    delfem2::ILU_CopyValue(ilu, mA);
    delfem2::ILU_Decompose(ilu);
    Eigen::VectorXd vecX1(vec_b.size());
  }
  // --------------
  delfem2::XPlusAY(aVal,
                   aBCFlag,
                   1.0, vec_x);
}

int main() {
  std::vector<unsigned int> tri_vtx;
  std::vector<double> vtx_xy;
  std::vector<int> vtx_bcflagxy;
  MakeMesh(
      vtx_xy, tri_vtx, vtx_bcflagxy,
      2);
  // ---
  std::vector<double> vtx_dispxy;
  {
    const unsigned int np = vtx_xy.size() / 2;
    vtx_dispxy.assign(np * 2, 0.0);
    Solve1(vtx_dispxy, vtx_xy, tri_vtx, vtx_bcflagxy);
  }
  // --------
  dfm2::glfw::CViewer3 viewer(1.5);
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  // ---------
  while (!::glfwWindowShouldClose(viewer.window)) {
    viewer.DrawBegin_oldGL();
    delfem2::opengl::DrawMeshTri2D_FaceDisp2D(
        vtx_xy.data(), vtx_xy.size() / 2,
        tri_vtx.data(), tri_vtx.size() / 3,
        vtx_dispxy.data(), 2);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
}
