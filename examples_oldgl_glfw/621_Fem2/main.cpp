/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/view_vectorx.h"
#include "delfem2/ls_ilu_block_sparse.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/ls_block_sparse.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/fem_navierstokes.h"
#include "delfem2/fem_stokes.h"
#include "delfem2/fem_solidlinear.h"
#include "delfem2/fem_poisson.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/dtri_topology.h"
#include "delfem2/color.h"
#include "delfem2/jagarray.h"
#include "delfem2/geo_polygon2.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// --------------------------------------------------------------

double AreaCGCurve(const std::vector<double> &aCV, double cg[2]) {
  const unsigned int nCV = (int) aCV.size() / 2;
  double area = 0;
  cg[0] = 0;
  cg[1] = 0;
  for (unsigned int idiv = 0; idiv < nCV; idiv++) {
    unsigned int ipo0 = idiv;
    unsigned int ipo1 = idiv + 1;
    if (idiv == nCV - 1) { ipo1 = 0; }
    const double p0[2] = {aCV[ipo0 * 2 + 0], aCV[ipo0 * 2 + 1]};
    const double p1[2] = {aCV[ipo1 * 2 + 0], aCV[ipo1 * 2 + 1]};
    const double p2[2] = {0, 0};
    double a0 = dfm2::Area_Tri2(p0, p1, p2);
    double cg0[2] = {(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0};
    cg[0] += cg0[0] * a0;
    cg[1] += cg0[1] * a0;
    area += a0;
  }
  cg[0] /= area;
  cg[1] /= area;
  return area;
}

void MakeRandomCV(unsigned int nCV, std::vector<double> &aCV) {
  std::random_device rd;
  std::mt19937 rdeng(rd());
  std::uniform_real_distribution<double> dist(0, 1.0);
  aCV.clear();
  for (unsigned int icv = 0; icv < nCV; icv++) {
    /*
     {
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     }
     */
    {
      double tht = icv * 3.1415 * 2.0 / nCV;
      double r = dist(rdeng);
      double px = r * sin(tht);
      double py = r * cos(tht);
      aCV.push_back(px);
      aCV.push_back(py);
    }
  }

  {
    double cnt[2];
    double area = AreaCGCurve(aCV, cnt);
    if (area < 0) {
      std::vector<double> aCV0;
      for (unsigned int idiv = 0; idiv < nCV; idiv++) {
        unsigned int idiv0 = nCV - 1 - idiv;
        aCV0.push_back(aCV[idiv0 * 2 + 0]);
        aCV0.push_back(aCV[idiv0 * 2 + 1]);
      }
      aCV = aCV0;
      area = -area;
    }
    if (area < 1.0e-5) {
      aCV.clear();
      return;
    }
  }
}

void MakeCurveSpline(
    const std::vector<double> &aCV,
    std::vector<double> &aVecCurve,
    unsigned int ndiv = 5) {
  aVecCurve.resize(0);
  const size_t nCV = aCV.size() / 2;
  for (unsigned int icv = 0; icv < nCV; icv++) {
    const unsigned int icv0 = (icv + 0) % nCV;
    const unsigned int icv1 = (icv + 1) % nCV;
    const unsigned int icv2 = (icv + 2) % nCV;
    const double p0[2] = {aCV[icv0 * 2 + 0], aCV[icv0 * 2 + 1]};
    const double p1[2] = {aCV[icv1 * 2 + 0], aCV[icv1 * 2 + 1]};
    const double p2[2] = {aCV[icv2 * 2 + 0], aCV[icv2 * 2 + 1]};
    for (unsigned int idiv = 0; idiv < ndiv; idiv++) {
      const double t = 1.0 - (double) idiv / ndiv;
      const double w[3] = {0.5 * t * t, -t * t + t + 0.5, 0.5 * (1 - t) * (1 - t)};
      const double px = w[0] * p0[0] + w[1] * p1[0] + w[2] * p2[0];
      const double py = w[0] * p0[1] + w[1] * p1[1] + w[2] * p2[1];
      aVecCurve.push_back(px);
      aVecCurve.push_back(py);
    }
  }
}

// ------------------------------------------

void MakeMesh(
    std::vector<double> &aXY1,
    std::vector<unsigned int> &aTri1,
    std::vector<int> &loopIP_ind,
    std::vector<int> &loopIP,
    double len) {
  std::vector<double> aCV0;
  MakeRandomCV(8, aCV0); // current cv
  std::vector<double> aXY0;
  MakeCurveSpline(aCV0, aXY0, 20); // current curve
  std::vector<std::vector<double> > aaXY;
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
  aaXY.push_back(aXY0);
  // ---------------------------------
  std::vector<dfm2::CVec2d> aVec2;
  constexpr double elen = 0.05;
  {
    JArray_FromVecVec_XY(loopIP_ind, loopIP, aVec2,
                         aaXY);
    if (!CheckInputBoundaryForTriangulation(loopIP_ind, aVec2)) {
      return;
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind, aVec2);
    if (elen > 10e-10) {
      ResamplingLoop(loopIP_ind, loopIP, aVec2,
                     elen);
    }
  }
  {
    std::vector<dfm2::CDynPntSur> aPo2D;
    std::vector<dfm2::CDynTri> aETri;
    Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                   loopIP_ind, loopIP);
    if constexpr (elen > 1.0e-10) {
      dfm2::CInputTriangulation_Uniform param(1.0);
      std::vector<unsigned int> aFlgPnt(aPo2D.size());
      std::vector<unsigned int> aFlgTri(aETri.size(), 0);
      MeshingInside(
          aPo2D, aETri, aVec2, aFlgPnt, aFlgTri,
          aVec2.size(), 0, elen, param);
    }
    MeshTri2D_Export(
        aXY1, aTri1,
        aVec2, aETri);
  }
  std::cout << "  ntri;" << aTri1.size() / 3 << "  nXY:" << aXY1.size() / 2 << std::endl;
}

// ----------------------------------
// iproblem: 0, 1
void InitializeProblem_Scalar(
    std::vector<int> &aBCFlag,
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    double len) {
  const size_t np = aXY1.size() / 2;
  aBCFlag.assign(np, 0);
  for (unsigned int ip = 0; ip < np; ++ip) {
    const double px = aXY1[ip * 2 + 0];
    const double py = aXY1[ip * 2 + 1];
    if (fabs(fabs(px) - len) < 0.0001 || fabs(fabs(py) - len) < 0.0001) {
      aBCFlag[ip] = 1;
    }
  }
  //
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind, psup,
      aTri1.data(), aTri1.size() / 3, 3,
      (int) aXY1.size() / 2);
  dfm2::JArray_Sort(psup_ind, psup);
  mat_A.Initialize(np, 1, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
  ilu_A.SetPattern0(mat_A);
}

// -----------------------------
// iproblem: 0
void SolveProblem_Poisson(
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    std::vector<double>& aVal,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &aBCFlag) {
  const size_t np = aXY1.size() / 2;
  const size_t nDoF = np;
  // -----------------------
  const double alpha = 1.0;
  const double source = 1.0;
  mat_A.setZero();
  std::vector<double> vec_b;
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_Poission_MeshTri2D(
      mat_A, vec_b.data(),
      alpha, source,
      aXY1.data(), aXY1.size() / 2,
      aTri1.data(), aTri1.size() / 3,
      aVal.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag, 0);
  // ------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.CopyValue(mat_A);
  ilu_A.Decompose();
  vec_x.resize(vec_b.size());
  /*
  Solve_PCG(vec_b.data(),vec_x.data(),
            conv_ratio,iteration, mat_A,ilu_A);
   */
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    Solve_CG(
        dfm2::ViewAsVectorXd(vec_b),
		dfm2::ViewAsVectorXd(vec_x), 
		dfm2::ViewAsVectorXd(tmp0), 
		dfm2::ViewAsVectorXd(tmp1),
        conv_ratio, iteration, mat_A);
  }
  // ----------------
  dfm2::XPlusAY(aVal, nDoF, aBCFlag,
                1.0, vec_x);
}

// -------------------------------
// iproblem: 1
void SolveProblem_Diffusion(
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    std::vector<double> &aVal,
    std::vector<double> &aVelo,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &aBCFlag) {
  const size_t np = aXY1.size() / 2;
  const size_t nDoF = np;
  // ------------------
  const double alpha = 1.0;
  const double rho = 1.0;
  const double source = 1.0;
  double dt_timestep = 0.01;
  double gamma_newmark = 0.6;
  mat_A.setZero();
  std::vector<double> vec_b;
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_Diffusion_MeshTri2D(
      mat_A, vec_b.data(),
      alpha, rho, source,
      dt_timestep, gamma_newmark,
      aXY1.data(), aXY1.size() / 2,
      aTri1.data(), aTri1.size() / 3,
      aVal.data(), aVelo.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag, 0);
  // ------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.CopyValue(mat_A);
  ilu_A.Decompose();
  vec_x.resize(vec_b.size());
  {
    const size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    Solve_PCG(
        dfm2::ViewAsVectorXd(vec_b),
        dfm2::ViewAsVectorXd(vec_x),
        dfm2::ViewAsVectorXd(tmp0),
        dfm2::ViewAsVectorXd(tmp1),
        conv_ratio, iteration, mat_A, ilu_A);
  }
//  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  // -----------------
  dfm2::XPlusAYBZ(
      aVal, nDoF, aBCFlag,
      dt_timestep * gamma_newmark, vec_x,
      dt_timestep, aVelo);
  dfm2::XPlusAY(
      aVelo, nDoF, aBCFlag,
      1.0, vec_x);
}

void DrawScalar(
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<double> &aVal) {
  ::glDisable(GL_LIGHTING);
  {
    std::vector<std::pair<double, delfem2::CColor> > colorMap;
    ColorMap_BlueGrayRed(colorMap, 0, +0.1f);
    delfem2::opengl::DrawMeshTri2D_ScalarP1(
        aXY1.data(), aXY1.size() / 2,
        aTri1.data(), aTri1.size() / 3,
        aVal.data(), 1, colorMap);
  }
  ::glColor3d(0, 0, 0);
  delfem2::opengl::DrawMeshTri2D_Edge(aTri1, aXY1);
  ::glPointSize(2);
  ::glColor3d(0, 0, 0);
  delfem2::opengl::DrawPoints2d_Points(aXY1);
}

void ProblemScalar(
    dfm2::glfw::CViewer3 &viewer,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    double len) {
  std::vector<double> aVal;
  const size_t np = aXY1.size() / 2;
  dfm2::CMatrixSparse<double> mat_A;
  dfm2::CPreconditionerILU<double> ilu_A;
  std::vector<int> aBCFlag;
  // ------------
  glfwSetWindowTitle(viewer.window, "Poisson");
  aVal.assign(np, 0.0);
  InitializeProblem_Scalar(
      aBCFlag, mat_A, ilu_A,
      aXY1, aTri1, len);
  SolveProblem_Poisson(
      mat_A, ilu_A, aVal,
      aXY1, aTri1, aBCFlag);
  for (unsigned int iframe = 0; iframe < 50; ++iframe) { // poisson
    viewer.DrawBegin_oldGL();
    DrawScalar(aXY1, aTri1, aVal);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  // ---------------------------
  glfwSetWindowTitle(viewer.window, "Diffusion");
  aVal.assign(np, 0.0);
  std::vector<double> aVelo;
  aVelo.assign(np, 0.0);
  SolveProblem_Diffusion(
      mat_A, ilu_A, aVal, aVelo,
      aXY1, aTri1, aBCFlag);
  for (unsigned int iframe = 0; iframe < 100; ++iframe) {
    SolveProblem_Diffusion(
        mat_A, ilu_A, aVal, aVelo,
        aXY1, aTri1, aBCFlag);
    // -------
    viewer.DrawBegin_oldGL();
    DrawScalar(aXY1, aTri1, aVal);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
}

// ---------------------------------------------------
// iproblem: 2, 3
void InitializeProblem_Solid(
    std::vector<int> &aBCFlag,
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    double len) {
  const size_t np = aXY1.size() / 2;
  const size_t nDoF = np * 2;
  // ----------------
  aBCFlag.assign(nDoF, 0);
  for (unsigned int ip = 0; ip < np; ++ip) {
//    const double px = aXY1[ip*2+0];
    const double py = aXY1[ip * 2 + 1];
    if (fabs(py - len) > 0.0001) { continue; }
    aBCFlag[ip * 2 + 0] = 1;
    aBCFlag[ip * 2 + 1] = 1;
  }
  // -----------
  std::vector<unsigned int> psup_ind0, psup0;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind0, psup0,
      aTri1.data(), aTri1.size() / 3, 3,
      aXY1.size() / 2);
  // -------------
  mat_A.Initialize(np, 2, true);
  mat_A.SetPattern(psup_ind0.data(), psup_ind0.size(), psup0.data(), psup0.size());
  ilu_A.SetPattern0(mat_A);
}

// -------------------------
// iproblem: 2
void SolveProblem_LinearSolid_Static(
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    std::vector<double> &aVal,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &aBCFlag) {
  const auto np = static_cast<unsigned int>(aXY1.size() / 2);
  const unsigned int nDoF = np * 2;
  // ----------------------
  double myu = 10.0;
  double lambda = 10.0;
  double rho = 1.0;
  double g_x = 0.0;
  double g_y = -3.0;
  mat_A.setZero();
  std::vector<double> vec_b;
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_SolidLinear_Static_MeshTri2D(
      mat_A, vec_b.data(),
      myu, lambda, rho, g_x, g_y,
      aXY1.data(), aXY1.size() / 2,
      aTri1.data(), aTri1.size() / 3,
      aVal.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag, 0);
  // ---------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.CopyValue(mat_A);
  ilu_A.Decompose();
  vec_x.resize(vec_b.size());
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    Solve_PCG(
        dfm2::ViewAsVectorXd(vec_b), dfm2::ViewAsVectorXd(vec_x), dfm2::ViewAsVectorXd(tmp0), dfm2::ViewAsVectorXd(tmp1),
        conv_ratio, iteration, mat_A, ilu_A);
  }
//  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  // --------------
  dfm2::XPlusAY(aVal, nDoF, aBCFlag,
                1.0, vec_x);
}

// -----------------------------------
// iproblem: 3
void SolveProblem_LinearSolid_Dynamic(
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    std::vector<double> &aVal,
    std::vector<double> &aVelo,
    std::vector<double> &aAcc,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &aBCFlag) {
  const auto np = static_cast<unsigned int>(aXY1.size() / 2);
  const unsigned int nDoF = np * 2;
  // ------------------
  double beta_newmark = 0.36;
  double myu = 10.0;
  double lambda = 10.0;
  double rho = 1.0;
  double g_x = 0.0;
  double g_y = -3.0;
  double dt_timestep = 0.01;
  double gamma_newmark = 0.6;
  mat_A.setZero();
  std::vector<double> vec_b;
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_SolidLinear_NewmarkBeta_MeshTri2D(
      mat_A, vec_b.data(),
      myu, lambda, rho, g_x, g_y,
      dt_timestep, gamma_newmark, beta_newmark,
      aXY1.data(), aXY1.size() / 2,
      aTri1.data(), aTri1.size() / 3,
      aVal.data(), aVelo.data(), aAcc.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag, 0);
  // -----------------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.CopyValue(mat_A);
  ilu_A.Decompose();
  vec_x.resize(vec_b.size());
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    Solve_PCG(
        dfm2::ViewAsVectorXd(vec_b), dfm2::ViewAsVectorXd(vec_x), dfm2::ViewAsVectorXd(tmp0), dfm2::ViewAsVectorXd(tmp1),
        conv_ratio, iteration, mat_A, ilu_A);
  }
//  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  // -----------------------------
  dfm2::XPlusAYBZCW(aVal, nDoF, aBCFlag,
                    dt_timestep, aVelo,
                    0.5 * dt_timestep * dt_timestep, aAcc,
                    dt_timestep * dt_timestep * beta_newmark, vec_x);
  dfm2::XPlusAYBZ(aVelo, nDoF, aBCFlag,
                  dt_timestep * gamma_newmark, vec_x,
                  dt_timestep, aAcc);
  dfm2::XPlusAY(aAcc, nDoF, aBCFlag,
                1.0, vec_x);
}

void ProblemSolid(
    dfm2::glfw::CViewer3 &viewer,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    double len) {
  std::vector<double> aVal;
  const auto np = static_cast<unsigned int>(aXY1.size() / 2);
  dfm2::CMatrixSparse<double> mat_A;
  dfm2::CPreconditionerILU<double> ilu_A;
  std::vector<int> aBCFlag;
  glfwSetWindowTitle(viewer.window, "Linear Elastic Static");
  aVal.assign(np * 2, 0.0);
  InitializeProblem_Solid(
      aBCFlag, mat_A, ilu_A,
      aXY1, aTri1, len);
  SolveProblem_LinearSolid_Static(
      mat_A, ilu_A, aVal,
      aXY1, aTri1, aBCFlag);
  for (unsigned int iframe = 0; iframe < 50; ++iframe) {
    viewer.DrawBegin_oldGL();
    delfem2::opengl::DrawMeshTri2D_FaceDisp2D(
        aXY1.data(), aXY1.size() / 2,
        aTri1.data(), aTri1.size() / 3,
        aVal.data(), 2);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  // ---------------------------
  glfwSetWindowTitle(viewer.window, "Linear Elastic Dynamic");
  aVal.assign(np * 2, 0.0);
  std::vector<double> aVelo;
  aVelo.assign(np * 2, 0.0);
  std::vector<double> aAcc;
  aAcc.assign(np * 2, 0.0);
  InitializeProblem_Solid(
      aBCFlag, mat_A, ilu_A,
      aXY1, aTri1, len);
  for (unsigned int iframe = 0; iframe < 100; ++iframe) {
    SolveProblem_LinearSolid_Dynamic(
        mat_A, ilu_A, aVal, aVelo, aAcc,
        aXY1, aTri1, aBCFlag);
    viewer.DrawBegin_oldGL();
    delfem2::opengl::DrawMeshTri2D_FaceDisp2D(
        aXY1.data(), aXY1.size() / 2,
        aTri1.data(), aTri1.size() / 3,
        aVal.data(), 2);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
}

// -----------------------------
// iproblem: 4, 5, 6
void InitializeProblem_Fluid(
    std::vector<int> &aBCFlag,
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    std::vector<double> &aVal,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &loopIP_ind,
    const std::vector<int> &loopIP,
    double len) {
  // set boundary condition
  const int np = (int) aXY1.size() / 2;
  const int nDoF = np * 3;
  // -----------
  aBCFlag.assign(nDoF, 0);
  for (int ip = 0; ip < np; ++ip) {
    const double px = aXY1[ip * 2 + 0];
    const double py = aXY1[ip * 2 + 1];
    if (fabs(fabs(px) - len) < 0.0001 || fabs(fabs(py) - len) < 0.0001) {
      aBCFlag[ip * 3 + 0] = 1;
      aBCFlag[ip * 3 + 1] = 1;
    }
  }
  for (int iip = loopIP_ind[1]; iip < loopIP_ind[2]; ++iip) {
    int ip0 = loopIP[iip];
    aBCFlag[ip0 * 3 + 0] = 1;
    aBCFlag[ip0 * 3 + 1] = 1;
  }
  aBCFlag[0 * 3 + 2] = 1;
  // -------
  for (int ip = 0; ip < np; ++ip) {
//    const double px = aXY1[ip*2+0];
    const double py = aXY1[ip * 2 + 1];
    if (fabs(py - len) < 0.0001 && aBCFlag[ip * 3 + 0] == 1) {
      aVal[ip * 3 + 0] = 10;
    }
  }
  //
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind, psup,
      aTri1.data(), aTri1.size() / 3, 3,
      aXY1.size() / 2);
  dfm2::JArray_Sort(psup_ind, psup);
  /*
   CJaggedArray crs;
   crs.SetEdgeOfElem(aTri1, (int)aTri1.size()/3, 3, (int)aXY1.size()/2, false);
   crs.Sort();
   */
  ////
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
  ilu_A.SetPattern0(mat_A);
  //  ilu_A.Initialize_ILUk(mat_A, 5);
}

// ---------------------------
// iproblem: 7,8,9
void InitializeProblem_Fluid2(
    std::vector<int> &aBCFlag,
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    std::vector<double> &aVal,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &loopIP_ind,
    const std::vector<int> &loopIP,
    double len) {
  const size_t np = aXY1.size() / 2;
  const size_t nDoF = np * 3;
  // set boundary condition
  aBCFlag.assign(nDoF, 0);
  for (unsigned int ip = 0; ip < np; ++ip) {
    const double px = aXY1[ip * 2 + 0];
    const double py = aXY1[ip * 2 + 1];
    if (fabs(fabs(py) - len) < 0.0001) {
      aBCFlag[ip * 3 + 0] = 1;
      aBCFlag[ip * 3 + 1] = 1;
      aVal[ip * 3 + 0] = 1;
    }
    if (fabs(px - len) < 0.0001) {
      aBCFlag[ip * 3 + 0] = 1;
      aBCFlag[ip * 3 + 1] = 1;
      aVal[ip * 3 + 0] = 1;
    }
  }
  for (int iip = loopIP_ind[1]; iip < loopIP_ind[2]; ++iip) {
    int ip0 = loopIP[iip];
    aBCFlag[ip0 * 3 + 0] = 1;
    aBCFlag[ip0 * 3 + 1] = 1;
  }
  //
  std::vector<unsigned int> psup_ind0, psup0;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind0, psup0,
      aTri1.data(), aTri1.size() / 3, 3,
      aXY1.size() / 2);
  //
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind0.data(), psup_ind0.size(), psup0.data(), psup0.size());
  ilu_A.SetPattern0(mat_A);
  //  ilu_A.Initialize_ILUk(mat_A, 5);
}

// -----------------------------------
// iproblem: 4
void SolveProblem_Stokes_Static(
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    std::vector<double> &aVal,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &aBCFlag) {
  const auto np = static_cast<unsigned int>(aXY1.size() / 2);
  const unsigned int nDoF = np * 3;
  // ---------------------
  double myu = 1.0;
  double g_x = 0.0;
  double g_y = -0.0;
  mat_A.setZero();
  std::vector<double> vec_b;
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_StokesStatic2D(
      mat_A, vec_b.data(),
      myu, g_x, g_y,
      aXY1.data(), aXY1.size() / 2,
      aTri1.data(), aTri1.size() / 3,
      aVal.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag, 0);
  // -----------------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.CopyValue(mat_A);
  ilu_A.Decompose();
  vec_x.resize(vec_b.size());
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    Solve_PCG(
        dfm2::ViewAsVectorXd(vec_b),
        dfm2::ViewAsVectorXd(vec_x),
        dfm2::ViewAsVectorXd(tmp0),
        dfm2::ViewAsVectorXd(tmp1),
        conv_ratio, iteration, mat_A, ilu_A);
  }
//  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  // ------------------------------
  dfm2::XPlusAY(aVal, nDoF, aBCFlag, 1.0, vec_x);
}

// ----------------------------------
// iprob:5
void SolveProblem_Stokes_Dynamic(
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    std::vector<double> &aVal,
    std::vector<double> &aVelo,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &aBCFlag) {
  const auto np = static_cast<unsigned int>(aXY1.size() / 2);
  const unsigned int nDoF = np * 3;
  // --------------------
  double myu = 1.0;
  double rho = 10;
  double g_x = 0.0;
  double g_y = -0.0;
  double dt_timestep = 0.01;
  double gamma_newmark = 0.6;
  mat_A.setZero();
  std::vector<double> vec_b;
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_StokesDynamic2D(
      mat_A, vec_b.data(),
      myu, rho, g_x, g_y,
      dt_timestep, gamma_newmark,
      aXY1.data(), aXY1.size() / 2,
      aTri1.data(), aTri1.size() / 3,
      aVal.data(), aVelo.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag, 0);
  // -------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.CopyValue(mat_A);
  ilu_A.Decompose();
  vec_x.resize(vec_b.size());
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    Solve_PCG(
        dfm2::ViewAsVectorXd(vec_b),
        dfm2::ViewAsVectorXd(vec_x),
        dfm2::ViewAsVectorXd(tmp0),
        dfm2::ViewAsVectorXd(tmp1),
        conv_ratio, iteration, mat_A, ilu_A);
  }
//  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  // --------------------
  dfm2::XPlusAYBZ(aVal, nDoF, aBCFlag,
                  dt_timestep * gamma_newmark, vec_x,
                  dt_timestep, aVelo);
  dfm2::XPlusAY(aVelo, nDoF, aBCFlag,
                1.0, vec_x);
}

// ------------------------
// iprob: 6
void SolveProblem_NavierStokes_Dynamic(
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    std::vector<double> &aVal,
    std::vector<double> &aVelo,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &aBCFlag) {
  const auto np = static_cast<unsigned int>(aXY1.size() / 2);
  const unsigned int nDoF = np * 3;
  // ----------------------
  double myu = 0.01;
  double rho = 1;
  double g_x = 0.0;
  double g_y = -0.0;
  double dt_timestep = 0.01;
  double gamma_newmark = 0.6;
  mat_A.setZero();
  std::vector<double> vec_b;
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_NavierStokes2D(
      mat_A, vec_b.data(),
      myu, rho, g_x, g_y,
      dt_timestep, gamma_newmark,
      aXY1.data(), aXY1.size() / 2,
      aTri1.data(), aTri1.size() / 3,
      aVal.data(), aVelo.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag, 0);
  // ----------------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.CopyValue(mat_A);
  ilu_A.Decompose();
  vec_x.resize(vec_b.size());
  Solve_PBiCGStab(vec_b.data(), vec_x.data(),
                  conv_ratio, iteration, mat_A, ilu_A);
//  SolveLinSys_BiCGStab(mat_A,vec_b,vec_x,ilu_A,
//                       conv_ratio, iteration);
  // -----------------------
  dfm2::XPlusAYBZ(aVal, nDoF, aBCFlag,
                  dt_timestep * gamma_newmark, vec_x,
                  dt_timestep, aVelo);
  dfm2::XPlusAY(aVelo, nDoF, aBCFlag,
                1.0, vec_x);
}

void DrawVelocityField(
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<double> &aVal) {
  std::vector<std::pair<double, delfem2::CColor> > colorMap;
  delfem2::ColorMap_BlueGrayRed(colorMap, -30, +30);
  delfem2::opengl::DrawMeshTri2D_ScalarP1(
      aXY1.data(), aXY1.size() / 2,
      aTri1.data(), aTri1.size() / 3,
      aVal.data() + 2, 3, colorMap);
  ::glColor3d(0, 0, 0);
  delfem2::opengl::DrawPoints2D_Vectors(
	  aXY1.data(), 
	  aXY1.size() / 2,
	  aVal.data(), 3, 0, 0.1);
  ::glPointSize(2);
  ::glColor3d(0, 0, 0);
  delfem2::opengl::DrawPoints2d_Points(aXY1);
}

void ProblemFluidCavity(
    dfm2::glfw::CViewer3 &viewer,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &loopIP_ind,
    const std::vector<int> &loopIP,
    double len) {
  std::vector<double> aVal;
  const size_t np = aXY1.size() / 2;
  dfm2::CMatrixSparse<double> mat_A;
  dfm2::CPreconditionerILU<double> ilu_A;
  std::vector<int> aBCFlag;
  glfwSetWindowTitle(viewer.window, "Stokes Static");
  aVal.assign(np * 3, 0.0);
  InitializeProblem_Fluid(
      aBCFlag, mat_A, ilu_A, aVal,
      aXY1, aTri1, loopIP_ind, loopIP, len);
  SolveProblem_Stokes_Static(
      mat_A, ilu_A, aVal,
      aXY1, aTri1, aBCFlag);
  for (unsigned int iframe = 0; iframe < 50; ++iframe) {
    viewer.DrawBegin_oldGL();
    DrawVelocityField(aXY1, aTri1, aVal);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  // ----------------------------
  glfwSetWindowTitle(viewer.window, "Stokes Dynamic");
  aVal.assign(np * 3, 0.0);
  std::vector<double> aVelo;
  aVelo.assign(np * 3, 0.0);
  InitializeProblem_Fluid(
      aBCFlag, mat_A, ilu_A, aVal,
      aXY1, aTri1, loopIP_ind, loopIP, len);
  for (unsigned int iframe = 0; iframe < 100; ++iframe) {
    SolveProblem_Stokes_Dynamic(
        mat_A, ilu_A, aVal, aVelo,
        aXY1, aTri1, aBCFlag);
    viewer.DrawBegin_oldGL();
    DrawVelocityField(aXY1, aTri1, aVal);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  // -----------------------------
  glfwSetWindowTitle(viewer.window, "Navier-Stokes");
  aVal.assign(np * 3, 0.0);
  aVelo.assign(np * 3, 0.0);
  InitializeProblem_Fluid(
      aBCFlag, mat_A, ilu_A, aVal,
      aXY1, aTri1, loopIP_ind, loopIP, len);
  for (unsigned int iframe = 0; iframe < 100; ++iframe) {
    SolveProblem_NavierStokes_Dynamic(
        mat_A, ilu_A, aVal, aVelo,
        aXY1, aTri1, aBCFlag);
    viewer.DrawBegin_oldGL();
    DrawVelocityField(aXY1, aTri1, aVal);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
}

void ProblemFluidTunnel(
    dfm2::glfw::CViewer3 &viewer,
    const std::vector<double> &aXY1,
    const std::vector<unsigned int> &aTri1,
    const std::vector<int> &loopIP_ind,
    const std::vector<int> &loopIP,
    double len) {
  const size_t np = aXY1.size() / 2;
  dfm2::CMatrixSparse<double> mat_A;
  dfm2::CPreconditionerILU<double> ilu_A;
  std::vector<int> aBCFlag;
  glfwSetWindowTitle(viewer.window, "Stokes Static");
  std::vector<double> aVal;
  aVal.assign(np * 3, 0.0);
  std::vector<double> aVelo;
  aVelo.assign(np * 3, 0.0);
  InitializeProblem_Fluid2(
      aBCFlag, mat_A, ilu_A, aVal,
      aXY1, aTri1, loopIP_ind, loopIP, len);
  SolveProblem_Stokes_Static(
      mat_A, ilu_A, aVal,
      aXY1, aTri1, aBCFlag);
  for (unsigned int iframe = 0; iframe < 100; ++iframe) {
    viewer.DrawBegin_oldGL();
    DrawVelocityField(aXY1, aTri1, aVal);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  // ----------------------------
  glfwSetWindowTitle(viewer.window, "Stokes Dynamic");
  aVal.assign(np * 3, 0.0);
  aVelo.assign(np * 3, 0.0);
  InitializeProblem_Fluid2(
      aBCFlag, mat_A, ilu_A, aVal,
      aXY1, aTri1, loopIP_ind, loopIP, len);
  for (unsigned int iframe = 0; iframe < 100; ++iframe) {
    SolveProblem_Stokes_Dynamic(
        mat_A, ilu_A, aVal, aVelo,
        aXY1, aTri1, aBCFlag);
    viewer.DrawBegin_oldGL();
    DrawVelocityField(aXY1, aTri1, aVal);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  // -----------------------------
  glfwSetWindowTitle(viewer.window, "Navier-Stokes");
  aVal.assign(np * 3, 0.0);
  aVelo.assign(np * 3, 0.0);
  InitializeProblem_Fluid2(
      aBCFlag, mat_A, ilu_A, aVal,
      aXY1, aTri1, loopIP_ind, loopIP, len);
  SolveProblem_NavierStokes_Dynamic(
      mat_A, ilu_A, aVal, aVelo,
      aXY1, aTri1, aBCFlag);
  for (unsigned int iframe = 0; iframe < 100; ++iframe) {
    SolveProblem_NavierStokes_Dynamic(
        mat_A, ilu_A, aVal, aVelo,
        aXY1, aTri1, aBCFlag);
    viewer.DrawBegin_oldGL();
    DrawVelocityField(aXY1, aTri1, aVal);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
}

int main() {
  dfm2::glfw::CViewer3 viewer(1.5);
  //
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  for (unsigned int itr = 0; itr < 10; ++itr) {
    std::vector<unsigned int> aTri1;
    std::vector<double> aXY1;
    std::vector<int> loopIP_ind, loopIP; // vtx on loop
    double len = 1.1;
    MakeMesh(
        aXY1, aTri1, loopIP_ind, loopIP,
        len);
    // ---------------------------
    ProblemScalar(viewer, aXY1, aTri1, len);
    ProblemSolid(viewer, aXY1, aTri1, len);
    ProblemFluidCavity(viewer, aXY1, aTri1, loopIP_ind, loopIP, len);
    ProblemFluidTunnel(viewer, aXY1, aTri1, loopIP_ind, loopIP, len);
  }
}
