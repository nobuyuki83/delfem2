/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

/*
#include "delfem2/ls_ilu_block_sparse.h"
#include "delfem2/ls_block_sparse.h"
#include "delfem2/vecxitrsol.h"
 */
#include "delfem2/ls_solver_block_sparse_ilu.h"
#include "delfem2/svd3.h"
#include "delfem2/mat3.h"
#include "delfem2/mat3_funcs.h"
#include "delfem2/fem_solidlinear.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/mshmisc.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/dtri_topology.h"
#include "delfem2/jagarray.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// ----------------------------------

void GenMesh(
    std::vector<dfm2::CVec2d> &aVec2,
    std::vector<dfm2::CDynPntSur> &aPo2D,
    std::vector<dfm2::CDynTri> &aETri,
    double elen,
    const std::vector<std::vector<double> > &aaXY) {
  std::vector<int> loopIP_ind, loopIP;
  {
    dfm2::JArray_FromVecVec_XY(
        loopIP_ind, loopIP, aVec2,
        aaXY);
    if (!dfm2::CheckInputBoundaryForTriangulation(loopIP_ind, aVec2)) {
      return;
    }
    dfm2::FixLoopOrientation(
        loopIP,
        loopIP_ind, aVec2);
    if (elen > 10e-10) {
      dfm2::ResamplingLoop(
          loopIP_ind, loopIP, aVec2,
          elen);
    }
  }
  // ------
  Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                 loopIP_ind, loopIP);
  if (elen > 1.0e-10) {
    dfm2::CInputTriangulation_Uniform param(1.0);
    std::vector<unsigned int> aFlgPnt(aVec2.size());
    std::vector<unsigned int> aFlgTri(aETri.size(), 0);
    MeshingInside(
        aPo2D, aETri, aVec2, aFlgPnt, aFlgTri,
        aVec2.size(), 0, elen, param);
  }
}

void RotationAtMeshPoints(
    std::vector<double> &aR,
    const std::vector<double> &aXYZ,
    const std::vector<double> &aDisp,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup) {
  const size_t np = aXYZ.size() / 3;
  aR.resize(np * 9);
  for (unsigned int ip = 0; ip < np; ++ip) {
    dfm2::CVec3d Pi(aXYZ[ip * 3 + 0], aXYZ[ip * 3 + 1], aXYZ[ip * 3 + 2]);
    dfm2::CVec3d pi(aXYZ[ip * 3 + 0] + aDisp[ip * 3 + 0],
                    aXYZ[ip * 3 + 1] + aDisp[ip * 3 + 1],
                    aXYZ[ip * 3 + 2] + aDisp[ip * 3 + 2]);
    dfm2::CMat3d A;
    A.setZero();
    for (unsigned int jjp = psup_ind[ip]; jjp < psup_ind[ip + 1]; ++jjp) {
      const unsigned int jp = psup[jjp];
      dfm2::CVec3d Pj(aXYZ[jp * 3 + 0], aXYZ[jp * 3 + 1], aXYZ[jp * 3 + 2]);
      dfm2::CVec3d pj(aXYZ[jp * 3 + 0] + aDisp[jp * 3 + 0],
                      aXYZ[jp * 3 + 1] + aDisp[jp * 3 + 1],
                      aXYZ[jp * 3 + 2] + aDisp[jp * 3 + 2]);
      A += dfm2::Mat3_OuterProduct(pj - pi, Pj - Pi);
    }
    dfm2::GetRotPolarDecomp(aR.data() + ip * 9,
                            A.p_, 100);
  }
}


// ---------------------------

// --------------------------------------------


void InitializeProblem_ShellEigenPB(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTet,
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A) {
  const size_t np = aXYZ.size() / 3;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind, psup,
      aTet.data(), aTet.size() / 4, 4,
      (int) aXYZ.size() / 3);
  dfm2::JArray_Sort(psup_ind, psup);
  mat_A.Initialize(static_cast<unsigned int>(np), 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(), psup.size());
  ilu_A.SetPattern0(mat_A);
}

// ------------------------------------------------------

void Solve_Linear(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTet,
    std::vector<double> &aDisp,
    std::vector<double> &aVelo,
    dfm2::CMatrixSparse<double> &mat_A,
    std::vector<double> &vec_b,
    dfm2::CPreconditionerILU<double> &ilu_A,
    double myu,
    double lambda,
    double rho,
    const double gravity[3],
    double dt,
    const std::vector<int> &aBCFlag) {
  mat_A.setZero();
  vec_b.assign(aXYZ.size(), 0.0);
  dfm2::MergeLinSys_SolidLinear_BEuler_MeshTet3D(
      mat_A, vec_b.data(),
      myu, lambda,
      rho, gravity,
      dt,
      aXYZ.data(), static_cast<unsigned int>(aXYZ.size() / 3),
      aTet.data(), static_cast<unsigned int>(aTet.size() / 4),
      aDisp.data(),
      aVelo.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag, 0);
  //
  ilu_A.CopyValue(mat_A);
  ilu_A.Decompose();
  const auto nDoF = static_cast<unsigned int>(aXYZ.size());
  std::vector<double> dv(nDoF, 0.0);
  std::vector<double> aConv = Solve_PBiCGStab(
      vec_b.data(), dv.data(),
      1.0e-4, 1000, mat_A, ilu_A);
  //
  dfm2::XPlusAYBZ(aDisp, nDoF, aBCFlag,
                  dt, dv,
                  dt, aVelo);
  dfm2::XPlusAY(aVelo, nDoF, aBCFlag,
                1.0, dv);
  std::cout << "conv; " << aConv.size() << std::endl;
}

void Solve_StiffnessWarping(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTet,
    std::vector<double> &aDisp,
    std::vector<double> &aVelo,
    std::vector<double> &aR,
    dfm2::CMatrixSparse<double> &mat_A,
    std::vector<double> &vec_b,
    dfm2::CPreconditionerILU<double> &ilu_A,
    double myu,
    double lambda,
    double rho,
    const double gravity[3],
    double dt,
    const std::vector<int> &aBCFlag,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup) {
  RotationAtMeshPoints(aR,
                       aXYZ, aDisp, psup_ind, psup);
  // ----------------------
  mat_A.setZero();
  vec_b.assign(aXYZ.size(), 0.0);
  dfm2::MergeLinSys_SolidStiffwarp_BEuler_MeshTet3D(
      mat_A, vec_b.data(),
      myu, lambda,
      rho, gravity,
      dt,
      aXYZ.data(), static_cast<unsigned int>(aXYZ.size() / 3),
      aTet.data(), static_cast<unsigned int>(aTet.size() / 4),
      aDisp.data(),
      aVelo.data(),
      aR);
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag, 0);
  // -------------------
  ilu_A.CopyValue(mat_A);
  ilu_A.Decompose();
  const auto nDoF = static_cast<unsigned int>(aXYZ.size());
  std::vector<double> dv(nDoF, 0.0);
  std::vector<double> aConv = Solve_PBiCGStab(
	  vec_b.data(), dv.data(),
	  1.0e-4, 1000, mat_A, ilu_A);
  dfm2::XPlusAYBZ(aDisp, nDoF, aBCFlag,
                  dt, dv,
                  dt, aVelo);
  dfm2::XPlusAY(aVelo, nDoF, aBCFlag,
                1.0, dv);
  std::cout << "conv; " << aConv.size() << std::endl;
}

// --------------------------------------------------------------

void myGlutDisplay(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTet,
    std::vector<double> &aDisp,
    std::vector<double> &aR ) {
  {
    float color[4] = {200.0 / 256.0, 200.0 / 256.0, 200.0 / 256.0, 1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
//    glShadeModel(GL_SMOOTH);
    glShadeModel(GL_FLAT);
  }
  ::glDisable(GL_LIGHTING);

  { // defomred edge
    ::glColor3d(0, 0, 0);
    delfem2::opengl::DrawMeshTet3D_EdgeDisp(
		aXYZ.data(),
		aTet.data(), static_cast<unsigned int>(aTet.size() / 4),
		aDisp.data(),
		1.0);
  }

  ::glDisable(GL_LIGHTING);
  for (std::size_t ip = 0; ip < aXYZ.size() / 3; ++ip) {
    dfm2::CVec3d pi(aXYZ[ip * 3 + 0] + aDisp[ip * 3 + 0],
                    aXYZ[ip * 3 + 1] + aDisp[ip * 3 + 1],
                    aXYZ[ip * 3 + 2] + aDisp[ip * 3 + 2]);
    dfm2::CVec3d ex(aR[ip * 9 + 0], aR[ip * 9 + 3], aR[ip * 9 + 6]);
    dfm2::CVec3d ey(aR[ip * 9 + 1], aR[ip * 9 + 4], aR[ip * 9 + 7]);
    dfm2::CVec3d ez(aR[ip * 9 + 2], aR[ip * 9 + 5], aR[ip * 9 + 8]);
    ::glBegin(GL_LINES);
    ::glColor3d(1, 0, 0);
    delfem2::opengl::myGlVertex(pi);
    delfem2::opengl::myGlVertex(pi + 0.04 * ex);
    ::glColor3d(0, 1, 0);
    delfem2::opengl::myGlVertex(pi);
    delfem2::opengl::myGlVertex(pi + 0.04 * ey);
    ::glColor3d(0, 0, 1);
    delfem2::opengl::myGlVertex(pi);
    delfem2::opengl::myGlVertex(pi + 0.04 * ez);
    ::glEnd();
  }

  {
    ::glEnable(GL_LIGHTING);
    {
      float color[4] = {180.0 / 256.0, 180.0 / 256.0, 130.0 / 256.0, 1.0f};
      ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
      ::glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
      glShadeModel(GL_FLAT);
    }
    delfem2::opengl::DrawMeshTet3D_FaceNorm(
        aXYZ.data(),
        aTet.data(), static_cast<unsigned int>(aTet.size() / 4));

  }

}

int main() {
  std::vector<unsigned int> tet_vtx;
  std::vector<double> vtx_xyz;
  std::vector<double> vtx_dispxyz;
  std::vector<double> vtx_veloxyz;
  std::vector<int> vtx_bcflag;
  double dt = 0.03;
  double myu = 200.0;
  double lambda = 1.0;
  double rho = 1.0;
  const double gravity[3] = {0.0, -4.0, 0.0};

  dfm2::LinearSystemSolver_BlockSparseILU solver;
  std::vector<unsigned int> psup_ind, psup;
  std::vector<double> aR;
  {
    std::vector<std::vector<double> > aaXY;
    {
      const double aXY[8] = {
          -1, -0.1,
          +1, -0.1,
          +1, +0.1,
          -1, +0.1};
      aaXY.emplace_back(aXY, aXY + 8);
    }
    std::vector<dfm2::CVec2d> aVec2;
    std::vector<dfm2::CDynPntSur> aPo2D;
    std::vector<dfm2::CDynTri> aETri;
    GenMesh(aVec2, aPo2D, aETri,
            0.075, aaXY);
    std::vector<double> aXY;
    std::vector<unsigned int> aTri;
    CMeshTri2D(aXY, aTri,
               aVec2, aETri);
    dfm2::ExtrudeTri2Tet(3, 0.075,
                         vtx_xyz, tet_vtx,
                         aXY, aTri);
  }
  vtx_dispxyz.assign(vtx_xyz.size(), 0.0);
  vtx_veloxyz.assign(vtx_xyz.size(), 0.0);
  vtx_bcflag.assign(vtx_xyz.size(), 0);
  for (std::size_t ip = 0; ip < vtx_xyz.size() / 3; ++ip) {
    double x0 = vtx_xyz[ip * 3 + 0];
    if (fabs(x0 + 1) < 1.0e-10) {
      vtx_bcflag[ip * 3 + 0] = 1;
      vtx_bcflag[ip * 3 + 1] = 1;
      vtx_bcflag[ip * 3 + 2] = 1;
    }
  }
  InitializeProblem_ShellEigenPB(
      vtx_xyz, tet_vtx, psup_ind, psup, solver.matrix, solver.ilu_sparse);
  RotationAtMeshPoints(
      aR,
      vtx_xyz, vtx_dispxyz, psup_ind, psup);

  delfem2::glfw::CViewer3 viewer(2.0);
  {
    viewer.view_rotation = std::make_unique<delfem2::ModelView_Ytop>();
  }
  bool is_stiffness_warping = true;
  viewer.keypress_callbacks.emplace_back([&is_stiffness_warping](int key, int){
    if( key == GLFW_KEY_SPACE ) {
      is_stiffness_warping = !is_stiffness_warping;
    }
  });
  //
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  //
  while (!glfwWindowShouldClose(viewer.window)) {
    if (is_stiffness_warping) {
      Solve_StiffnessWarping(
          vtx_xyz, tet_vtx, vtx_dispxyz, vtx_veloxyz, aR,
          solver.matrix, solver.vec_r, solver.ilu_sparse,
          myu, lambda, rho, gravity, dt,
          vtx_bcflag, psup_ind, psup);
    }
    else {
      Solve_Linear(
          vtx_xyz, tet_vtx, vtx_dispxyz, vtx_veloxyz,
          solver.matrix, solver.vec_r, solver.ilu_sparse,
          myu, lambda, rho, gravity, dt,
          vtx_bcflag);
    }
    // -----
    viewer.DrawBegin_oldGL();
    myGlutDisplay(
        vtx_xyz, tet_vtx, vtx_dispxyz, aR );
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);

}
