/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cstdlib>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/ls_solver_block_sparse_ilu.h"
#include "delfem2/fem_solidhyper.h"
#include "delfem2/jagarray.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/msh_primitive.h"
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
    dfm2::LinearSystemSolver_BlockSparseILU &solver,
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aHex,
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
  solver.BeginMerge();
  for (unsigned int ih = 0; ih < aHex.size() / 8; ++ih) {
    double aP0[8][3];
    delfem2::FetchData<8, 3>(aP0, aHex.data() + ih * 8, aXYZ0.data());
    double aU[8][3];
    delfem2::FetchData<8, 3>(aU, aHex.data() + ih * 8, aDisp.data());
    //
    double W = 0.0, dW[8][3], ddW[8][8][3][3];
    std::fill_n(&dW[0][0], 8 * 3, 0.0);
    std::fill_n(&ddW[0][0][0][0], 8*8*3*3, 0.0);
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
      solver.vec_r[ip0 * 3 + 0] -= dW[ino][0];
      solver.vec_r[ip0 * 3 + 1] -= dW[ino][1];
      solver.vec_r[ip0 * 3 + 2] -= dW[ino][2];
    }
    delfem2::Merge<8, 8, 3, 3, double>(
        solver.matrix, aIP, aIP, ddW, solver.merge_buffer);
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    solver.matrix.val_dia_[ip * 9 + 0] += mass / (dt * dt);
    solver.matrix.val_dia_[ip * 9 + 4] += mass / (dt * dt);
    solver.matrix.val_dia_[ip * 9 + 8] += mass / (dt * dt);
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    solver.vec_r[ip * 3 + 0] += mass * gravity[0];
    solver.vec_r[ip * 3 + 1] += mass * gravity[1];
    solver.vec_r[ip * 3 + 2] += mass * gravity[2];
  }
  solver.Solve_PcgIlu();
  dfm2::XPlusAY(
      aDisp,
      nDoF, solver.dof_bcflag, 1.0, solver.vec_x);
  for (unsigned int ip = 0; ip < np; ++ip) {
    aVelo[ip * 3 + 0] += solver.vec_x[ip * 3 + 0] / dt;
    aVelo[ip * 3 + 1] += solver.vec_x[ip * 3 + 1] / dt;
    aVelo[ip * 3 + 2] += solver.vec_x[ip * 3 + 2] / dt;
  }
}

int main() {
  std::vector<double> vtx_xyz0;
  std::vector<unsigned int> hex_vtx;
  dfm2::MeshHex3_Grid(
      vtx_xyz0, hex_vtx,
      10, 10, 1, 0.2);
  dfm2::LinearSystemSolver_BlockSparseILU solver;
  InitializeMatrix(
      solver.matrix, solver.ilu_sparse,
      hex_vtx, vtx_xyz0);
  solver.dof_bcflag.resize(vtx_xyz0.size(), 0);
  for (unsigned int ip = 0; ip < vtx_xyz0.size() / 3; ++ip) {
    double x0 = vtx_xyz0[ip * 3 + 0];
    if (x0 > 1.0e-10) { continue; }
    solver.dof_bcflag[ip * 3 + 0] = 1;
    solver.dof_bcflag[ip * 3 + 1] = 1;
    solver.dof_bcflag[ip * 3 + 2] = 1;
  }
  std::vector<double> vtx_dispxyz(vtx_xyz0.size(), 0.0);
  std::vector<double> vtx_veloxyz(vtx_xyz0.size(), 0.0);
  const double dt = 0.01;
  const double mass = 0.5;
  const double stiff_c1 = 1000.0;
  const double stiff_c2 = 1000.0;
  const double stiff_comp = 10000.0;
  const double gravity[3] = {0, 0, -10};

  // ----------------------
  delfem2::glfw::CViewer3 viewer(1.5);
  viewer.view_rotation = std::make_unique<dfm2::ModelView_Ztop>();
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();

  delfem2::opengl::setSomeLighting();
  while (!glfwWindowShouldClose(viewer.window)) {
    Simulation(
        vtx_dispxyz, vtx_veloxyz, solver,
        vtx_xyz0, hex_vtx,
        dt, mass, stiff_c1, stiff_c2, stiff_comp, gravity);
    // -----
    viewer.DrawBegin_oldGL();
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0, 0, 0);
    delfem2::opengl::DrawMeshHex3D_EdgeDisp(
        vtx_xyz0.data(), vtx_xyz0.size() / 3,
        hex_vtx.data(), hex_vtx.size() / 8,
        vtx_dispxyz.data());
    //
    ::glEnable(GL_LIGHTING);
//    dfm2::opengl::DrawMeshHex3D_FaceNorm(aXYZ0.data(), aHex.data(), aHex.size() / 8);
    delfem2::opengl::DrawHex3D_FaceNormDisp(vtx_xyz0, hex_vtx, vtx_dispxyz);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
