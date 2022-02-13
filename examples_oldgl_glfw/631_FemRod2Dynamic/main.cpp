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
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/vec2.h"
#include "delfem2/femutil.h"
#include "delfem2/fem_rod2.h"
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
    std::vector<double> &vtx_xy_deformed,
    const std::vector<double> &vtx_xy_initial){
  ::glLineWidth(1);
  ::glColor3d(0, 0, 0);
  DrawPolyline(vtx_xy_initial);
  ::glPointSize(5);
  delfem2::opengl::DrawPoints2d_Points(vtx_xy_initial);
  //
  ::glColor3d(1, 0, 0);
  DrawPolyline(vtx_xy_deformed);
  delfem2::opengl::DrawPoints2d_Points(vtx_xy_deformed);
}


template<class LIN_SYS_SOLVER>
void StepTime(
    std::vector<double> &vtx_xy_deformd,
    std::vector<double> &vtx_uv,
    const std::vector<double> &vtx_xy_initial,
    double dt,
    double stiff_stretch,
    double stiff_bend,
    double mass_point,
    const double gravity[2],
    LIN_SYS_SOLVER &sparse) {
  auto np = static_cast<unsigned int>(vtx_xy_initial.size() / 2);
  assert(np >= 3);
  assert(sparse.nblk() == np && sparse.ndim() == 2);
  for (unsigned int ip = 0; ip < np; ++ip) {
    vtx_xy_deformd[ip * 2 + 0] += dt * vtx_uv[ip * 2 + 0];
    vtx_xy_deformd[ip * 2 + 1] += dt * vtx_uv[ip * 2 + 1];
  }
  sparse.BeginMerge();
  double W = 0.0;
  for (unsigned int ihinge = 0; ihinge < np - 2; ++ihinge) {
    const unsigned int aIP[3] = {ihinge, ihinge + 1, ihinge + 2};
    double aP[3][2];
    dfm2::FetchData<3, 2>(aP, aIP, vtx_xy_initial.data());
    double ap[3][2];
    dfm2::FetchData<3, 2>(ap, aIP, vtx_xy_deformd.data());
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
  for (unsigned int ip = 0; ip < np; ++ip) {
    sparse.AddValueToDiagonal(ip, 0, mass_point / (dt * dt));
    sparse.AddValueToDiagonal(ip, 1, mass_point / (dt * dt));
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    sparse.vec_r[ip * 2 + 0] += mass_point * gravity[0];
    sparse.vec_r[ip * 2 + 1] += mass_point * gravity[1];
  }
  std::cout << W << std::endl;
  sparse.Solve();
  for (unsigned int ip = 0; ip < np; ++ip) {
    vtx_xy_deformd[ip * 2 + 0] += sparse.vec_x[ip * 2 + 0];
    vtx_xy_deformd[ip * 2 + 1] += sparse.vec_x[ip * 2 + 1];
    vtx_uv[ip * 2 + 0] += sparse.vec_x[ip * 2 + 0] / dt;
    vtx_uv[ip * 2 + 1] += sparse.vec_x[ip * 2 + 1] / dt;
  }
}

int main() {
  const double dt = 1.0 / 60.0;
  const double stiff_stretch = 10.0;
  const double stiff_bend = 0.001;
  const double mass_point = 1.0e-5;
  const double gravity[2] = {0, -10};
  std::vector<double> vtx_xy_deformed;  // vertices of polyline in the initial config
  {
    const unsigned int num_points = 11;
    for (unsigned int i = 0; i < num_points; ++i) {
      vtx_xy_deformed.push_back(double(i) / num_points);
      vtx_xy_deformed.push_back(0.0);
    }
  }
  const std::vector<double> vtx_xy_initial = vtx_xy_deformed;
  std::mt19937 rndeng(std::random_device{}());
  std::uniform_real_distribution<double> distm1p1(-1, +1);

  // opengl set up
  dfm2::glfw::InitGLOld();
  dfm2::glfw::CViewer2 viewer;
  {
    viewer.view_height = 0.7f;
    viewer.trans[0] = -0.2f;
    viewer.trans[1] = +0.3f;
  }
  viewer.OpenWindow();  // opengl start here

  while (true) {
    {
      glfwSetWindowTitle(viewer.window, "penta-diagonal solver free rotation");
      delfem2::LinearSystemSolver_BlockPentaDiagonal<2> ls_solver;
      ls_solver.Initialize(vtx_xy_deformed.size() / 2);
      std::cout << "initial" << std::endl;
      ls_solver.dof_bcflag[0] = 1;
      ls_solver.dof_bcflag[1] = 1;
      vtx_xy_deformed = vtx_xy_initial;
      for (unsigned int i = 0; i < vtx_xy_deformed.size(); ++i) {
        if (ls_solver.dof_bcflag[i] != 0) { continue; }
        vtx_xy_deformed[i] += 0.1*distm1p1(rndeng);
      }
      std::vector<double> vtx_uv(vtx_xy_deformed.size(), 0.);
      for(unsigned int iframe=0;iframe<100;++iframe){
        StepTime(
            vtx_xy_deformed, vtx_uv,
            vtx_xy_initial,
            dt, stiff_stretch, stiff_bend, mass_point, gravity,
            ls_solver);
        viewer.DrawBegin_oldGL();
        Draw(vtx_xy_deformed, vtx_xy_initial);
        viewer.SwapBuffers();
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) { break; }
      }
      if (glfwWindowShouldClose(viewer.window)) { break; }
    }
    {
      glfwSetWindowTitle(viewer.window, "penta-diagonal solver fix rotation");
      delfem2::LinearSystemSolver_BlockPentaDiagonal<2> ls_solver;
      ls_solver.Initialize(vtx_xy_deformed.size() / 2);
      ls_solver.dof_bcflag[0] = 1;
      ls_solver.dof_bcflag[1] = 1;
      ls_solver.dof_bcflag[2] = 1;
      ls_solver.dof_bcflag[3] = 1;
      vtx_xy_deformed = vtx_xy_initial;
      for (unsigned int i = 0; i < vtx_xy_deformed.size(); ++i) {
        if (ls_solver.dof_bcflag[i] != 0) { continue; }
        vtx_xy_deformed[i] += 0.1*distm1p1(rndeng);
      }
      std::vector<double> vtx_uv(vtx_xy_deformed.size(), 0.);
      for(unsigned int iframe=0;iframe<100;++iframe){
        StepTime(
            vtx_xy_deformed, vtx_uv,
            vtx_xy_initial,
            dt, stiff_stretch, stiff_bend, mass_point, gravity,
            ls_solver);
        viewer.DrawBegin_oldGL();
        Draw(vtx_xy_deformed, vtx_xy_initial);
        viewer.SwapBuffers();
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) { break; }
      }
      if (glfwWindowShouldClose(viewer.window)) { break; }
    }
    {
      glfwSetWindowTitle(viewer.window, "block sparse solver free rotation");
      delfem2::LinearSystemSolver_BlockSparse ls_solver;
      {
        std::vector<unsigned int> psup_ind, psup;
        delfem2::JArray_PSuP_Hair(
            psup_ind, psup,
            {0, static_cast<unsigned int>(vtx_xy_deformed.size() / 2)});
        ls_solver.Initialize(
            vtx_xy_deformed.size() / 2, 2,
            psup_ind, psup);
      }
      ls_solver.dof_bcflag[0] = 1;
      ls_solver.dof_bcflag[1] = 1;
      vtx_xy_deformed = vtx_xy_initial;
      for (unsigned int i = 0; i < vtx_xy_deformed.size(); ++i) {
        if (ls_solver.dof_bcflag[i] != 0) { continue; }
        vtx_xy_deformed[i] += 0.1*distm1p1(rndeng);
      }
      std::vector<double> vtx_uv(vtx_xy_deformed.size(), 0.);
      for(unsigned int iframe=0;iframe<100;++iframe){
        StepTime(
            vtx_xy_deformed, vtx_uv,
            vtx_xy_initial,
            dt, stiff_stretch, stiff_bend, mass_point, gravity,
            ls_solver);
        viewer.DrawBegin_oldGL();
        Draw(vtx_xy_deformed, vtx_xy_initial);
        viewer.SwapBuffers();
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) { break; }
      }
      if (glfwWindowShouldClose(viewer.window)) { break; }
    }
  }
  viewer.ExitIfClosed();

  return 0;
}


