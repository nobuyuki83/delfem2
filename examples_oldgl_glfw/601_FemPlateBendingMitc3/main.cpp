
#include <iostream>
#include <cmath>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/ls_ilu_block_sparse.h"
#include "delfem2/ls_block_sparse.h"
#include "delfem2/ls_solver_block_sparse.h"
#include "delfem2/view_vectorx.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/fem_mitc3.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/vec2.h"
#include "delfem2/jagarray.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// ------------------------

void InitializeProblem_PlateBendingMITC3(
    dfm2::LinearSystemSolver_BlockSparse &ls,
    dfm2::CPreconditionerILU<double> &ilu_sparse,
    const double lenx,
    [[maybe_unused]] const double leny,
    const std::vector<double> &vtx_xy_initial,
    const std::vector<unsigned int> &tri_vtx) {
  const std::size_t np = vtx_xy_initial.size() / 2;
  {  // initialize linear system
    std::vector<unsigned int> psup_ind, psup;
    dfm2::JArray_PSuP_MeshElem(
        psup_ind, psup,
        tri_vtx.data(), tri_vtx.size() / 3, 3,
        vtx_xy_initial.size() / 2);
    dfm2::JArray_Sort(psup_ind, psup);
    ls.Initialize(np, 3, psup_ind, psup);
  }
  ls.vec_x.resize(vtx_xy_initial.size() / 2 * 3, 0.0);
  // set boundary condition flag
  for (unsigned int ip = 0; ip < np; ++ip) {
    const double px = vtx_xy_initial[ip * 2 + 0];
    if (fabs(px - (-lenx * 0.5)) < 0.0001) {
      ls.dof_bcflag[ip * 3 + 0] = 1;
      ls.dof_bcflag[ip * 3 + 1] = 1;
      ls.dof_bcflag[ip * 3 + 2] = 1;
    }
  }
  // initialize sparse solver
  ilu_sparse.Initialize_ILUk(ls.matrix, 0);
}

/**
 *
 * @param vec_x
 * @param vec_r residual
 * @param mat_A
 * @param ilu_A
 * @param aBCFlag
 * @param thickness
 * @param myu
 * @param lambda
 * @param rho
 * @param gravity_z
 * @param aXY0
 * @param aTri
 */
void SolveProblem_PlateBendingMITC3(
    std::vector<double> &vec_x,
    std::vector<double> &vec_r,
    dfm2::CMatrixSparse<double> &mat_A,
    dfm2::CPreconditionerILU<double> &ilu_A,
    const std::vector<int> &aBCFlag,
    const double thickness,
    const double myu,
    const double lambda,
    const double rho,
    const double gravity_z,
    const std::vector<double> &aXY0,
    const std::vector<unsigned int> &aTri) {
  const std::size_t np = aXY0.size() / 2;
  const std::size_t nDoF = np * 3;
  //
  mat_A.setZero();
  vec_r.assign(nDoF, 0.0);
  assert(vec_x.size() == nDoF);
  dfm2::MergeLinSys_ShellStaticPlateBendingMITC3_MeshTri2D(
      mat_A, vec_r.data(),
      thickness, lambda, myu,
      rho, gravity_z,
      aXY0.data(), aXY0.size() / 2,
      aTri.data(), aTri.size() / 3,
      vec_x.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_r, aBCFlag, 0);
  //
  std::vector<double> vec_u;
  {
    ilu_A.CopyValue(mat_A);
    ilu_A.Decompose();
    vec_u.resize(vec_r.size());
    {
      const std::size_t n = vec_r.size();
      std::vector<double> tmp0(n), tmp1(n);
      auto vr = dfm2::ViewAsVectorXd(vec_r);
      auto vu = dfm2::ViewAsVectorXd(vec_u);
      auto vt = dfm2::ViewAsVectorXd(tmp0);
      auto vs = dfm2::ViewAsVectorXd(tmp1);
      std::vector<double> conv = dfm2::Solve_PCG(
          vr, vu, vt, vs,
          1.0e-5, 1000, mat_A, ilu_A);
      std::cout << "convergence   nitr:" << conv.size() << "    res:" << conv[conv.size() - 1] << std::endl;
    }
  }
  //
  dfm2::XPlusAY(vec_x, nDoF, aBCFlag,
                1.0, vec_u);
}

void MyDisplay(
    const std::vector<double> &vtx_xy_initial,
    const std::vector<unsigned int> &tri_vtx,
    const std::vector<double> &vtx_value) {
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0, 0, 0);
  delfem2::opengl::DrawMeshTri2D_Edge(tri_vtx, vtx_xy_initial);
  {
    assert(vtx_value.size() / 3 == vtx_xy_initial.size() / 2);
    ::glColor3d(1, 0, 0);
    ::glBegin(GL_LINES);
    for (size_t itri = 0; itri < tri_vtx.size() / 3; itri++) {
      const unsigned int i0 = tri_vtx[itri * 3 + 0];
      const unsigned int i1 = tri_vtx[itri * 3 + 1];
      const unsigned int i2 = tri_vtx[itri * 3 + 2];
      const double p0[3] = {vtx_xy_initial[i0 * 2 + 0], vtx_xy_initial[i0 * 2 + 1], vtx_value[i0 * 3 + 0]};
      const double p1[3] = {vtx_xy_initial[i1 * 2 + 0], vtx_xy_initial[i1 * 2 + 1], vtx_value[i1 * 3 + 0]};
      const double p2[3] = {vtx_xy_initial[i2 * 2 + 0], vtx_xy_initial[i2 * 2 + 1], vtx_value[i2 * 3 + 0]};
      ::glVertex3dv(p0);
      ::glVertex3dv(p1);
      ::glVertex3dv(p1);
      ::glVertex3dv(p2);
      ::glVertex3dv(p2);
      ::glVertex3dv(p0);
    }
    ::glEnd();
  }
}

int main() {
  const double lenx = 1.0;
  const double leny = 0.2;
  std::vector<unsigned int> tri_vtx;
  std::vector<double> vtx_xy_initial;
  {
    std::vector<std::vector<double> > aaXY;
    {
      aaXY.resize(1);
      aaXY[0].push_back(-lenx * 0.5);
      aaXY[0].push_back(-leny * 0.5);
      aaXY[0].push_back(+lenx * 0.5);
      aaXY[0].push_back(-leny * 0.5);
      aaXY[0].push_back(+lenx * 0.5);
      aaXY[0].push_back(+leny * 0.5);
      aaXY[0].push_back(-lenx * 0.5);
      aaXY[0].push_back(+leny * 0.5);
    }
    std::vector<dfm2::CDynPntSur> aPo2D;
    std::vector<dfm2::CDynTri> aETri;
    std::vector<dfm2::CVec2d> aVec2;
    GenMesh(aPo2D, aETri, aVec2,
            aaXY, 0.03, 0.03);
    MeshTri2D_Export(vtx_xy_initial, tri_vtx,
                     aVec2, aETri);
    std::cout << "  ntri;" << tri_vtx.size() / 3 << "  nXY:" << vtx_xy_initial.size() / 2 << std::endl;
  }
  // -------------
  dfm2::LinearSystemSolver_BlockSparse ls;
  dfm2::CPreconditionerILU<double> ilu_A;
  InitializeProblem_PlateBendingMITC3(
      ls, ilu_A,
      lenx, leny, vtx_xy_initial, tri_vtx);
  // -----------
  const double thickness = 0.05;
  const double myu = 10000.0;
  const double lambda = 0.0;
  const double rho = 1.0;
  const double gravity_z = -10.0;
  SolveProblem_PlateBendingMITC3(
      ls.vec_x, ls.vec_r, ls.matrix, ilu_A, ls.dof_bcflag,
      thickness, myu, lambda, rho, gravity_z,
      vtx_xy_initial, tri_vtx);

  delfem2::glfw::CViewer3 viewer(0.8);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  viewer.view_rotation = std::make_unique<dfm2::ModelView_Ztop>();
  delfem2::opengl::setSomeLighting();

  {
    assert(fabs(lambda) < 1.0e-10);
    const double E = myu * 2.0;
    const double I = thickness * thickness * thickness * leny / 12.0;
    const double W = thickness * lenx * leny * rho * gravity_z;
    const double w = W / lenx;
    const double disp = w * (lenx * lenx * lenx * lenx) / (8.0 * E * I);
    std::cout << "disp:" << disp << std::endl;
    for (size_t ip = 0; ip < vtx_xy_initial.size() / 2; ++ip) {
      const double px = vtx_xy_initial[ip * 2 + 0];
      if (fabs(px - (+lenx * 0.5)) > 0.0001) { continue; }
      std::cout << ls.vec_x[ip * 3 + 0] << std::endl;
    }
  }

  while (!glfwWindowShouldClose(viewer.window)) {
    viewer.DrawBegin_oldGL();
    MyDisplay(vtx_xy_initial, tri_vtx, ls.vec_x);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


