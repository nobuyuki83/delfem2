/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <cstdlib>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/fem_invertiblefem.h"
#include "delfem2/sampling.h"
#include "delfem2/ls_solver_block_sparse.h"
#include "delfem2/mshuni.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// ---------------------------------------

int main() {
  std::mt19937 rndeng(0);
  std::uniform_real_distribution<double> dist_01(0, 1);

  auto neohook = [](
      double dW[3], double ddW[3][3],
      double l0, double l1, double l2) -> double {
    dW[0] = (l0 - 1);
    dW[1] = (l1 - 1);
    dW[2] = (l2 - 1);
    ddW[0][0] = 1;
    ddW[0][1] = 0;
    ddW[0][2] = 0;
    ddW[1][0] = 0;
    ddW[1][1] = 1;
    ddW[1][2] = 0;
    ddW[2][0] = 0;
    ddW[2][1] = 0;
    ddW[2][2] = 1;
    return (l0 - 1) * (l0 - 1) * 0.5
        + (l1 - 1) * (l1 - 1) * 0.5
        + (l2 - 1) * (l2 - 1) * 0.5;
  };

  std::vector<double> Pos0;
  Pos0 = {
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1};
  std::vector<unsigned int> aTet = {0, 1, 2, 3};
  std::vector<double> pos0(Pos0.size());
  for(auto& v: pos0){
    v = dist_01(rndeng);
  }

  dfm2::LinearSystemSolver_BlockSparse sparse;
  {
    std::vector<unsigned int> psup_ind, psup;
    dfm2::JArray_PSuP_MeshElem(
        psup_ind, psup,
        aTet.data(), aTet.size() / 4, 4, 4);
    sparse.Initialize(4, 3, psup_ind, psup);
  }

  delfem2::glfw::CViewer3 viewer(1.5);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();

  double time_prev = glfwGetTime();
  while (!glfwWindowShouldClose(viewer.window)) {
    if( glfwGetTime() - time_prev > 1 ){
      for(auto& v: pos0){
        v = dist_01(rndeng);
      }
      time_prev = glfwGetTime();
    }
    {
      sparse.BeginMerge();
      for(unsigned int it=0;it<aTet.size()/4;++it) {
        unsigned int i0 = aTet[it*4+0];
        unsigned int i1 = aTet[it*4+1];
        unsigned int i2 = aTet[it*4+2];
        unsigned int i3 = aTet[it*4+3];
        const double P0[4][3] = {
            {Pos0[i0*3+0],Pos0[i0*3+1],Pos0[i0*3+2]},
            {Pos0[i1*3+0],Pos0[i1*3+1],Pos0[i1*3+2]},
            {Pos0[i2*3+0],Pos0[i2*3+1],Pos0[i2*3+2]},
            {Pos0[i3*3+0],Pos0[i3*3+1],Pos0[i3*3+2]} };
        const double p0[4][3] = {
            {pos0[i0*3+0],pos0[i0*3+1],pos0[i0*3+2]},
            {pos0[i1*3+0],pos0[i1*3+1],pos0[i1*3+2]},
            {pos0[i2*3+0],pos0[i2*3+1],pos0[i2*3+2]},
            {pos0[i3*3+0],pos0[i3*3+1],pos0[i3*3+2]} };
        double dW0[4][3], ddW0[4][4][3][3];
        const double W0 = delfem2::WdWddW_InvertibleFEM(
            dW0, ddW0,
            P0, p0, neohook);
        std::cout << W0 << std::endl;
        {
          const unsigned int *aIP = aTet.data() + it * 4;
          sparse.Merge<4, 4, 3, 3>(aIP, aIP, ddW0);
        }
        for (unsigned int iblk = 0; iblk < sparse.nblk(); ++iblk) {
          for (unsigned int idim = 0; idim < 3; ++idim) {
            sparse.vec_r[iblk * 3 + idim] -= dW0[iblk][idim];
          }
        }
      }
      for (unsigned int iblk = 0; iblk < sparse.nblk(); ++iblk) {
        for (unsigned int idim = 0; idim < 3; ++idim) {
          sparse.AddValueToDiagonal(iblk, idim, 1.0);
        }
      }
      sparse.Solve();
      std::cout << sparse.conv_hist.size() << std::endl;
      unsigned int ndof = sparse.ndof();
      for (unsigned int i = 0; i < ndof; ++i) {
        pos0[i] += sparse.vec_x[i];
      }
    }
    //
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    delfem2::opengl::DrawMeshTet3D_FaceNorm(pos0.data(), aTet.data(), aTet.size()/4);
    ::glColor3b(0,0,0);
    ::glDisable(GL_LIGHTING);
    delfem2::opengl::DrawMeshTet3D_Edge(pos0.data(), pos0.size()/3, aTet.data(), aTet.size()/4);
    //
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
