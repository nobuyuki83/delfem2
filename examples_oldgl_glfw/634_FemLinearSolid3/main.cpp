/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <random>
#include <vector>
#include <cstdlib>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/mshprimitive.h"
#include "delfem2/lsilu_mats.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/lsmats.h"
#include "delfem2/lsvecx.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/femsolidlinear.h"
#include "delfem2/mshuni.h"
#include "delfem2/femutil.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"



namespace dfm2 = delfem2;

// --------------------------------------------------------------

void InitializeMatrix(
    dfm2::CMatrixSparse<double>& smat,
    delfem2::CPreconditionerILU<double>& ilu,
    const std::vector<unsigned int>& aHex,
    const std::vector<double>& aXYZ)
{
  const unsigned int np = aXYZ.size()/3;
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind, psup,
      aHex.data(), aHex.size()/8, 8,
      (int)aXYZ.size()/3);
  smat.Initialize(np, 3, true);
  smat.SetPattern(
      psup_ind.data(), psup_ind.size(),
      psup.data(), psup.size());
  //ilu.Initialize_ILU0(smat);
  ilu.Initialize_ILUk(smat,2);
}

void Simulation(
    std::vector<double>& aDisp,
    dfm2::CMatrixSparse<double>& smat,
    delfem2::CPreconditionerILU<double>& silu,
    //
    const std::vector<double>& aXYZ0,
    const std::vector<unsigned int>& aHex,
    const std::vector<int>& aBCFlag,
    //
    double mass,
    double myu,
    double lambda,
    const double gravity[3])
{
  const unsigned int np = aXYZ0.size()/3;
  const unsigned int nDoF = np*3;
  smat.setZero();
  {
    double ddW[8][8][3][3];
    {
      const double gravity_zero[3] = {0, 0, 0};
      double aP0[8][3], aU[8][3];
      delfem2::FetchData<8, 3>(aP0, aHex.data(), aXYZ0.data());
      delfem2::FetchData<8, 3>(aU, aHex.data(), aDisp.data());
      double dW[8][3];
      std::fill_n(&dW[0][0], 8 * 3, 0.0);
      std::fill_n(&ddW[0][0][0][0], 8 * 8 * 3 * 3, 0.0);
      delfem2::elemMatRes_LinearSolidGravity3_Static_Q1(
          myu, lambda,
          0, gravity_zero,
          aP0, aU, ddW, dW);
    }
    std::vector<unsigned int> tmp_buffer;
    for (unsigned int ih = 0; ih < aHex.size() / 8; ++ih) {
      const unsigned int *aIP = aHex.data() + ih * 8;
      delfem2::Merge<8, 8, 3, 3, double>(smat, aIP, aIP, ddW, tmp_buffer);
    }
  }
  std::vector<double> vecb(nDoF, 0.0);
  { // comput rhs vectors
    auto vb = dfm2::CVecXd(vecb);
    auto vd = dfm2::CVecXd(aDisp);
    AddMatVec(vb, 0.0, -1.0, smat, vd);
    std::cout << "energy" << vb.dot(vd) << std::endl;
  }
  for(unsigned int ip=0;ip<np;++ip) {
    vecb[ip*3+0] += mass*gravity[0];
    vecb[ip*3+1] += mass*gravity[1];
    vecb[ip*3+2] += mass*gravity[2];
  }
  smat.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vecb, aBCFlag,0);
  // --------------------------------
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  silu.CopyValue(smat);
  silu.Decompose();
  std::vector<double> vecx(vecb.size());
  {
    const std::size_t n = vecb.size();
    std::vector<double> tmp0(n), tmp1(n);
    std::vector<double> aHist = Solve_PCG(
        dfm2::CVecXd(vecb),dfm2::CVecXd(vecx),dfm2::CVecXd(tmp0),dfm2::CVecXd(tmp1),
        conv_ratio, iteration,
        smat, silu);
    std::cout << "nconv:" << aHist.size() << std::endl;
  }
  // ------------------------------
  dfm2::XPlusAY(
      aDisp,
      nDoF, aBCFlag,1.0, vecx);
}

int main(
    [[maybe_unused]] int argc,
    [[maybe_unused]] char* argv[])
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aHex;
  dfm2::MeshHex3_Grid(
      aXYZ0, aHex,
      40, 20, 20, 0.1);
  std::vector<double> aMass(aXYZ0.size()/3);

  dfm2::CMatrixSparse<double> smat;
  delfem2::CPreconditionerILU<double> silu;
  InitializeMatrix(
      smat,silu,
      aHex,aXYZ0);

  std::vector<double> aDisp(aXYZ0.size(), 0.0);
  std::vector<int> aBCFlag(aXYZ0.size(), 0.0); // 0: free, 1: fix BC
  {
    for(unsigned int ip=0;ip<aXYZ0.size()/3;++ip){
      double x0 = aXYZ0[ip*3+0];
      if( x0 > 1.0e-10 ){ continue; }
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aBCFlag[ip*3+2] = 1;
    }
  }
  const double mass = 0.5;
  const double gravity[3] = {0,0,-10};
  Simulation(
      aDisp, smat, silu,
      aXYZ0, aHex, aBCFlag,
      mass, 1.0e+5,1.e+5, gravity);

  // ----------------------
  delfem2::glfw::CViewer3 viewer;
  viewer.camera.view_height = 1.5;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::ZTOP;
  viewer.camera.theta = 0.1;
  viewer.camera.psi = 0.1;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();

  delfem2::opengl::setSomeLighting();
  while(!glfwWindowShouldClose(viewer.window)){
    // -----
    viewer.DrawBegin_oldGL();
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,0,0);
    delfem2::opengl::DrawMeshHex3D_EdgeDisp(
        aXYZ0.data(),aXYZ0.size()/3,
        aHex.data(),aHex.size()/8,
        aDisp.data());
    //
    ::glEnable(GL_LIGHTING);
//    dfm2::opengl::DrawMeshHex3D_FaceNorm(aXYZ0.data(), aHex.data(), aHex.size() / 8);
    delfem2::opengl::DrawHex3D_FaceNormDisp(aXYZ0,aHex,aDisp);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
