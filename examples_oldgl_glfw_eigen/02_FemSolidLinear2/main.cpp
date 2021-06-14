/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/lsvecx.h"
#include "delfem2/lsmats.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/femsolidlinear.h"
#include "delfem2/mshuni.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/dtri.h"
#include "delfem2/eigen/ls_dense.h"
#include "delfem2/eigen/ls_sparse.h"
#include "delfem2/lsitrsol.h"
#include <GLFW/glfw3.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Core>

namespace dfm2 = delfem2;

// ------------------------------------------


// ------------------------------

void MakeMesh(
    std::vector<double>& aXY1,
    std::vector<unsigned int>& aTri1,
    std::vector<int>& loopIP_ind,
    std::vector<int>& loopIP,
    double len)
{
  std::vector< std::vector<double> > aaXY;
  {
    aaXY.resize(1);
    aaXY[0].push_back(-len); aaXY[0].push_back(-len);
    aaXY[0].push_back(-len); aaXY[0].push_back(+len);
    aaXY[0].push_back(+len); aaXY[0].push_back(+len);
    aaXY[0].push_back(+len); aaXY[0].push_back(-len);
  }
  // ---------------------------------
  std::vector<dfm2::CVec2d> aVec2;
  const double elen = 0.05;
  {
    JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                         aaXY);
    if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
      return;
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind,aVec2);
    if( elen > 10e-10 ){
      ResamplingLoop(loopIP_ind,loopIP,aVec2,
                     elen );
    }
  }
  {
    std::vector<dfm2::CDynPntSur> aPo2D;
    std::vector<dfm2::CDynTri> aETri;
    Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                   loopIP_ind,loopIP);
    if( elen > 1.0e-10 ){
      dfm2::CInputTriangulation_Uniform param(1.0);
      std::vector<int> aFlgPnt(aPo2D.size());
      std::vector<unsigned int> aFlgTri(aETri.size(),0);
      MeshingInside(
          aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
          aVec2.size(),0,elen, param);
    }
    MeshTri2D_Export(
        aXY1,aTri1,
        aVec2,aETri);
  }
  std::cout<<"  ntri;"<<aTri1.size()/3<<"  nXY:"<<aXY1.size()/2<<std::endl;
}

// -------------------------
void Solve0(
    std::vector<double>& aVal,
    const std::vector<double>& aXY1,
    const std::vector<unsigned int>& aTri1,
    const std::vector<int>& aBCFlag)
{
  const unsigned int np = aXY1.size()/2;
  const unsigned int nDoF = np*2;
  // -----------
  std::vector<unsigned int> psup_ind0, psup0;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind0, psup0,
      aTri1.data(), aTri1.size()/3, 3,
      aXY1.size()/2);
  // -------------
  dfm2::CMatrixSparse<double> mat_A;
  mat_A.Initialize(np, 2, true);
  mat_A.SetPattern(psup_ind0.data(), psup_ind0.size(), psup0.data(), psup0.size());
  // ----------------------
  double myu = 10.0;
  double lambda = 10.0;
  double rho = 1.0;
  double g_x = 0.0;
  double g_y = -3.0;
  mat_A.setZero();
  std::vector<double> vec_b(nDoF, 0.0);
  dfm2::MergeLinSys_SolidLinear_Static_MeshTri2D(
      mat_A,vec_b.data(),
      myu,lambda,rho,g_x,g_y,
      aXY1.data(), aXY1.size()/2,
      aTri1.data(), aTri1.size()/3,
      aVal.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag,0);
  // ---------------
  std::vector<double> vec_x(vec_b.size());
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    std::vector<double> aConv = Solve_CG(
        dfm2::CVecXd(vec_b), dfm2::CVecXd(vec_x), dfm2::CVecXd(tmp0), dfm2::CVecXd(tmp1),
        conv_ratio, iteration, mat_A);
    std::cout << aConv.size() << std::endl;
  }
//  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  // --------------
  dfm2::XPlusAY(aVal,nDoF,aBCFlag,
          1.0,vec_x);
}

// -------------------------
void Solve1(
    std::vector<double>& aVal,
    const std::vector<double>& aXY1,
    const std::vector<unsigned int>& aTri1,
    const std::vector<int>& aBCFlag)
{
  const unsigned int np = aXY1.size()/2;
  const unsigned int nDoF = np*2;
  // -----------
  std::vector<unsigned int> psup_ind0, psup0;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind0, psup0,
      aTri1.data(), aTri1.size()/3, 3,
      aXY1.size()/2);
  // -------------
  delfem2::CMatrixSparseBlock<Eigen::Matrix2d,Eigen::aligned_allocator<Eigen::Matrix2d>> mA;
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
      mA,vec_b.data(),
      myu,lambda,rho,g_x,g_y,
      aXY1.data(), aXY1.size()/2,
      aTri1.data(), aTri1.size()/3,
      aVal.data());
  SetFixedBC_Dia(mA, aBCFlag.data(), 1.f);
  SetFixedBC_Col(mA, aBCFlag.data());
  SetFixedBC_Row(mA, aBCFlag.data());
  delfem2::setZero_Flag(vec_b, aBCFlag,0);
  // ---------------
  Eigen::VectorXd vec_x(vec_b.size());
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  {
    const std::size_t n = vec_b.size();
    Eigen::VectorXd tmp0(n), tmp1(n);
    std::vector<double> aConv = delfem2::Solve_CG(
        vec_b, vec_x, tmp0, tmp1,
        conv_ratio, iteration, mA);
    std::cout << aConv.size() << std::endl;
  }
//  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  // --------------
  delfem2::XPlusAY(aVal,nDoF,aBCFlag,
      1.0,vec_x);
}

int main(int argc,char* argv[])
{
  std::vector<unsigned int> aTri1;
  std::vector<double> aXY1;
  std::vector<int> loopIP_ind, loopIP; // vtx on loop
  double len = 1.1;
  MakeMesh(
      aXY1, aTri1, loopIP_ind, loopIP,
      len);
  // ---
  std::vector<double> aVal;
  {
    const unsigned int np = aXY1.size()/2;
    const unsigned int nDoF = np*2;
    std::vector<int> aBCFlag; // master slave flag
    aBCFlag.assign(nDoF, 0);
    for(unsigned int ip=0;ip<np;++ip){
//    const double px = aXY1[ip*2+0];
      const double py = aXY1[ip*2+1];
      if( fabs(py-len) > 0.0001 ){ continue; }
      aBCFlag[ip*2+0] = 1;
      aBCFlag[ip*2+1] = 1;
    }
    aVal.assign(np * 2, 0.0);
    Solve1(aVal,aXY1,aTri1,aBCFlag);
    aVal.assign(np * 2, 0.0);
    Solve0(aVal,aXY1,aTri1,aBCFlag);
  }
  // --------
  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 1.5;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  // ---------
  while(!::glfwWindowShouldClose(viewer.window)){
    viewer.DrawBegin_oldGL();
    delfem2::opengl::DrawMeshTri2D_FaceDisp2D(
        aXY1.data(), aXY1.size()/2,
        aTri1.data(), aTri1.size()/3,
        aVal.data(), 2);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
}
