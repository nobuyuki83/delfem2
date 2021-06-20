/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/eigen/ls_dense.h"
#include "delfem2/eigen/ls_sparse.h"
#include "delfem2/femsolidlinear.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/mshuni.h"
#include "delfem2/femutil.h"
#include <random>
#include <vector>
#include <cstdlib>
#include <GLFW/glfw3.h>
#include <Eigen/Core>

namespace dfm2 = delfem2;

// --------------------------------------------------------------

void Simulation_Mat3(
    std::vector<double>& aDisp,
    delfem2::CMatrixSparseBlock<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d>>& mA,
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
  mA.setZero();
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
      delfem2::Merge<8, 8, 3, 3, double>(mA, aIP, aIP, ddW, tmp_buffer);
    }
  }
  Eigen::VectorXd vec_b(nDoF);
  vec_b.setZero();
  for(unsigned int ip=0;ip<np;++ip) {
    vec_b[ip*3+0] += mass*gravity[0];
    vec_b[ip*3+1] += mass*gravity[1];
    vec_b[ip*3+2] += mass*gravity[2];
  }
  { // comput rhs vectors
    const Eigen::VectorXd& vd = Eigen::Map<const Eigen::VectorXd>(aDisp.data(),nDoF);
    AddMatVec(vec_b, 1.0, -1.0, mA, vd);
    std::cout << "energy" << vec_b.dot(vd) << std::endl;
  }
  SetFixedBC_Dia(mA, aBCFlag.data(), 1.f);
  SetFixedBC_Col(mA, aBCFlag.data());
  SetFixedBC_Row(mA, aBCFlag.data());
  delfem2::setZero_Flag(vec_b, aBCFlag,0);
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
      aDisp,
      aBCFlag,1.0, vec_x);
}

void Simulation_Mat4(
    std::vector<double>& aDisp,
    delfem2::CMatrixSparseBlock<Eigen::Matrix4d,Eigen::aligned_allocator<Eigen::Matrix4d>, 3>& mA,
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
  mA.setZero();
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
      delfem2::Merge<8, 8, 3, 3, double>(mA, aIP, aIP, ddW, tmp_buffer);
    }
  }
  Eigen::Matrix<double,-1,4,Eigen::RowMajor> vec_b(np,4);
  vec_b.setZero();
  for(unsigned int ip=0;ip<np;++ip) {
    vec_b(ip,0) += mass*gravity[0];
    vec_b(ip,1) += mass*gravity[1];
    vec_b(ip,2) += mass*gravity[2];
  }
  { // comput rhs vectors
    Eigen::Matrix<double,-1,4,Eigen::RowMajor>  vd(np, 4);
    for(unsigned int ip=0;ip<np;++ip){
      vd(ip, 0) = aDisp[ip*3+0];
      vd(ip, 1) = aDisp[ip*3+1];
      vd(ip, 2) = aDisp[ip*3+2];
      vd(ip, 3) = 0.0;
    }
    AddMatVec(vec_b, 1.0, -1.0, mA, vd);
    std::cout << "energy" << delfem2::Dot(vec_b, vd) << std::endl;
  }
  SetFixedBC_Dia(mA, aBCFlag.data(), 1.f);
  SetFixedBC_Col(mA, aBCFlag.data());
  SetFixedBC_Row(mA, aBCFlag.data());
  delfem2::setZero_Flag(vec_b, np,aBCFlag,0);
  // --------------------------------
  Eigen::Matrix<double,-1,4,Eigen::RowMajor> vec_x(np, 4);
  {
    double conv_ratio = 1.0e-6;
    int iteration = 1000;
    Eigen::Matrix<double,-1,4,Eigen::RowMajor> tmp0(np,4), tmp1(np,4);
    std::vector<double> aConv = delfem2::Solve_CG(
        vec_b, vec_x, tmp0, tmp1,
        conv_ratio, iteration, mA);
    std::cout << aConv.size() << std::endl;
  }
  // ------------------------------
  dfm2::XPlusAY(
      aDisp,
      np, aBCFlag,1.0, vec_x);
}


int main(int argc,char* argv[])
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aHex;
  dfm2::MeshHex3_Grid(
      aXYZ0, aHex,
      50, 25, 25, 0.1);
  std::vector<double> aMass(aXYZ0.size()/3);

  delfem2::CMatrixSparseBlock<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d>> A3;
  delfem2::CMatrixSparseBlock<Eigen::Matrix4d,Eigen::aligned_allocator<Eigen::Matrix4d>, 3> A4;
  {
    const unsigned int np = aXYZ0.size() / 3;
    std::vector<unsigned int> psup_ind, psup;
    dfm2::JArray_PSuP_MeshElem(
        psup_ind, psup,
        aHex.data(), aHex.size() / 8, 8,
        (int) aXYZ0.size() / 3);
    A3.Initialize(np);
    A3.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
    A4.Initialize(np);
    A4.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
  }

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
  aDisp.assign(aXYZ0.size(), 0.0);
  for(int i=0;i<aDisp.size();++i){
    if( aBCFlag[i] != 0 ){ continue; }
    aDisp[i] = (i%10)*1.0e-4;
  }
  Simulation_Mat3(
      aDisp, A3,
      aXYZ0, aHex, aBCFlag,
      mass, 1.0e+5, 1.e+5, gravity);

  aDisp.assign(aXYZ0.size(), 0.0);
  for(int i=0;i<aDisp.size();++i){
    if( aBCFlag[i] != 0 ){ continue; }
    aDisp[i] = (i%10)*1.0e-4;
  }
  Simulation_Mat4(
      aDisp, A4,
      aXYZ0, aHex, aBCFlag,
      mass, 1.0e+5, 1.e+5, gravity);

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
