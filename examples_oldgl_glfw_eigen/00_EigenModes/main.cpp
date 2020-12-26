/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief demo for visualizing eigenmodes
 */

#include "delfem2/eigen/lsitrsol.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/mshuni.h"
#include "delfem2/points.h"
#include "delfem2/femsolidlinear.h"
namespace dfm2 = delfem2;

#include <Eigen/Eigenvalues>
#include <GLFW/glfw3.h>

void ShowEigen_SolidLinear_MeshQuad2(
    const std::vector<double>& aXY,
    const std::vector<unsigned int>& aQuad,
    double elen) {
  const unsigned int np = aXY.size() / 2;
  Eigen::MatrixXd A(np * 2, np * 2);
  A.setZero();

  double emat[4][4][2][2];
  dfm2::EMat_SolidLinear2_QuadOrth_GaussInt(emat, elen, elen, 1.0, 1.0, 2);

  for (unsigned int iq = 0; iq < aQuad.size() / 4; ++iq) {
    const unsigned int aIp[4] = {
        aQuad[iq * 4 + 0], aQuad[iq * 4 + 1], aQuad[iq * 4 + 2], aQuad[iq * 4 + 3]};
    for (unsigned int in = 0; in < 4; ++in) {
      for (unsigned int jn = 0; jn < 4; ++jn) {
        const unsigned int ip = aIp[in];
        const unsigned int jp = aIp[jn];
        for (unsigned int idim = 0; idim < 2; ++idim) {
          for (unsigned int jdim = 0; jdim < 2; ++jdim) {
            A(ip * 2 + 0, jp * 2 + 0) += emat[in][jn][0][0];
            A(ip * 2 + 0, jp * 2 + 1) += emat[in][jn][0][1];
            A(ip * 2 + 1, jp * 2 + 0) += emat[in][jn][1][0];
            A(ip * 2 + 1, jp * 2 + 1) += emat[in][jn][1][1];
          }
        }
      }
    }
  }

  {
    auto B = A;
    B(0,0) += 1.0;
    B(1,1) += 1.0;
    B(2,2) += 1.0;
    B(3,3) += 1.0;
    Eigen::VectorXd r(np*2), u(np*2), Ap(np*2), p(np*2);
    r.setRandom();
    std::vector<double> aConv = dfm2::Solve_CG(r,u,Ap,p,1.0e-5,300,B);
    std::cout << aConv.size() << std::endl;
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> s(A);
  const auto& evec = s.eigenvectors(); // somehow this takes time

  //
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();

  std::vector<double> aDisp(np*2,0.0);

  for(unsigned int iframe=0;iframe<10;++iframe){
    std::cout << s.eigenvalues()(iframe) << std::endl;
    for(unsigned int i=0;i<np*2;++i){
      aDisp[i] = evec(i,iframe);
    }
    for(unsigned int i=0;i<100;++i) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshQuad2D_EdgeDisp(
          aXY.data(), aXY.size()/2,
          aQuad.data(), aQuad.size()/4,
          aDisp.data(),sin(i*0.1));
      viewer.SwapBuffers();
      glfwPollEvents();
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
}

void ShowEigen_SolidLinear_MeshTri2(
    const std::vector<double>& aXY,
    const std::vector<unsigned int>& aTri)
{
  const unsigned int np = aXY.size()/2;
  Eigen::MatrixXd A(np*2, np*2);
  A.setZero();

  std::vector<double> aDisp(np*2,0.0);

  for(unsigned int it=0;it<aTri.size()/3;++it){
    const unsigned int aIp[3] = {aTri[it*3+0], aTri[it*3+1], aTri[it*3+2]};
    double disp[3][2]; dfm2::FetchData(&disp[0][0],3,2,aIp,aDisp.data());
    double coord[3][2]; dfm2::FetchData(&coord[0][0],3,2,aIp,aXY.data());
    double emat[3][3][2][2], eres[3][2];
    dfm2::EMat_SolidStaticLinear_Tri2D(eres,emat,1,1,0,0,0,disp,coord);
    for(unsigned int in=0;in<3;++in){
      for(unsigned int jn=0;jn<3;++jn) {
        const unsigned int ip = aIp[in];
        const unsigned int jp = aIp[jn];
        for(unsigned int idim=0;idim<2;++idim){
          for(unsigned int jdim=0;jdim<2;++jdim) {
            A(ip*2+0,jp*2+0) += emat[in][jn][0][0];
            A(ip*2+0,jp*2+1) += emat[in][jn][0][1];
            A(ip*2+1,jp*2+0) += emat[in][jn][1][0];
            A(ip*2+1,jp*2+1) += emat[in][jn][1][1];
          }
        }
      }
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> s(A);
  const auto& evec = s.eigenvectors();
  //
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  for(unsigned int iframe=0;iframe<10;++iframe){
    std::cout << s.eigenvalues()(iframe) << std::endl;
    for(unsigned int i=0;i<np*2;++i){
      aDisp[i] = evec(i,iframe);
    }
    for(unsigned int i=0;i<100;++i) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshTri2D_EdgeDisp(
          aXY.data(), aXY.size()/2,
          aTri.data(), aTri.size()/3,
          aDisp.data(),sin(i*0.1));
      viewer.SwapBuffers();
      glfwPollEvents();
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
}

int main()
{
  std::vector<double> aXY;
  std::vector<unsigned int> aQuad;
  unsigned int nx = 10;
  unsigned int ny = 5;
  dfm2::MeshQuad2D_Grid(
      aXY,aQuad,
      nx,ny);
  const double elen = 0.1;
  {
    dfm2::Scale_Points(aXY.data(), aXY.size() / 2, 2, elen);
    dfm2::Translate_Points2(aXY,-0.5*nx*elen, -0.5*ny*elen);
  }
  ShowEigen_SolidLinear_MeshQuad2(aXY, aQuad, elen);
  //
  std::vector<unsigned int> aTri;
  dfm2::convert2Tri_Quad(
      aTri,
      aQuad);
  ShowEigen_SolidLinear_MeshTri2(aXY,aTri);
}
