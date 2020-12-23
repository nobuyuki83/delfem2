/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief demo for view navigation in 3D
 * @details this demo is for showing CViewer_GLFW funtinoalities
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/primitive.h"
#include "delfem2/points.h"
#include "delfem2/femem2.h"
namespace dfm2 = delfem2;

#include <Eigen/Eigenvalues>
#include <GLFW/glfw3.h>

void ShowEigen_Laplace(
    const std::vector<double>& aXY,
    const std::vector<unsigned int>& aQuad,
    double elen)
{
  const unsigned int np = aXY.size()/2;
  Eigen::MatrixXd A(np, np);
  A.setZero();

  double emat[4][4];
  dfm2::EMat_Poission2_QuadOrth(emat,elen,elen);

  for(unsigned int iq=0;iq<aQuad.size()/4;++iq){
    const unsigned int aIp[4] = {aQuad[iq*4+0], aQuad[iq*4+1], aQuad[iq*4+2], aQuad[iq*4+3]};
    for(unsigned int in=0;in<4;++in){
      for(unsigned int jn=0;jn<4;++jn) {
        unsigned int ip = aIp[in];
        unsigned int jp = aIp[jn];
        A(ip,jp) += emat[in][jn];
      }
    }
  }

  Eigen::EigenSolver<Eigen::MatrixXd> s(A);
  std::cout << s.eigenvalues() << std::endl;
  std::cout << s.eigenvectors() << std::endl;

  //
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();

  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_oldGL();
    dfm2::opengl::DrawMeshQuad2D_Edge(aXY,aQuad);
//    dfm2::opengl::DrawMeshQuad2_ScalarP1(aXY,aQuad,)
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
}

void ShowEigen_SolidLinear(
    const std::vector<double>& aXY,
    const std::vector<unsigned int>& aQuad,
    double elen)
{
  const unsigned int np = aXY.size()/2;
  Eigen::MatrixXd A(np*2, np*2);
  A.setZero();

  double emat[4][4][2][2];
  dfm2::EMat_SolidLinear2_QuadOrth_GaussInt(emat,elen,elen, 1.0, 1.0, 1);

  for(unsigned int iq=0;iq<aQuad.size()/4;++iq){
    const unsigned int aIp[4] = {aQuad[iq*4+0], aQuad[iq*4+1], aQuad[iq*4+2], aQuad[iq*4+3]};
    for(unsigned int in=0;in<4;++in){
      for(unsigned int jn=0;jn<4;++jn) {
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

  Eigen::EigenSolver<Eigen::MatrixXd> s(A);
  std::cout << s.eigenvalues() << std::endl;
  std::cout << s.eigenvectors() << std::endl;
  const auto& evec = s.eigenvectors();

  //
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();

  std::vector<double> aDisp(np*2,0.0);

  unsigned int iframe = 0;
  while (!glfwWindowShouldClose(viewer.window))
  {
    iframe = (iframe+1)%(np*2);
    std::cout << s.eigenvalues()(iframe).real() << std::endl;
    for(unsigned int i=0;i<np*2;++i){
      aDisp[i] = evec(i,iframe).real();
    }
    for(unsigned int i=0;i<10;++i) {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawMeshQuad2D_EdgeDisp(
          aXY.data(), aXY.size()/2,
          aQuad.data(), aQuad.size()/4,
          aDisp.data());
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
  dfm2::MeshQuad2D_Grid(aXY,aQuad,
      10,5);
  const double elen = 0.1;
  dfm2::Scale_Points(aXY.data(), aXY.size()/2,2, elen);
  ShowEigen_SolidLinear(aXY,aQuad,elen);
  ShowEigen_Laplace(aXY,aQuad,elen);
}
