/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include "delfem2/deflap.h"
#include "delfem2/mat4.h"
#include "delfem2/mshmisc.h"
#include "delfem2/primitive.h"
// ----------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// --------------------------------------------------

void SetPositionAtFixedBoundary
(std::vector<double>& aRhs,
 unsigned int iframe,
 const std::vector<double>& aXYZ0,
 const std::vector<int>& aBCFlag)
{
  double A[16];
  {
    dfm2::Mat4_Identity(A);
    const double trans0[3] = {0, -0.8, 0};
    dfm2::Translate_Mat4Affine(A,
                               trans0);
    const double axis0[3] = {0, +2.0*sin(0.03*iframe), 1.0*sin(0.07*iframe)};
    dfm2::Rotate_Mat4AffineRodriguez(A,
                                     axis0);
    const double trans1[3] = {0.2*sin(0.03*iframe), +0.6+0.2*cos(0.05*iframe), 0};
    dfm2::Translate_Mat4Affine(A,
                               trans1);
  }
  const unsigned int np = aRhs.size()/3;
  for(unsigned int ip=0;ip<np;++ip){
    if( aBCFlag[ip*3+0] == 0 ){ continue; }
    if( aBCFlag[ip*3+0] == 1 ){
      aRhs[ip*3+0] = aXYZ0[ip*3+0];
      aRhs[ip*3+1] = aXYZ0[ip*3+1];
      aRhs[ip*3+2] = aXYZ0[ip*3+2];
    }
    if( aBCFlag[ip*3+0] == 2 ) {
      dfm2::Vec3_Mat4Vec3_AffineProjection(aRhs.data()+ip*3, A, aXYZ0.data()+ip*3);
    }
  }
}

// ==================================

void myGlutDisplay
 (const std::vector<double>& aXYZ0,
  const std::vector<double>& aXYZ1,
  const std::vector<unsigned int>& aTri,
  const std::vector<int>& aBCFlag)
{
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1,aTri);
  ::glColor3d(0.8,0.8,0.8);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ0.data(),aXYZ0.size()/3,
                                   aTri.data(),aTri.size()/3);
  ::glColor3d(0,0,0);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ1.data(),aXYZ1.size()/3,
                                   aTri.data(),aTri.size()/3);
  
  { // draw bc as a point
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    for(unsigned int ip=0;ip<aXYZ0.size()/3;++ip){
      if( aBCFlag[ip*3+0] == 0 ){ continue; }
      if( aBCFlag[ip*3+0] == 1 ){ ::glColor3d(0,0,1); }
      if( aBCFlag[ip*3+0] == 2 ){ ::glColor3d(0,1,0); }
      ::glVertex3dv(aXYZ1.data()+ip*3);
    }
    ::glEnd();
  }
}

int main(int argc,char* argv[])
{
  std::vector<int> aBCFlag;
  std::vector<unsigned int> aTri;
  std::vector<double> aXYZ0;
  { // set problem
    dfm2::MeshTri3D_CylinderClosed(aXYZ0, aTri,
                                   0.2, 1.6,
                                   24, 24);
    aBCFlag.assign(aXYZ0.size(), 0);
    for(unsigned int ip=0;ip<aXYZ0.size()/3;++ip) {
      double y0 = aXYZ0[ip*3+1];
      if( y0 < -0.65 ){
        aBCFlag[ip*3+0] = 1;
        aBCFlag[ip*3+1] = 1;
        aBCFlag[ip*3+2] = 1;
      }
      if( y0 > +0.65 ){
        aBCFlag[ip*3+0] = 2;
        aBCFlag[ip*3+1] = 2;
        aBCFlag[ip*3+2] = 2;
      }
    }
  }
  std::vector<double> aXYZ1 = aXYZ0;
    
  // ------------------------
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  
  int iframe = 0;
  while (!glfwWindowShouldClose(viewer.window))
  {
    {
      glfwSetWindowTitle(viewer.window, "Def_LaplacianLinearDegenrate without preconditioner");
      dfm2::CDef_LaplacianLinearDegenerate def;
      def.Init(aXYZ0, aTri, false);
      def.aBCFlag = aBCFlag;
      for(;iframe<100;++iframe){
        SetPositionAtFixedBoundary(aXYZ1,
                                   iframe, aXYZ0, aBCFlag);
        def.Deform(aXYZ1,
                   aXYZ0);
        std::cout << "   laplacian linear degenerate without preconditioner : " << def.aConvHist.size() << std::endl;
        // -----
        viewer.DrawBegin_oldGL();
        myGlutDisplay(aXYZ0,aXYZ1,aTri,aBCFlag);
        viewer.SwapBuffers();
        glfwPollEvents();
        if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
      }
    }
    {
      glfwSetWindowTitle(viewer.window, "Def_LaplacianLinearDegenrate with preconditioner");
      dfm2::CDef_LaplacianLinearDegenerate def;
      def.Init(aXYZ0, aTri, true);
      def.aBCFlag = aBCFlag;
      def.SetBoundaryConditionToPreconditioner();
      for(;iframe<200;++iframe){
        SetPositionAtFixedBoundary(aXYZ1,
                                   iframe, aXYZ0, aBCFlag);
        def.Deform(aXYZ1,
                   aXYZ0);
        std::cout << "   laplacian linear degenerate with preconditioner: " << def.aConvHist.size() << std::endl;
        // -----
        viewer.DrawBegin_oldGL();
        myGlutDisplay(aXYZ0,aXYZ1,aTri,aBCFlag);
        viewer.SwapBuffers();
        glfwPollEvents();
        if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
      }
    }
    {
      glfwSetWindowTitle(viewer.window, "Def_LaplacianLinearGram without Preconditioner");
      dfm2::CDef_LaplacianLinearGram def;
      def.Init(aXYZ0, aTri, false);
      def.aBCFlag = aBCFlag;
      for(;iframe<300;++iframe){
        SetPositionAtFixedBoundary(aXYZ1,
                                   iframe, aXYZ0, aBCFlag);
        def.Deform(aXYZ1,
                   aXYZ0);
        std::cout << "  cg nitr:" << def.aConvHist.size() << std::endl;
        // -------
        viewer.DrawBegin_oldGL();
        myGlutDisplay(aXYZ0,aXYZ1,aTri,aBCFlag);
        viewer.SwapBuffers();
        glfwPollEvents();
        if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
      }
    }
    {
      dfm2::CDef_LaplacianLinearGram def;
      def.Init(aXYZ0, aTri, true);
      def.aBCFlag = aBCFlag;
      def.SetBoundaryConditionToPreconditioner();
      glfwSetWindowTitle(viewer.window, "Def_LaplacianLinearGram with Preconditioner");
      for(;iframe<400;++iframe){
        SetPositionAtFixedBoundary(aXYZ1,
                                   iframe, aXYZ0, aBCFlag);
        def.Deform(aXYZ1,
                   aXYZ0);
        std::cout << "  pcg nitr:" << def.aConvHist.size() << std::endl;
        // -------
        viewer.DrawBegin_oldGL();
        myGlutDisplay(aXYZ0,aXYZ1,aTri,aBCFlag);
        viewer.SwapBuffers();
        glfwPollEvents();
        if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
      }
    }
    {
      glfwSetWindowTitle(viewer.window, "Def_LaplacianLinear without Preconditioner");
      dfm2::CDef_LaplacianLinear def;
      def.Init(aXYZ0, aTri, false);
      def.aBCFlag = aBCFlag;
      def.SetValueToPreconditioner();
      for(;iframe<500;++iframe){
        SetPositionAtFixedBoundary(aXYZ1,
                                   iframe, aXYZ0, aBCFlag);
        def.Deform(aXYZ1,
                   aXYZ0);
        std::cout << "  cg nitr:" << def.aConvHist.size() << std::endl;
        // -------
        viewer.DrawBegin_oldGL();
        myGlutDisplay(aXYZ0,aXYZ1,aTri,aBCFlag);
        viewer.SwapBuffers();
        glfwPollEvents();
        if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
      }
    }
    {
      glfwSetWindowTitle(viewer.window, "Def_LaplacianLinear with Preconditioner");
      dfm2::CDef_LaplacianLinear def;
      def.Init(aXYZ0, aTri, true);
      def.aBCFlag = aBCFlag;
      def.SetValueToPreconditioner();
      for(;iframe<600;++iframe){
        SetPositionAtFixedBoundary(aXYZ1,
                                   iframe, aXYZ0, aBCFlag);
        def.Deform(aXYZ1,
                   aXYZ0);
        std::cout << "  pcg nitr:" << def.aConvHist.size() << std::endl;
        // -------
        viewer.DrawBegin_oldGL();
        myGlutDisplay(aXYZ0,aXYZ1,aTri,aBCFlag);
        viewer.SwapBuffers();
        glfwPollEvents();
        if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
      }
    }
    {
      glfwSetWindowTitle(viewer.window, "Direct Constraint");
      dfm2::CDef_LaplacianLinearAsym def;
      def.Init(aXYZ0, aTri);
      for(;iframe<700;++iframe){
        SetPositionAtFixedBoundary(aXYZ1,
                                   iframe, aXYZ0, aBCFlag);
        def.Deform(aXYZ1,
                   aXYZ0, aBCFlag);
        std::cout << "   bicgstab" << def.aHistConv.size() << std::endl;
        // -----
        viewer.DrawBegin_oldGL();
        myGlutDisplay(aXYZ0,aXYZ1,aTri,aBCFlag);
        viewer.SwapBuffers();
        glfwPollEvents();
        if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
      }
    }
    // ------
    iframe = 0;
  }
  
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


