/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/defarap.h"
#include "delfem2/mat4.h"
#include "delfem2/primitive.h"
#include <GLFW/glfw3.h>
#include <cmath>

namespace dfm2 = delfem2;

// -------------------------------------

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
    const double trans1[3] = {0.2*sin(0.03*iframe), +0.5+0.1*cos(0.05*iframe), 0};
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

void myGlutDisplay_Mesh
(const std::vector<double>& aXYZ0,
const std::vector<double>& aXYZ1,
const std::vector<unsigned int>& aTri)
{
  ::glLineWidth(1);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1,aTri);
  ::glColor3d(0.8,0.8,0.8);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ0.data(),aXYZ0.size()/3,
                                   aTri.data(),aTri.size()/3);
  ::glColor3d(0,0,0);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ1.data(),aXYZ1.size()/3,
                                   aTri.data(),aTri.size()/3);
  
}

void Draw_BCFlag(const std::vector<double>& aXYZ1,
                 const std::vector<int>& aBCFlag)
{ // draw bc as a point
  ::glPointSize(10);
  ::glBegin(GL_POINTS);
  for(unsigned int ip=0;ip<aXYZ1.size()/3;++ip){
    if( aBCFlag[ip*3+0] == 0 ){ continue; }
    if( aBCFlag[ip*3+0] == 1 ){ ::glColor3d(0,0,1); }
    if( aBCFlag[ip*3+0] == 2 ){ ::glColor3d(0,1,0); }
    ::glVertex3dv(aXYZ1.data()+ip*3);
  }
  ::glEnd();
}


// --------------------------------------------------

int main(int argc,char* argv[])
{
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.camera.view_height = 1.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  
  std::vector<unsigned int> aTri;
  std::vector<double> aXYZ0;
  std::vector<int> aBCFlag;
  {
    dfm2::MeshTri3D_CylinderClosed(
        aXYZ0, aTri,
        0.2, 1.6,
        16, 16);
    {
      const unsigned int np = aXYZ0.size() / 3;
      aBCFlag.assign(np * 3, 0);
      for(unsigned int ip=0;ip<np;++ip) {
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
  }
  std::vector<double> aXYZ1 = aXYZ0;

  for(unsigned int itr=0;itr<3;++itr){
    const double weight_bc = 100.0;
    int iframe = 0;
    { // arap edge linear disponly
      dfm2::CDef_ArapEdgeLinearDisponly def0(aXYZ0, aTri, weight_bc, aBCFlag);
      glfwSetWindowTitle(viewer.window, "(1) ARAP Edge Linear Disponly");
      for(;iframe<50;++iframe)
      {
        SetPositionAtFixedBoundary(aXYZ1,
                                   iframe,aXYZ0,aBCFlag);
        def0.Deform(aXYZ1,
                    aXYZ0);
        // --------------------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(aXYZ0,aXYZ1,aTri);
        Draw_BCFlag(aXYZ1,aBCFlag);
        viewer.SwapBuffers();
        glfwPollEvents();
        viewer.ExitIfClosed();
      }
    } // end linear disponly
    // -------------------------------------------------------
    { // begin lienar disprot without preconditioner
      glfwSetWindowTitle(viewer.window, "(2) Arap Edge Linear Disprot without Prec");
      unsigned int np = aXYZ0.size()/3;
      dfm2::CDef_ArapEdge def1;
      def1.Init(aXYZ0, aTri, weight_bc, aBCFlag, false);
      std::vector<double> aQuat(np*4); // array of quaternion
      for(;iframe<100;++iframe){
        for(unsigned int ip=0;ip<np;++ip){ dfm2::Quat_Identity(aQuat.data()+ip*4); } // Initialize
        SetPositionAtFixedBoundary(aXYZ1,
                                   iframe,aXYZ0,aBCFlag);
        def1.Deform(aXYZ1, aQuat,
                    aXYZ0);
        // ------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(aXYZ0,aXYZ1, aTri);
        Draw_BCFlag(aXYZ1,aBCFlag);
        dfm2::opengl::Draw_QuaternionsCoordinateAxes(aXYZ1,aQuat,0.04);
        viewer.SwapBuffers();
        glfwPollEvents();
        viewer.ExitIfClosed();
      } // end of frame loop
    } // end linear disprot without preconditioner
    // -------------------------------
    { // begin lienar disprot with preconditioner
      glfwSetWindowTitle(viewer.window, "(3) Arap Edge Linear Disprot with Prec");
      const unsigned int np = aXYZ0.size()/3;
      dfm2::CDef_ArapEdge def1;
      def1.Init(aXYZ0, aTri, weight_bc, aBCFlag, true);
      std::vector<double> aQuat(np*4);
      for(;iframe<200;++iframe){
        for(unsigned int ip=0;ip<np;++ip){ dfm2::Quat_Identity(aQuat.data()+ip*4); } // initialize
        for(unsigned int i=0;i<np*3;++i){ // adding noise for debuggng purpose
          aXYZ1[i] += 0.02*(double)rand()/(RAND_MAX+1.0)-0.01;
        }
        SetPositionAtFixedBoundary(
            aXYZ1,
            iframe,aXYZ0,aBCFlag);
        def1.Deform(
            aXYZ1, aQuat,
            aXYZ0);
        // ------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(aXYZ0,aXYZ1, aTri);
        Draw_BCFlag(aXYZ1,aBCFlag);
        dfm2::opengl::Draw_QuaternionsCoordinateAxes(aXYZ1,aQuat,0.04);
        viewer.SwapBuffers();
        glfwPollEvents();
        viewer.ExitIfClosed();
      } // end of frame loop
    } // end linear disprot with preconditioner
    // -------------------------------
    { // begin nonlienar disprot with preconditioner
      glfwSetWindowTitle(viewer.window, "(4) Arap Edge NonLinear Disprot with Prec");
      const unsigned int np = aXYZ0.size()/3;
      dfm2::CDef_ArapEdge def1;
      def1.Init(aXYZ0, aTri, weight_bc, aBCFlag, true);
      std::vector<double> aQuat(np*4);
      for(unsigned int ip=0;ip<np;++ip){ dfm2::Quat_Identity(aQuat.data()+ip*4); }
      for(;iframe<400;++iframe){
        for(unsigned int i=0;i<np*3;++i){ // adding noise for debuggng purpose
          aXYZ1[i] += 0.02*(double)rand()/(RAND_MAX+1.0)-0.01;
        }
        SetPositionAtFixedBoundary(
            aXYZ1,
            iframe,aXYZ0,aBCFlag);
        def1.Deform(
            aXYZ1, aQuat,
            aXYZ0);
        // ------
        viewer.DrawBegin_oldGL();
        myGlutDisplay_Mesh(aXYZ0,aXYZ1, aTri);
        Draw_BCFlag(aXYZ1,aBCFlag);
        dfm2::opengl::Draw_QuaternionsCoordinateAxes(aXYZ1,aQuat,0.04);
        viewer.SwapBuffers();
        glfwPollEvents();
        viewer.ExitIfClosed();
      } // end of frame loop
    } // end linear disprot with preconditioner
  }
}


