/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <iostream>
#include <cmath>
#include "delfem2/def.h"
#include "delfem2/gizmo_geo3.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/primitive.h"
#include "delfem2/vec3.h"
#include "delfem2/quat.h"
#include "delfem2/mat4.h"
#include "delfem2/geo3_v23m34q.h"
// ---
#include <GLFW/glfw3.h>
#include "delfem2/opengl/gizmo_glold.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// -------------------------------------------

int main(int argc,char* argv[])
{
  class CMyViewer : public delfem2::opengl::CViewer_GLFW {
  public:
    CMyViewer(){
      dfm2::MeshTri3D_CylinderClosed(aXYZ0, aTri,
                                     0.2, 1.6,
                                     32, 32);
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
      { // initialize gizmo
        unsigned int nbc = 2;
        std::vector<dfm2::CVec3d> aCntBC;
        aCntBC.assign(nbc, dfm2::CVec3d(0,0,0) );
        std::vector<unsigned int> aW(nbc, 0);
        for(unsigned int ip=0;ip<np;++ip){
          if( aBCFlag[ip*3+0] <= 0 ){ continue; }
          unsigned int ibc = aBCFlag[ip*3+0]-1;
          aCntBC[ibc] += dfm2::CVec3d(aXYZ0.data()+ip*3);
          aW[ibc] += 1;
        }
        for(unsigned int ibc=0;ibc<nbc;++ibc){
          aCntBC[ibc] /= (double)aW[ibc];
        }
        giz1.pivot0 = aCntBC[1].Float();
        giz1.gizmo_rot.pos = aCntBC[1].Float();
        giz1.gizmo_trnsl.pos = aCntBC[1].Float();
        giz1.gizmo_rot.size = 0.3;
        giz1.gizmo_trnsl.size = 0.3;
      }
      aXYZ1 = aXYZ0;
      aQuat1.resize(np*4);
      for(unsigned int ip=0;ip<np;++ip){
        dfm2::Quat_Identity(aQuat1.data()+4*ip);
      }
      def0.Init(aXYZ0, aTri, true);
    }
    virtual void mouse_press(const float src[3], const float dir[3]){
      giz1.Pick(src, dir);
    }
    virtual void mouse_drag(const float src0[3], const float src1[3], const float dir[3]){
      giz1.Drag(src0, src1, dir);
    }
    virtual void key_release(int key, int mods){
    }
    virtual void key_press(int key, int mods){
      if( key == GLFW_KEY_R ){ giz1.igizmo_mode = 1; }
      if( key == GLFW_KEY_G ){ giz1.igizmo_mode = 0; }
    }
    //
    void Draw(){
      { // set boundary condition
        const dfm2::CMat4<double> aff1 = giz1.Affine().Double();
        for(unsigned int ip=0;ip<aXYZ0.size()/3;++ip){
          if( aBCFlag[ip*3+0] == 0 ){ continue; }
          if( aBCFlag[ip*3+0] == 1 ){
            aXYZ1[ip*3+0] = aXYZ0[ip*3+0];
            aXYZ1[ip*3+1] = aXYZ0[ip*3+1];
            aXYZ1[ip*3+2] = aXYZ0[ip*3+2];
          }
          if( aBCFlag[ip*3+0] == 2 ) {
            dfm2::Vec3_Mat4Vec3_AffineProjection(aXYZ1.data()+ip*3,
                                                 aff1.mat,
                                                 aXYZ0.data()+ip*3);
          }
        }
        for(int itr=0;itr<2;++itr){
          dfm2::UpdateRotationsByMatchingCluster_Linear(aQuat1,
                                                 aXYZ0,aXYZ1,
                                                 def0.psup_ind,def0.psup);
        }
      }
      def0.Deform(aXYZ1,aQuat1,
                  aXYZ0,aBCFlag);
      // -------------------------------
      DrawBegin_oldGL();
      delfem2::opengl::DrawAxis(1);
      { // mesh
        ::glEnable(GL_LIGHTING);
        ::glColor3d(0,0,0);
        delfem2::opengl::DrawMeshTri3D_Edge(aXYZ1.data(), aXYZ1.size()/3,
                                            aTri.data(), aTri.size()/3);
        delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(),
                                                aTri.data(), aTri.size()/3);
      }
      { // draw bc
        ::glDisable(GL_LIGHTING);
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
      /*
      for(int ibc=0;ibc<aCntBC.size();++ibc){
        dfm2::CVec3d p0 = aCntBC[ibc];
        dfm2::opengl::DrawSphereAt(32, 32, 0.1, p0.x(), p0.y(), p0.z());
        
      }
       */
      delfem2::opengl::Draw(giz1);
      DrawEnd_oldGL();
    }
  public:
    delfem2::CGizmo_Affine<float> giz1; // bcflag==2
    std::vector<unsigned int> aTri;
    std::vector<double> aXYZ0, aXYZ1;
    std::vector<double> aQuat1;
    std::vector<int> aBCFlag;
    delfem2::CDef_Arap def0;
  } viewer;
  // --------------------
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  // --------------------
  while(true){
    viewer.Draw();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


