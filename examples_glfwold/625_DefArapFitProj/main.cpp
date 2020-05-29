/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshio.h"
#include "delfem2/def.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/r2tglo_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

bool GetProjectedPoint
(dfm2::CVec3d& p0,
 const dfm2::CVec3d& ps,
 const dfm2::opengl::CRender2Tex_DrawOldGL& smplr)
{
  const unsigned int nx = smplr.nResX;
  const dfm2::CVec3d& dx = smplr.x_axis;
  const dfm2::CVec3d& dz = smplr.z_axis;
  const dfm2::CVec3d& dy = dfm2::Cross(dz,dx);
  double lx = (ps-dfm2::CVec3d(smplr.origin))*dx;
  double ly = (ps-dfm2::CVec3d(smplr.origin))*dy;
  int ix0 = (int)floor(lx/smplr.lengrid-0.5);
  int iy0 = (int)floor(ly/smplr.lengrid-0.5);
  int ix1 = ix0+1;
  int iy1 = iy0+1;
  if( smplr.aZ[iy0*nx+ix0] < -smplr.z_range*0.99 ) return false;
  if( smplr.aZ[iy0*nx+ix1] < -smplr.z_range*0.99 ) return false;
  if( smplr.aZ[iy1*nx+ix0] < -smplr.z_range*0.99 ) return false;
  if( smplr.aZ[iy1*nx+ix1] < -smplr.z_range*0.99 ) return false;
  double rx = lx/smplr.lengrid-ix0;
  double ry = ly/smplr.lengrid-iy0;
  double z01 = (1-rx)*(1-ry)*smplr.aZ[iy0*nx+ix0]
  + rx*(1-ry)*smplr.aZ[iy0*nx+ix1]
  + (1-rx)*ry*smplr.aZ[iy1*nx+ix0]
  + rx*ry*smplr.aZ[iy1*nx+ix1];
  p0 = dfm2::CVec3d(smplr.origin) + lx*dx+ly*dy+dz*z01;
  return true;
}


int main(int argc,char* argv[])
{
  dfm2::opengl::CRender2Tex_DrawOldGL_BOX sampler_box;
  sampler_box.Initialize(128, 128, 128, 0.01);
  
  
  double pmin0[3]={1,0,0}, pmax0[3]={-1,0,0};
  sampler_box.BoundingBox3(pmin0, pmax0);
  
  // ---------------------------------------
  

  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  {
    {
      std::vector<double> aXYZ0;
      std::vector<unsigned int> aTri0;
      dfm2::Read_Obj(std::string(PATH_INPUT_DIR)+"/../car_0000_full_1.obj",
                     aXYZ0,aTri0);
      std::vector<int> aMap01;
      dfm2::RemoveUnreferencedPoints_MeshElem(aXYZ, aTri,
                                              aMap01, 3, aXYZ0, aTri0);
    }
    double pmin1[3], pmax1[3];
    dfm2::BoundingBox3_Points3(pmin1,pmax1,
                               aXYZ.data(), aXYZ.size()/3);
    dfm2::Translate_Points3(aXYZ,
                            -(pmin1[0]+pmax1[0])*0.5,
                            -(pmin1[1]+pmax1[1])*0.5,
                            -(pmin1[2]+pmax1[2])*0.5);
    double len0 = (pmax0[0]-pmin0[0]) + (pmax0[1]-pmin0[1]) + (pmax0[2]-pmin0[2]);
    double len1 = (pmax1[0]-pmin1[0]) + (pmax1[1]-pmin1[1]) + (pmax1[2]-pmin1[2]);
    double scale = len0/len1;
    dfm2::Scale_Points3(aXYZ.data(),
                        aXYZ.size()/3, scale);
  }
 
  // ---------------------------------------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 2.0;
  viewer.nav.camera.camera_rot_mode = dfm2::CCamera<double>::CAMERA_ROT_MODE::TBALL;
  viewer.nav.camera.Rot_Camera(+0.2, -0.2);
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }
  dfm2::opengl::setSomeLighting();
  ::glEnable(GL_DEPTH_TEST);
  
  std::vector<double> aXYZ0 = aXYZ;
  dfm2::CDef_LaplacianLinearGram def;
  def.Init(aXYZ0, aTri, true);
  std::vector<int> aBCFlag;

  // ---------------------

  while (true)
  {
    aBCFlag.assign(aXYZ0.size(),0);
    for(unsigned int ip=0;ip<aXYZ.size()/3;++ip){
      const dfm2::CVec3d& ps = dfm2::CVec3d(aXYZ.data()+ip*3);
      bool is_valid = false;
      dfm2::CVec3d pmin;
      for(unsigned int ismplr=0;ismplr<sampler_box.aSampler.size();++ismplr){
        const auto& smplr = sampler_box.aSampler[ismplr];
        dfm2::CVec3d p0;
        bool res = GetProjectedPoint(p0, ps,smplr);
        if( !res ){ continue; }
        double len = (p0-ps).Length();
        if( len > 0.1 ){ continue; }
        if( !is_valid ){ is_valid = true; pmin = p0; continue; }
        if( (p0-ps).Length() < (pmin-ps).Length() ){
          pmin = p0;
        }
      }
      if( !is_valid ){ continue; }
      aXYZ[ip*3+0] = pmin.x();
      aXYZ[ip*3+1] = pmin.y();
      aXYZ[ip*3+2] = pmin.z();
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aBCFlag[ip*3+2] = 1;
    }
    def.SetBoundaryCondition(aBCFlag);
    def.Deform(aXYZ, aXYZ0);
    // ----------------------------
    for(int iframe=0;iframe<100;++iframe){
      // ----
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawBackground( dfm2::CColor(0.2,0.7,0.7) );
      ::glEnable(GL_LIGHTING);
      ::glColor3d(1,1,1);
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
      ::glDisable(GL_LIGHTING);
      ::glColor3d(0,0,0);
      ::glLineWidth(1);
      dfm2::opengl::DrawMeshTri3D_Edge(aXYZ,aTri);
      dfm2::opengl::DrawBox3_Edge(pmin0, pmax0);
      glPointSize(3);
      sampler_box.Draw();
      glfwSwapBuffers(viewer.window);
      glfwPollEvents();
      if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
    }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


