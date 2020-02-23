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
#include "delfem2/srchuni_v3.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_color.h"
#include "delfem2/opengl/glold_v23.h"
#include "delfem2/opengl/glold_smplr.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

std::vector<double> aXYZ;
std::vector<unsigned int> aTri;

// ------------------------------------------------------

void DrawObject(){
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
}

int main(int argc,char* argv[])
{
  dfm2::Read_Obj(std::string(PATH_INPUT_DIR)+"/rollsRoyce.obj",
                 aXYZ,aTri);
  dfm2::Normalize_Points3(aXYZ,4.0);
  // ---------------------------------------
  
  dfm2::opengl::CRender2Tex_DrawOldGL sampler;
  {
    unsigned int nresX = 128;
    unsigned int nresY = 128;
    unsigned int nresZ = 256;
    double elen = 0.02;
    sampler.SetTextureProperty(nresX, nresZ, true);
    sampler.SetCoord(elen, elen*nresY,
                         dfm2::CVec3d(-0.5*elen*nresX,+0.5*elen*nresY,+0.5*elen*nresZ).stlvec(),
                         dfm2::CVec3d(0,+1,0).stlvec(),
                         dfm2::CVec3d(1,+0,0).stlvec() );
    sampler.SetPointColor(0.0, 1.0, 0.0);
    sampler.draw_len_axis = 0.2;
    sampler.isDrawTex = false;
    sampler.isDrawOnlyHitPoints = true;
  }
  // ---------------------------------------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 2.0;
  viewer.nav.camera.camera_rot_mode = dfm2::CAMERA_ROT_TBALL;
  viewer.nav.camera.Rot_Camera(+0.2, -0.2);
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }

  dfm2::opengl::setSomeLighting();
  ::glEnable(GL_DEPTH_TEST);
  
  {
    sampler.InitGL(); // move the sampled image to a texture
    sampler.Start();
    ::glClearColor(1.0, 1.0, 1.0, 1.0 );
    ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    ::glEnable(GL_DEPTH_TEST);
    ::glDisable(GL_BLEND);
    ::glEnable(GL_LIGHTING);
    DrawObject();
    sampler.End();
    sampler.GetDepth();
    sampler.GetColor();
  }

  for(int iframe=0;;iframe++){
    dfm2::CVec3d dir(cos(iframe*0.003),sin(iframe*0.005),sin(iframe*0.007));
    dfm2::CVec3d src(1.5*cos(iframe*0.005), 1.5*sin(iframe*0.007), 1.5*sin(iframe*0.009) );
    // ----
    std::vector<double> aXYZ1;
    {
      const dfm2::CVec3d ex(sampler.x_axis);
      const dfm2::CVec3d ez(sampler.z_axis);
      const dfm2::CVec3d ey = dfm2::Cross(ez,ex);
      const double lsrc[3] = {
        (src-dfm2::CVec3d(sampler.origin))*ex - 0.5*sampler.lengrid,
        (src-dfm2::CVec3d(sampler.origin))*ey - 0.5*sampler.lengrid,
        (src-dfm2::CVec3d(sampler.origin))*ez };
      const double ldir[3] = { dir*ex, dir*ey, dir*ez };
      std::vector<dfm2::CPointElemSurf<double>> aPES;
      dfm2::IntersectionLine_Hightfield(aPES,
          -sampler.z_range*0.999, -sampler.z_range*0.001,
          lsrc, ldir,
          sampler.nResX, sampler.nResY,
          sampler.lengrid,
          sampler.aZ);
      for(auto & pes : aPES){
        dfm2::CVec3d lpos = pes.Pos_Grid(sampler.nResX, sampler.nResY, sampler.lengrid, sampler.aZ);
        dfm2::CVec3d offset = dfm2::CVec3d(sampler.origin) + ex * sampler.lengrid*0.5 + ey * sampler.lengrid*0.5;
        dfm2::CVec3d gp = ex*lpos.x() + ey*lpos.y() + ez*lpos.z() + offset;
        dfm2::CVec3d gp1 = dfm2::nearest_Line_Point(gp, src, dir);
        aXYZ1.push_back(gp.x());
        aXYZ1.push_back(gp.y());
        aXYZ1.push_back(gp.z());
      }
    }
    {
      viewer.DrawBegin_oldGL();
      dfm2::opengl::DrawBackground( dfm2::CColor(0.2,0.7,0.7) );
      // ----------
      glPointSize(5);
      {
        ::glBegin(GL_POINTS);
        ::glColor3d(1,0,0);
        ::glVertex3d(src.x(), src.y(), src.z());
        ::glEnd();
      }
      for(int ixyz=0;ixyz<aXYZ1.size()/3;++ixyz){
        ::glBegin(GL_POINTS);
        ::glColor3d(0,0,1);
        ::glVertex3d(aXYZ1[ixyz*3+0], aXYZ1[ixyz*3+1], aXYZ1[ixyz*3+2]);
        ::glEnd();
      }
      {
        dfm2::CVec3d p0 = src + 10*dir;
        dfm2::CVec3d p1 = src - 10*dir;
        ::glLineWidth(1);
        ::glBegin(GL_LINES);
        ::glVertex3d(p0.x(), p0.y(), p0.z());
        ::glVertex3d(p1.x(), p1.y(), p1.z());
        ::glEnd();
      }
      // ----------
      ::glEnable(GL_LIGHTING);
      DrawObject();
      glPointSize(1);
      sampler.Draw();
      viewer.DrawEnd_oldGL();
    }
    if( glfwWindowShouldClose(viewer.window) ){ break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


