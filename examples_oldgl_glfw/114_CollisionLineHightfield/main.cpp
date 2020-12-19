/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <glad/glad.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/r2tglo_glold.h"
#include "delfem2/srchuni_v3.h"
#include "delfem2/mshio.h"
#include "delfem2/points.h"
#include <GLFW/glfw3.h>
#include <cmath>

namespace dfm2 = delfem2;

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  dfm2::Read_Obj(
      std::string(PATH_INPUT_DIR)+"/rollsRoyce.obj",
      aXYZ,aTri);
  dfm2::Normalize_Points3(aXYZ,4.0);
  // ---------------------------------------
  
  dfm2::opengl::CDrawerOldGL_Render2Tex draw_sampler;
  dfm2::opengl::CRender2Tex sampler;
  {
    unsigned int nresX = 128;
    unsigned int nresY = 128;
    unsigned int nresZ = 256;
    double elen = 0.02;
    dfm2::CVec3d origin = dfm2::CVec3d(-0.5 * elen * nresX, +0.5 * elen * nresY, +0.5 * elen * nresZ);
    dfm2::CVec3d ez = dfm2::CVec3d(0, +1, 0);
    dfm2::CVec3d ex = dfm2::CVec3d(1, +0, 0);
    //
    sampler.SetTextureProperty(nresX, nresZ, true);
    ::delfem2::Mat4_OrthongoalProjection_AffineTrans(
        sampler.mMV, sampler.mP,
        origin.p, ez.p, ex.p,
        nresX, nresZ, elen, elen * nresY);
    draw_sampler.SetPointColor(0.0, 1.0, 0.0);
    draw_sampler.draw_len_axis = 0.2;
    draw_sampler.isDrawTex = false;
    draw_sampler.isDrawOnlyHitPoints = true;
  }
  // ---------------------------------------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.camera.view_height = 2.0;
  viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  viewer.camera.Rot_Camera(+0.2, -0.2);
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }

  dfm2::opengl::setSomeLighting();
  ::glEnable(GL_DEPTH_TEST);
  
  {
    sampler.InitGL(); // move the sampled image to a texture
    sampler.Start();
    dfm2::opengl::SetView(sampler);
    ::glClearColor(1.0, 1.0, 1.0, 1.0 );
    ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    ::glEnable(GL_DEPTH_TEST);
    ::glDisable(GL_BLEND);
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
    sampler.End();
    sampler.CopyToCPU_Depth();
    sampler.CopyToCPU_RGBA8UI();
  }

  for(int iframe=0;;iframe++){
    dfm2::CVec3d dir(cos(iframe*0.003),sin(iframe*0.005),sin(iframe*0.007));
    dfm2::CVec3d src(1.5*cos(iframe*0.005), 1.5*sin(iframe*0.007), 1.5*sin(iframe*0.009) );
    // ----

    std::vector<double> aXYZ1;
    {
      dfm2::CMat4d mMVPG; sampler.GetMVPG(mMVPG.mat);
      const dfm2::CMat4d mMVPGinv = mMVPG.Inverse();
      dfm2::CVec3d p0; dfm2::Vec3_Vec3Mat4_AffineProjection(p0.p,src.p,mMVPG.mat);
      dfm2::CVec3d p1; dfm2::Vec3Mat4(p1.p,dir.p,mMVPG.mat);
      std::vector<dfm2::CPtElm2<double>> aPES;
      dfm2::IntersectionLine_Hightfield(aPES,
          p0.p, p1.Normalize().p,
          sampler.nResX, sampler.nResY,
          sampler.aZ);
      for(const auto & pes : aPES){
        dfm2::CVec3d lpos = pes.Pos_Grid(sampler.nResX, sampler.nResY, 1.0, sampler.aZ);
        dfm2::CVec3d q2; dfm2::Vec3_Vec3Mat4_AffineProjection(q2.p, lpos.p, mMVPGinv.mat);
        aXYZ1.push_back(q2.x());
        aXYZ1.push_back(q2.y());
        aXYZ1.push_back(q2.z());
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
//        ::glVertex3dv((src-dir.Normalize()).p);
//        ::glColor3d(0,1,0);
//        ::glVertex3dv((src+dir.Normalize()).p);
        ::glVertex3dv(src.p);
        ::glEnd();
      }
      {
        ::glColor3d(1,0,0);
        dfm2::CVec3d p0 = src + 10.0*dir;
        dfm2::CVec3d p1 = src - 10.0*dir;
        ::glLineWidth(1);
        ::glBegin(GL_LINES);
        ::glVertex3d(p0.x(), p0.y(), p0.z());
        ::glVertex3d(p1.x(), p1.y(), p1.z());
        ::glEnd();
      }
      for(unsigned int ixyz=0;ixyz<aXYZ1.size()/3;++ixyz){
        ::glBegin(GL_POINTS);
        ::glColor3d(0,0,1);
        ::glVertex3d(aXYZ1[ixyz*3+0], aXYZ1[ixyz*3+1], aXYZ1[ixyz*3+2]);
        ::glEnd();
      }
      // ----------
      ::glEnable(GL_LIGHTING);
      dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
      glPointSize(1);
      draw_sampler.Draw(sampler);
      viewer.SwapBuffers();
      glfwPollEvents();
    }
    if( glfwWindowShouldClose(viewer.window) ){ break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


