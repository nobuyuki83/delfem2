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

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_color.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/render2tex_glold.h"
//
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

std::vector<double> aXYZ;
std::vector<unsigned int> aTri;

// ------------------------------------------------------

void DrawObject(){
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
}

void myGlutDisplay(const std::vector<dfm2::opengl::CRender2Tex_DrawOldGL>& aSampler)
{
  dfm2::opengl::DrawBackground( dfm2::CColor(0.2,0.7,0.7) );
  ::glEnable(GL_LIGHTING);
  ::glColor3d(1,1,1);
  DrawObject();

  glPointSize(3);
  float mMV[16]; glGetFloatv(GL_MODELVIEW, mMV);
  float mP[16]; glGetFloatv(GL_PROJECTION, mP);
  for(auto& smplr: aSampler){
    smplr.Draw();
  }
}

int main(int argc,char* argv[])
{
  dfm2::Read_Obj(std::string(PATH_INPUT_DIR)+"/rollsRoyce.obj",
    aXYZ,aTri);
  dfm2::Normalize_Points3(aXYZ,4.0);
  // ---------------------------------------
  
  std::vector<dfm2::opengl::CRender2Tex_DrawOldGL> aSampler;
  {
    unsigned int nresX = 128;
    unsigned int nresY = 128;
    unsigned int nresZ = 256;
    double elen = 0.02;
    
    aSampler.resize(6);
    aSampler[0].SetTextureProperty(nresY, nresZ, true);
    aSampler[0].SetCoord(elen, elen*nresX,
                         dfm2::CVec3d(+0.5*elen*nresX,-0.5*elen*nresY,-0.5*elen*nresZ).stlvec(),
                         dfm2::CVec3d(+1,  0, 0).stlvec(),
                         dfm2::CVec3d( 0, +1, 0).stlvec() );
    aSampler[0].SetPointColor(1.0, 0.0, 0.0);
    //
    aSampler[1].SetTextureProperty(nresY, nresZ, true);
    aSampler[1].SetCoord(elen, elen*nresX,
                         dfm2::CVec3d(-0.5*elen*nresX,-0.5*elen*nresY,+0.5*elen*nresZ).stlvec(),
                         dfm2::CVec3d(-1,  0, 0).stlvec(),
                         dfm2::CVec3d( 0, +1, 0).stlvec() );
    aSampler[1].SetPointColor(1.0, 0.5, 0.5);
    //
    aSampler[2].SetTextureProperty(nresX, nresZ, true);
    aSampler[2].SetCoord(elen, elen*nresY,
                         dfm2::CVec3d(-0.5*elen*nresX,+0.5*elen*nresY,+0.5*elen*nresZ).stlvec(),
                         dfm2::CVec3d(0,+1,0).stlvec(),
                         dfm2::CVec3d(1,+0,0).stlvec() );
    aSampler[2].SetPointColor(0.0, 1.0, 0.0);
    //
    aSampler[3].SetTextureProperty(nresX, nresZ, true);
    aSampler[3].SetCoord(elen, elen*nresY,
                         dfm2::CVec3d(-0.5*elen*nresX,-0.5*elen*nresY,-0.5*elen*nresZ).stlvec(),
                         dfm2::CVec3d(0,-1,0).stlvec(),
                         dfm2::CVec3d(1,+0,0).stlvec() );
    aSampler[3].SetPointColor(0.5, 1.0, 0.5);
    //
    aSampler[4].SetTextureProperty(nresX, nresY, true);
    aSampler[4].SetCoord(elen, elen*nresZ,
                         dfm2::CVec3d(-0.5*elen*nresX,-0.5*elen*nresY,+0.5*elen*nresZ).stlvec(),
                         dfm2::CVec3d(0,0,+1).stlvec(),
                         dfm2::CVec3d(1,0,0).stlvec() );
    aSampler[4].SetPointColor(0.0, 0.0, 1.0);
    //
    aSampler[5].SetTextureProperty(nresX, nresY, true);
    aSampler[5].SetCoord(elen, elen*nresZ,
                         dfm2::CVec3d(-0.5*elen*nresX,+0.5*elen*nresY,-0.5*elen*nresZ).stlvec(),
                         dfm2::CVec3d(0,0,-1).stlvec(),
                         dfm2::CVec3d(1,0,0).stlvec() );
    aSampler[5].SetPointColor(0.5, 0.5, 1.0);
  }
  for(auto& smplr : aSampler){
    smplr.draw_len_axis = 0.2;
    smplr.isDrawTex = false;
    smplr.isDrawOnlyHitPoints = true;
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
  
  for(auto& smplr: aSampler){
    smplr.InitGL(); // move the sampled image to a texture
    smplr.Start();
    ::glClearColor(1.0, 1.0, 1.0, 1.0 );
    ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    ::glEnable(GL_DEPTH_TEST);
    ::glDisable(GL_BLEND);
    ::glEnable(GL_LIGHTING);
    DrawObject();
    smplr.End();
    smplr.GetDepth();
    smplr.GetColor();
  }

  while (!glfwWindowShouldClose(viewer.window))
  {
    // ----
    viewer.DrawBegin_oldGL();
    myGlutDisplay(aSampler);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


