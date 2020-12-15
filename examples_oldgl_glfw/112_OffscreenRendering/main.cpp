/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <glad/glad.h>
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/color_glold.h"
#include "delfem2/opengl/v3q_glold.h"
#include "delfem2/opengl/r2tglo_glold.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include <GLFW/glfw3.h>

namespace dfm2 = delfem2;

// ------------------------------------------------------

void DrawObject(
    double cur_time,
    std::vector<double>& aXYZ,
    std::vector<unsigned int>& aTri)
{
  ::glRotated(+cur_time, 1,0,0);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
  ::glRotated(-cur_time, 1,0,0);
}

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  dfm2::Read_Obj(
      std::string(PATH_INPUT_DIR)+"/bunny_1k.obj",
    aXYZ,aTri);
  dfm2::Normalize_Points3(
      aXYZ,
      1.0);
  // ---------------------------------------
  int nres = 100;
  double elen = 0.02;
  dfm2::opengl::CRender2Tex_DrawOldGL sampler;
  sampler.SetTextureProperty(nres, nres, true);

/*
  dfm2::Mat4_OrthongoalProjection_AffineTrans(
      sampler.mMV, sampler.mP,
      dfm2::CVec3d(-nres*elen*0.5,nres*elen*0.5,-2).p,
      dfm2::CVec3d(0,0,-1).p,
      dfm2::CVec3d(1,0,0).p,
      nres, nres, elen, 4);
      */

  ::delfem2::Mat4_OrthongoalProjection_AffineTrans(
      sampler.mMV, sampler.mP,
      dfm2::CVec3d(+0.5 * elen * nres, -0.5 * elen * nres, -0.5 * elen * nres).p,
      dfm2::CVec3d(+1, 0, 0).p,
      dfm2::CVec3d(0, +1, 0).p,
      nres, nres, elen, 2);

  sampler.SetPointColor(1, 0, 0);
  sampler.draw_len_axis = 1.0;
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
  
  sampler.InitGL(); // move the sampled image to a texture

  double cur_time = 0.0;
  while (!glfwWindowShouldClose(viewer.window))
  {
    sampler.Start();
    ::glClearColor(1.0, 1.0, 1.0, 1.0 );
    ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    ::glEnable(GL_DEPTH_TEST);
    ::glDisable(GL_BLEND);
    ::glEnable(GL_LIGHTING);
    DrawObject(cur_time,aXYZ,aTri);
    sampler.End();
    sampler.GetDepth();
    sampler.GetColor();
    cur_time += 1.0;
    // ----
    viewer.DrawBegin_oldGL();
    dfm2::opengl::DrawBackground( dfm2::CColor(0.2f,0.7f,0.7f) );
    ::glEnable(GL_LIGHTING);
    ::glColor3d(1,1,1);
    DrawObject(cur_time,aXYZ,aTri);
    sampler.Draw();
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();    
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


