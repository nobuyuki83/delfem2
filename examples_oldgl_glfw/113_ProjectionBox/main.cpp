/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <glad/glad.h> // glad need to be defiend in the begenning
#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/r2tglo_glold.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include <GLFW/glfw3.h>

namespace dfm2 = delfem2;

// ------------------------------------------------------

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  dfm2::Read_Obj(
      std::string(PATH_INPUT_DIR)+"/rollsRoyce.obj",
    aXYZ,aTri);
  dfm2::Normalize_Points3(aXYZ,4.0);
  // ---------------------------------------

  dfm2::opengl::CRender2Tex_DrawOldGL_BOX sampler_box;
  sampler_box.Initialize(128, 128, 256, 0.02);

  for(auto& smplr : sampler_box.aSampler){
    smplr.draw_len_axis = 0.2;
    smplr.isDrawTex = false;
    smplr.isDrawOnlyHitPoints = true;
  }
  // ---------------------------------------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.camera.view_height = 2.0;
  viewer.camera.camera_rot_mode = dfm2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
//  viewer.camera.Rot_Camera(+0.2, -0.2);
  if(!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }

  dfm2::opengl::setSomeLighting();
  ::glEnable(GL_DEPTH_TEST);


  sampler_box.Draw();
  for(auto& smplr: sampler_box.aSampler){
    smplr.InitGL(); // move the sampled image to a texture
    smplr.Start();
    ::glClearColor(1.0, 1.0, 1.0, 1.0 );
    ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    ::glEnable(GL_DEPTH_TEST);
    ::glDisable(GL_BLEND);
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
    smplr.End();
    smplr.GetDepth();
    smplr.GetColor();
  }

  while (!glfwWindowShouldClose(viewer.window))
  {
    viewer.DrawBegin_oldGL();
    sampler_box.Draw();
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


