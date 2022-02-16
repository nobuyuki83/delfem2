/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should be before glfw3.h
#endif
#include <glad/glad.h> // glad need to be defiend in the begenning
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/r2tglo.h"
#include "delfem2/opengl/funcs.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/mshmisc.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/file.h"

namespace dfm2 = delfem2;

// ------------------------------------------------------

void DrawObject(
    double cur_time,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri) {
  ::glRotated(+cur_time, 1, 0, 0);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aTri);
  ::glRotated(-cur_time, 1, 0, 0);
}

int main() {
  std::vector<double> vtx_xyz;
  std::vector<unsigned int> tri_vtx;
  dfm2::Read_Obj(
      vtx_xyz, tri_vtx,
      std::filesystem::path(PATH_ASSET_DIR) / "bunny_1k.obj");
  dfm2::Normalize_Points3(
      vtx_xyz,
      1.0);
  // ---------------------------------------
  int nres = 100;
  double elen = 0.02;
  dfm2::opengl::CRender2Tex smpl;
  smpl.SetTextureProperty(nres, nres, true);

  dfm2::CMat4d::AffineAxisTransform(
      {0,0,1},
      {1,0,0},
      {0,1,0}
      ).CopyTo(smpl.mat_modelview);
  dfm2::CMat4d::AffineOrthogonalProjection(
      -elen*nres*0.5, elen*nres*0.5,
      -elen*nres*0.5, elen*nres*0.5,
      -1, 1
      ).CopyTo(smpl.mat_projection);

  dfm2::opengl::CDrawerOldGL_Render2Tex drawer_r2t;
  drawer_r2t.SetPointColor(1, 0, 0);
  drawer_r2t.draw_len_axis = 1.0;
  // ---------------------------------------
  dfm2::glfw::CViewer3 viewer(2.0);
  dfm2::glfw::InitGLOld();
  viewer.OpenWindow();
  if (!gladLoadGL()) {     // glad: load all OpenGL function pointers
    printf("Something went wrong in loading OpenGL functions!\n");
    exit(-1);
  }

  int shaderProgram;
  {
    std::string vrt_path = std::string(PATH_SOURCE_DIR) + "/glsl120_normalmap.vert";
    std::string frg_path = std::string(PATH_SOURCE_DIR) + "/glsl120_normalmap.frag";
    std::string vrt = dfm2::LoadFile(vrt_path);
    std::string frg = dfm2::LoadFile(frg_path);
    shaderProgram = dfm2::opengl::setUpGLSL(vrt, frg);
  }

  dfm2::opengl::setSomeLighting();
  ::glEnable(GL_DEPTH_TEST);

  smpl.InitGL(); // move the sampled image to a texture

  double cur_time = 0.0;
  while (!glfwWindowShouldClose(viewer.window)) {
    smpl.Start();
    SetView(smpl);
    ::glClearColor(1.0, 1.0, 1.0, 1.0);
    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    ::glEnable(GL_DEPTH_TEST);
    ::glDisable(GL_BLEND);
    ::glEnable(GL_LIGHTING);
    ::glUseProgram(shaderProgram);
    DrawObject(cur_time, vtx_xyz, tri_vtx);
    smpl.End();
    cur_time += 1.0;
    // ----
    viewer.DrawBegin_oldGL();
    ::glUseProgram(0);
    dfm2::opengl::DrawBackground(dfm2::CColor(0.2f, 0.7f, 0.7f));
    ::glEnable(GL_LIGHTING);
    ::glColor3d(1, 1, 1);
    ::glUseProgram(shaderProgram);
    DrawObject(cur_time, vtx_xyz, tri_vtx);
    ::glUseProgram(0);
    drawer_r2t.Draw(smpl);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


