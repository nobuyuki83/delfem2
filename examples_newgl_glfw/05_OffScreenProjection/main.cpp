/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#else
#  include <glad/glad.h>
#endif
#if defined(_MSC_VER)
#  include <windows.h>
#endif
#include <GLFW/glfw3.h>

#include "delfem2/mshprimitive.h"
#include "delfem2/vec3.h"
#include "delfem2/opengl/new/shdr_points.h"
#include "delfem2/opengl/new/shdr_mshtri.h"
#include "delfem2/opengl/new/r2tgln.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

delfem2::glfw::CViewer3 viewer(2);
dfm2::opengl::CRender2Tex r2t;
dfm2::opengl::CShader_TriMesh shdr_trimsh;
dfm2::opengl::CShader_Points shdr_points;
dfm2::opengl::CRender2Tex_DrawNewGL draw_r2t;

void draw() {
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);

  dfm2::CMat4f mP = viewer.GetProjectionMatrix();
  mP = mP.transpose() * dfm2::CMat4f::ScaleXYZ(1,1,-1);
  dfm2::CMat4f mMV = viewer.GetModelViewMatrix();
  mMV = mMV.transpose();
  shdr_points.Draw(GL_POINTS, mP.data(), mMV.data());
  shdr_trimsh.Draw(mP.data(), mMV.data());
  draw_r2t.Draw(r2t, mP.data(), mMV.data());
  viewer.SwapBuffers();
  glfwPollEvents();
}

int main() {
  {
    int nres = 200;
    double elen = 0.01;
    r2t.SetTextureProperty(nres, nres, true);
    dfm2::Mat4_OrthongoalProjection_AffineTrans(
        r2t.mat_modelview_colmajor, r2t.mat_projection_colmajor,
        dfm2::CVec3d(-nres * elen * 0.5, nres * elen * 0.5, -2).p,
        dfm2::CVec3d(0, 0, -1).p,
        dfm2::CVec3d(1, 0, 0).p,
        nres, nres, elen, 4.0);
    draw_r2t.draw_len_axis = 1.0;
  }

  dfm2::glfw::InitGLNew();
  viewer.InitGL();

#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif

  r2t.InitGL();
  draw_r2t.InitGL();

  {
    std::vector<double> aXYZ;
    std::vector<unsigned int> aTri;
    dfm2::MeshTri3_Torus(aXYZ, aTri, 0.8, 0.1, 8, 8);
    shdr_trimsh.Compile();
    shdr_trimsh.Initialize(aXYZ, 3, aTri);
  }
  r2t.Start();
  ::glClearColor(1.0, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glDisable(GL_BLEND);
  ::glEnable(GL_DEPTH_TEST);
  {
    float mMVf[16];
    dfm2::Copy_Mat4(mMVf, r2t.mat_modelview_colmajor);
    float mPf[16];
    dfm2::Copy_Mat4(mPf, r2t.mat_projection_colmajor);
    shdr_trimsh.Draw(mPf, mMVf);
  }
  r2t.End();
  draw_r2t.SetDepth(r2t);

#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(); }
#endif

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

