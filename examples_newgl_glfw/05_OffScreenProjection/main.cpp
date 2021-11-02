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
#include "delfem2/points.h"
#include "delfem2/opengl/new/shdr_points.h"
#include "delfem2/opengl/new/shdr_mshtri.h"
#include "delfem2/opengl/new/r2tgln.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

delfem2::glfw::CViewer3 viewer(2);
dfm2::opengl::CRender2Tex r2t;
dfm2::opengl::CShader_TriMesh shdr_trimsh;
dfm2::opengl::CRender2Tex_DrawNewGL drawer_r2t;

void draw() {
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);
  shdr_trimsh.Draw(
      viewer.GetProjectionMatrix().data(),
      viewer.GetModelViewMatrix().data());
  drawer_r2t.Draw(
      r2t,
      viewer.GetProjectionMatrix().data(),
      viewer.GetModelViewMatrix().data());
  viewer.SwapBuffers();
  glfwPollEvents();
}

int main() {
  {
    int nres = 200;
    r2t.SetTextureProperty(nres, nres, true);
    dfm2::Mat4_Identity(r2t.mat_modelview);
    dfm2::Mat4_Identity(r2t.mat_projection);
    drawer_r2t.draw_len_axis = 1.0;
  }

  dfm2::glfw::InitGLNew();
  viewer.InitGL();

#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cerr << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif

  r2t.InitGL();
  drawer_r2t.InitGL();
  shdr_trimsh.Compile();

  {
    std::vector<float> vtx_xyz;
    std::vector<unsigned int> tri_vtx;
    dfm2::MeshTri3_Torus(
        vtx_xyz, tri_vtx,
        0.8f, 0.1f,
        8, 8);
    dfm2::Rotate_Points3(
        vtx_xyz,
        0.1f, 0.2f, 0.3f);
    shdr_trimsh.Initialize(
        vtx_xyz, 3, tri_vtx);
  }
  r2t.Start();
  ::glClearColor(1.0, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glDisable(GL_BLEND);
  ::glEnable(GL_DEPTH_TEST);
  shdr_trimsh.Draw(
      dfm2::CMat4f(r2t.mat_modelview).data(),
      dfm2::CMat4f(r2t.mat_projection).data());
  r2t.End();
  drawer_r2t.SetDepth(r2t);

#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(); }
#endif

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

