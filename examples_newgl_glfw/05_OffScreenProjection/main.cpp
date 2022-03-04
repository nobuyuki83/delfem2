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

#include "delfem2/msh_primitive.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/opengl/new/drawer_mshtri.h"
#include "delfem2/opengl/new/drawer_render2tex.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

struct MyData {
  dfm2::opengl::CRender2Tex r2t;
  dfm2::opengl::CShader_TriMesh shdr_trimsh;
  dfm2::opengl::CRender2Tex_DrawNewGL drawer_r2t;
  delfem2::glfw::CViewer3 viewer;
};

void draw(MyData* data) {
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);
  data->shdr_trimsh.Draw(
      data->viewer.GetProjectionMatrix().data(),
      data->viewer.GetModelViewMatrix().data());
  data->drawer_r2t.Draw(
      data->r2t,
      data->viewer.GetProjectionMatrix().data(),
      data->viewer.GetModelViewMatrix().data());
  data->viewer.SwapBuffers();
  glfwPollEvents();
}

int main() {
  MyData data;
  {
    data.viewer.projection = std::make_unique<delfem2::Projection_LookOriginFromZplus>(2, false);
    int nres = 200;
    data.r2t.SetTextureProperty(nres, nres, true);
    dfm2::Mat4_Identity(data.r2t.mat_modelview);
    dfm2::Mat4_Identity(data.r2t.mat_projection);
    data.drawer_r2t.draw_len_axis = 1.0;
  }

  dfm2::glfw::InitGLNew();
  data.viewer.OpenWindow();
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cerr << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif

  data.r2t.InitGL();
  data.drawer_r2t.InitGL();
  data.shdr_trimsh.InitGL();

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
    data.shdr_trimsh.Initialize(
        vtx_xyz, 3, tri_vtx);
  }
  data.r2t.Start();
  ::glClearColor(1.0, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glDisable(GL_BLEND);
  ::glEnable(GL_DEPTH_TEST);
  data.shdr_trimsh.Draw(
      dfm2::CMat4f(data.r2t.mat_modelview).data(),
      dfm2::CMat4f(data.r2t.mat_projection).data());
  data.r2t.End();
  data.drawer_r2t.SetDepth(data.r2t);

#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, &data, 0, 1);
#else
  while (!glfwWindowShouldClose(data.viewer.window)) { draw(&data); }
#endif

  glfwDestroyWindow(data.viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

