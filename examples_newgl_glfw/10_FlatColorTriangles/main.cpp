

#include <cstdio>
#include <iostream>
#include <vector>
#if defined(_MSC_VER)
#  include <windows.h>
#endif
//
#define GL_SILENCE_DEPRECATION
#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#  define GL_GLEXT_PROTOTYPES
#  define EGL_EGLEXT_PROTOTYPES
#else
#  include <glad/glad.h>
#endif
#include <GLFW/glfw3.h>

#include "delfem2/vec3.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_unindexed.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/opengl/new/drawer_mshunindex.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

// ---------------------------------------------------------------

delfem2::glfw::CViewer3 viewer;
dfm2::opengl::Drawer_MeshUnIndexed drawer;
std::vector<double> vtx_xyz;
std::vector<unsigned int> tri_vtx;
std::vector<int> tri_flg;
std::vector<double> tri_rgb;
const double flg_rgb[2][3] = {
  {1, 1, 1},
  {1, 0, 0}};

void draw(GLFWwindow *window) {
  const double time0 = glfwGetTime();
  { // make flag
    for (unsigned int it = 0; it < tri_vtx.size() / 3; ++it) {
      const double *p0 = vtx_xyz.data() + tri_vtx[it * 3 + 0] * 3;
      const double *p1 = vtx_xyz.data() + tri_vtx[it * 3 + 1] * 3;
      const double *p2 = vtx_xyz.data() + tri_vtx[it * 3 + 2] * 3;
      double y0 = (p0[0] + p1[0] + p2[0]) / 3.0;
      if (y0 > 0.3*std::sin(time0) ) { tri_flg[it] = 1; }
      else{ tri_flg[it] = 0; }
    }
    dfm2::UnindexedColorTriMesh3(
      tri_rgb,
      flg_rgb, tri_flg, tri_vtx);
    drawer.UpdateTriRgb(tri_rgb);
  }

  // ==========

  glfwPollEvents();
  ::glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);
  {
    ::glEnable(GL_DEPTH_TEST);
    drawer.Draw(
      viewer.GetProjectionMatrix().data(),
      viewer.GetModelViewMatrix().data());
  }

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);

  glfwSwapBuffers(window);
}

int main(int, char **) {

  delfem2::glfw::InitGLNew();
  viewer.OpenWindow();
#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif
  drawer.InitGL();
  {
    std::string filePathName = std::string(PATH_SOURCE_DIR) + "/../../test_inputs/bunny_1k.obj";
    delfem2::Read_Obj3(
      vtx_xyz, tri_vtx,
      filePathName);
    delfem2::Normalize_Points3(
      vtx_xyz,
      1.);
    tri_flg.assign(tri_vtx.size() / 3, 0);
    std::vector<double> tri_xyz, tri_nrm;
    dfm2::UnidexedVertexDataTriMesh(
      tri_xyz,
      vtx_xyz, tri_vtx);
    dfm2::UnindexedNormalTriMesh3(
      tri_nrm,
      vtx_xyz, tri_vtx);
    dfm2::UnindexedColorTriMesh3(
      tri_rgb,
      flg_rgb, tri_flg, tri_vtx);
    drawer.Initialize(tri_xyz, tri_nrm, tri_rgb, 3);
  }

#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  return 0;
}
