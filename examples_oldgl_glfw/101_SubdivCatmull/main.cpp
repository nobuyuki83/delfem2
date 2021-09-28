/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @brief simple demo of subdivision surface
 */

#include <cstdlib>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/mshmisc.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/mshsubdiv.h"
#include "delfem2/mshprimitive.h"

namespace dfm2 = delfem2;

// --------------------------------------------------

int main() {
  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 2.0;
  delfem2::opengl::setSomeLighting();

  for (unsigned itr = 0;; ++itr) {
    std::vector<std::vector<unsigned int> > array_quad_vtx;
    std::vector<std::vector<double> > array_vtx_xyz;
    unsigned int nlevel_subdiv;
    array_quad_vtx.resize(1);
    array_vtx_xyz.resize(1);
    if (itr % 2 == 0) {
      const double bbmin[3] = {-1, -1, -1};
      const double bbmax[3] = {+1, +1, +1};
      delfem2::MeshQuad3_CubeVox(
          array_vtx_xyz[0], array_quad_vtx[0],
          bbmin, bbmax);
      nlevel_subdiv = 5;
    } else {
      dfm2::Read_Obj_MeshQuad3(
          array_vtx_xyz[0], array_quad_vtx[0],
          std::filesystem::path(PATH_INPUT_DIR) / "basemesh_hand.obj");
      nlevel_subdiv = 3;
    }
    array_vtx_xyz.resize(nlevel_subdiv + 1);
    array_quad_vtx.resize(nlevel_subdiv + 1);
    for (unsigned int il = 0; il < nlevel_subdiv; ++il) {
      const std::vector<double> &vtx_xyz0 = array_vtx_xyz[il];
      const std::vector<unsigned int> &quads0 = array_quad_vtx[il];
      std::vector<unsigned int> &quads1 = array_quad_vtx[il + 1];
      std::vector<unsigned int> aEdgeFace0;
      std::vector<unsigned int> psupIndQuad0, psupQuad0;
      dfm2::SubdivTopo_MeshQuad(
          quads1,
          psupIndQuad0, psupQuad0, aEdgeFace0,
          quads0.data(), quads0.size() / 4,
          vtx_xyz0.size() / 3);
      std::vector<double> &vtx_xyz1 = array_vtx_xyz[il + 1];
      delfem2::SubdivisionPoints_QuadCatmullClark(
          vtx_xyz1,
          quads1, aEdgeFace0, psupIndQuad0, psupQuad0,
          quads0.data(), quads0.size() / 4,
          vtx_xyz0.data(), vtx_xyz0.size() / 3);
    }
    // -------
    for (unsigned int isub = 0; isub <= nlevel_subdiv; ++isub) {
      for (unsigned int iframe = 0; iframe < 30; ++iframe) {
        viewer.DrawBegin_oldGL();
        ::glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
        ::glEnable(GL_LIGHTING);
        delfem2::opengl::DrawMeshQuad3D_FaceNorm(array_vtx_xyz[isub], array_quad_vtx[isub]);
        ::glDisable(GL_LIGHTING);
        ::glColor3d(0, 0, 0);
        delfem2::opengl::DrawMeshQuad3D_Edge(array_vtx_xyz[isub], array_quad_vtx[isub]);
        viewer.SwapBuffers();
        glfwPollEvents();
        if (glfwWindowShouldClose(viewer.window)) break;
      }
      if (glfwWindowShouldClose(viewer.window)) break;
    }
    if (glfwWindowShouldClose(viewer.window)) break;
  }
  // ----------------------
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
