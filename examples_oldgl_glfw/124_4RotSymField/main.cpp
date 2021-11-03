/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief implementation of 4 rotatoinal symetry field
 * @details implementation is based on "Wenzel Jakob, Marco Tarini, Daniele Panozzo, and Olga Sorkine-Hornung. Instant field-aligned meshes. Siggraph Asia 2015"
 */

#include <cstdlib>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/msh_io_ply.h"
#include "delfem2/mshmisc.h"
#include "delfem2/points.h"
#include "delfem2/mshuni.h"
#include "delfem2/vec3.h"
#include "delfem2/4rotsym.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/v3q.h"

namespace dfm2 = delfem2;

// ------------------------------------------------
// TODO: Add random permutation version in the demo

int main() {
  std::vector<double> vtx_xyz;
  std::vector<unsigned int> tri_vtx;
  delfem2::Read_Ply(
      vtx_xyz, tri_vtx,
      std::filesystem::path(PATH_INPUT_DIR) / "bunny_1k.ply");
  delfem2::Normalize_Points3(vtx_xyz);
  std::vector<double> vtx_normal(vtx_xyz.size());
  dfm2::Normal_MeshTri3D(
      vtx_normal.data(),
      vtx_xyz.data(), vtx_xyz.size() / 3,
      tri_vtx.data(), tri_vtx.size() / 3);

  std::vector<double> aOdir;
  {
    const double minCoords[3] = {-1., -1., -1.};
    const double maxCoords[3] = {+1., +1., +1.};
    aOdir.resize(vtx_xyz.size());
    dfm2::Points_RandomUniform(aOdir.data(),
                               vtx_xyz.size() / 3, 3,
                               minCoords, maxCoords);
    dfm2::TangentVector_Points3(aOdir,
                                vtx_normal);
  }

  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind, psup,
      tri_vtx.data(), tri_vtx.size() / 3, 3,
      vtx_xyz.size() / 3);

  // ------------------
  delfem2::glfw::CViewer3 viewer;
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  dfm2::opengl::setSomeLighting();
  unsigned int iframe = 0;
  while (true) {
    if (iframe == 0) {
      const double minCoords[3] = {-1., -1., -1.};
      const double maxCoords[3] = {+1., +1., +1.};
      aOdir.resize(vtx_xyz.size());
      dfm2::Points_RandomUniform(
          aOdir.data(),
          vtx_xyz.size() / 3, 3,
          minCoords, maxCoords);
      dfm2::TangentVector_Points3(aOdir, vtx_normal);
    }
    if (iframe > 30) {
      dfm2::Smooth4RotSym(aOdir,
                          vtx_normal, psup_ind, psup);
    }
    // --------------------
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(
        vtx_xyz.data(),
        tri_vtx.data(), tri_vtx.size() / 3);
    {
      ::glDisable(GL_LIGHTING);
      double len = 0.03;
      ::glLineWidth(3);
      const size_t np = vtx_xyz.size() / 3;
      for (unsigned int ip = 0; ip < np; ++ip) {
        const dfm2::CVec3d p = dfm2::CVec3d(vtx_xyz.data() + ip * 3);
        const dfm2::CVec3d n = dfm2::CVec3d(vtx_normal.data() + ip * 3).normalized();
        const dfm2::CVec3d o = dfm2::CVec3d(aOdir.data() + ip * 3).normalized();
        const dfm2::CVec3d q = dfm2::Cross(n, o);
        ::glBegin(GL_LINES);
        ::glColor3d(0, 0, 0);
        dfm2::opengl::myGlVertex(p);
        dfm2::opengl::myGlVertex(p + len * n);
        ::glColor3d(0, 0, 1);
        dfm2::opengl::myGlVertex(p - len * o);
        dfm2::opengl::myGlVertex(p);
        ::glColor3d(1, 0, 0);
        dfm2::opengl::myGlVertex(p);
        dfm2::opengl::myGlVertex(p + len * o);
        dfm2::opengl::myGlVertex(p - len * q);
        dfm2::opengl::myGlVertex(p + len * q);
        ::glEnd();
      }
    }
    iframe = (iframe + 1) % 100;
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if (glfwWindowShouldClose(viewer.window)) { goto EXIT; }
  }
  EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
