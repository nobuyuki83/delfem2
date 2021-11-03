/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <cmath>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // should be before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/gizmo_geo3.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/points.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/gizmo.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/v3q.h"

namespace dfm2 = delfem2;

// -------------------------------------------

int main() {
  class CMyViewer : public delfem2::glfw::CViewer3 {
   public:
    CMyViewer() {
      delfem2::Read_Ply(
          aXYZ, aTri,
          std::filesystem::path(PATH_INPUT_DIR) / "bunny_1k.ply");
      delfem2::Normalize_Points3(aXYZ);
    }
    //
    void mouse_press(const float src[3], const float dir[3]) override {
      ga.Pick(src, dir);
    }
    void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) override {
      ga.Drag(src0, src1, dir);
    }
    void key_release([[maybe_unused]] int key, [[maybe_unused]] int mods) override {
    }
    void key_press(int key, int mods) override {
      delfem2::glfw::CViewer3::key_press(key, mods);
      if (key == GLFW_KEY_R) { ga.igizmo_mode = 1; }
      if (key == GLFW_KEY_G) { ga.igizmo_mode = 0; }
    }
    //
    void Draw() {
      DrawBegin_oldGL();
      delfem2::opengl::DrawAxis(1);
      {
        ::glMatrixMode(GL_MODELVIEW);
        ::glPushMatrix();
        const auto m0 = ga.Affine();
        const auto m1 = m0.transpose();
        delfem2::opengl::MyGlMultMat(m1);
        // ------
        ::glEnable(GL_LIGHTING);
        ::glColor3d(0, 0, 0);
        delfem2::opengl::DrawMeshTri3D_Edge(
            aXYZ.data(), aXYZ.size() / 3,
            aTri.data(), aTri.size() / 3);
        delfem2::opengl::DrawMeshTri3D_FaceNorm(
            aXYZ.data(),
            aTri.data(), aTri.size() / 3);
        // -------
        ::glMatrixMode(GL_MODELVIEW);
        ::glPopMatrix();
      }
      delfem2::opengl::Draw(ga);
      SwapBuffers();
    }
   public:
    delfem2::CGizmo_Affine<float> ga;
    std::vector<double> aXYZ;
    std::vector<unsigned int> aTri;
  } viewer;
  // ------------------------
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  // --------------------
  while (!glfwWindowShouldClose(viewer.window)) {
    viewer.Draw();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


