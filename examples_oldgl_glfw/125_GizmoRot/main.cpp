/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/gizmo_geo3.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/mat4.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/gizmo.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"

namespace dfm2 = delfem2;

// -------------------------------------------

int main() {
  class CMyViewer : public delfem2::glfw::CViewer3 {
   public:
    CMyViewer() {
      delfem2::Read_Ply(
          vtx_xyz, tri_vtx,
          std::filesystem::path(PATH_INPUT_DIR) / "bunny_1k.ply");
      delfem2::Normalize_Points3(vtx_xyz);
      gizmo_rot.size = 0.7f;
    }
    //
    void mouse_press(const float src[3], const float dir[3]) override {
      gizmo_rot.Pick(true, src, dir, 0.1f);
    }
    void mouse_drag(const float src0[3], const float src1[3], const float dir[3]) override {
      gizmo_rot.Drag(src0, src1, dir);
    }
    void mouse_release() override {
      gizmo_rot.ielem_picked = -1;
    }
    //
    void Draw() {
      DrawBegin_oldGL();
      {
        float r[16];
        dfm2::Mat4_AffineQuaternion(r, gizmo_rot.quat);
        float r0[16];
        dfm2::Transpose_Mat4(r0, r);
        ::glMatrixMode(GL_MODELVIEW);
        ::glPushMatrix();
        ::glMultMatrixf(r0);
        ::glEnable(GL_LIGHTING);
        ::glColor3d(0, 0, 0);
        delfem2::opengl::DrawMeshTri3D_Edge(
            vtx_xyz.data(), vtx_xyz.size() / 3,
            tri_vtx.data(), tri_vtx.size() / 3);
        delfem2::opengl::DrawMeshTri3D_FaceNorm(
            vtx_xyz.data(),
            tri_vtx.data(), tri_vtx.size() / 3);
        ::glMatrixMode(GL_MODELVIEW);
        ::glPopMatrix();
      }
      delfem2::opengl::Draw(gizmo_rot);
      SwapBuffers();
    }
   public:
    dfm2::CGizmo_Rotation<float> gizmo_rot;
    std::vector<double> vtx_xyz;
    std::vector<unsigned int> tri_vtx;
  } viewer;
  // --------------------
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


