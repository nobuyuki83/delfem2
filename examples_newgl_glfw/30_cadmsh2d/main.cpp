/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#if defined(_MSC_VER)
#  include <windows.h>
#endif

#ifdef EMSCRIPTEN
#  include <emscripten/emscripten.h>
#  define GLFW_INCLUDE_ES3
#else
#  include <glad/glad.h>
#endif
#include <GLFW/glfw3.h>

#include "delfem2/cagedef.h"
#include "delfem2/glfw/util.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/opengl/new/v23dtricad.h"
#include "delfem2/cad2_dtri2.h"

namespace dfm2 = delfem2;

// -------------------------------------

class CCadDtri_Viewer : public delfem2::glfw::CViewer2 {
 public:
  CCadDtri_Viewer() {
    {
      std::vector<double> aXY = {-1, -1, +1, -1, +1, +1, -1, +1};
      cad.AddPolygon(aXY);
    }
  }
  void InitShader() {
    shdr_cad.Compile();
    shdr_dmsh.Compile();
    {
      shdr_cad.MakeBuffer(cad);
      {
        std::vector<int> aFlgPnt, aFlgTri;
        delfem2::CMesher_Cad2D mesher;
        mesher.edge_length = 0.08;
        mesher.Meshing(dmsh, cad);
        shdr_dmsh.MakeBuffer(dmsh.aVec2, dmsh.aETri);
      }
      {
        std::vector<double> aXYVtx = cad.XY_VtxCtrl_Face(0);
        const unsigned int nxy = dmsh.aVec2.size();
        const auto nv = static_cast<unsigned int>(aXYVtx.size() / 2);
        aW.resize(nxy * nv);
        for (int ixy = 0; ixy < nxy; ++ixy) {
          dfm2::MeanValueCoordinate_Polygon2<dfm2::CVec2d>(
              aW.data() + nv * ixy,
              dmsh.aVec2[ixy].x, dmsh.aVec2[ixy].y,
              aXYVtx.data(), aXYVtx.size() / 2);
          double sum = 0.0;
          for (int iv = 0; iv < nv; ++iv) {
            sum += aW[nv * ixy + iv];
          }
          assert(fabs(sum - 1) < 1.0e-10);
        }
      }
    }

    shdr_cad.is_show_face = false;

    view_height = 1.5;
  }
  void mouse_press(const float src[2]) override {
    cad.Pick(src[0], src[1], view_height);
  }
  void mouse_drag(const float src0[2], const float src1[2]) override {
    if (nav.ibutton == 0) {
      cad.DragPicked(src1[0], src1[1], src0[0], src0[1]);
      shdr_cad.MakeBuffer(cad);
      // --
      std::vector<double> aXYVtx = cad.XY_VtxCtrl_Face(0);
      unsigned int nv = aXYVtx.size() / 2;
      const unsigned int np = dmsh.aVec2.size();
      for (unsigned int ip = 0; ip < np; ++ip) {
        dmsh.aVec2[ip].p[0] = 0.0;
        dmsh.aVec2[ip].p[1] = 0.0;
        for (int iv = 0; iv < nv; ++iv) {
          dmsh.aVec2[ip].p[0] += aW[ip * nv + iv] * aXYVtx[iv * 2 + 0];
          dmsh.aVec2[ip].p[1] += aW[ip * nv + iv] * aXYVtx[iv * 2 + 1];
        }
      }
      shdr_dmsh.MakeBuffer(dmsh.aVec2, dmsh.aETri);
    }
  }
 public:
  delfem2::CCad2D cad;
  delfem2::CMeshDynTri2D dmsh;
  std::vector<double> aW;

  delfem2::opengl::CShader_Cad2D shdr_cad;
  delfem2::opengl::CShader_MeshDTri2D shdr_dmsh;
};

CCadDtri_Viewer viewer;

// -----------------------------------

void draw(GLFWwindow *window) {
  ::glClearColor(0.8, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LESS);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);

  int nw, nh;
  glfwGetFramebufferSize(window, &nw, &nh);
  const float asp = static_cast<float>(nw) / static_cast<float>(nh);
  float mMV[16], mP[16];
  viewer.Mat4_MVP_OpenGL(mMV, mP, asp);
  viewer.shdr_cad.Draw(mP, mMV, viewer.cad);
  viewer.shdr_dmsh.Draw(mP, mMV);

  glfwSwapBuffers(window);
  glfwPollEvents();
}

/*
void callback_key(GLFWwindow *window, int key, int scancode, int action, int mods) {
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, GL_TRUE);
  }
}

void callback_resize(GLFWwindow *window, int width, int height) {
  glViewport(0, 0, width, height);
}
 */

int main() {
  dfm2::glfw::InitGLNew();
  viewer.InitGL();

#ifndef EMSCRIPTEN
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }
#endif

  viewer.InitShader();

#ifdef EMSCRIPTEN
  emscripten_set_main_loop_arg((em_arg_callback_func) draw, viewer.window, 60, 1);
#else
  while (!glfwWindowShouldClose(viewer.window)) { draw(viewer.window); }
#endif

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

