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

#include "delfem2/cad2_mesh_deformation.h"
#include "delfem2/cagedef.h"
#include "delfem2/glfw/util.h"
#include "delfem2/glfw/viewer2.h"
#include "delfem2/opengl/new/drawer_cad2.h"
#include "delfem2/opengl/new/drawer_dynmesh.h"
#include "delfem2/cad2_dtri2.h"

namespace dfm2 = delfem2;

// -------------------------------------

class CCadDtri_Viewer : public delfem2::glfw::CViewer2 {
 public:
  CCadDtri_Viewer() {
    // std::vector<double> aXY = {-1, -1, +1, -1, +1, +1, -1, +1};
    // cad.AddPolygon(aXY);
    cad.AddPolygon({
      -0.84, +1.20,
      -1.00, +0.25,
      -1.00, -1.40,
      +1.00, -1.40,
      +1.00, +0.25,
      +0.84, +1.20,
      +0.30, +1.40,
      -0.30, +1.40});
    cad.aEdge[0].SetCubicBezierCurve({-0.74,+1.0}, {-0.7, +0.45});
    cad.aEdge[4].SetCubicBezierCurve({+0.7,+0.45},{0.74,+1.0});
    cad.aEdge[6].SetCubicBezierCurve({+0.1,+1.2},{-0.1,+1.2});
    cad.Tessellation();    
  }
  
  void InitShader() {
    shdr_cad.Compile();
    shdr_dmsh.Compile();
    shdr_cad.MakeBuffer(cad);
    mesher.edge_length = 0.08;
    mesher.Meshing(dmsh, cad);
    shdr_dmsh.MakeBuffer(dmsh.aVec2, dmsh.aETri);
    shdr_cad.is_show_face = false;
    view_height = 1.5;
  }
  void mouse_press(const float src[2]) override {
    cad.Pick(src[0], src[1], view_height);
    picked_pos = {src[0], src[1]};
    {
      std::vector<unsigned int> tri_vtx;
      dmsh.Export_StlVectors(vtx_xy_when_picked, tri_vtx);
    }
    if (cad.ivtx_picked != UINT_MAX) { // vertex is picked
      picked_pos = {
        static_cast<float>(cad.aVtx[cad.ivtx_picked].pos.x),
        static_cast<float>(cad.aVtx[cad.ivtx_picked].pos.y)};
    }
    SetCadMeshDeformationWeight(vtx_w,
                                cad, mesher, vtx_xy_when_picked);
  }
  void mouse_drag(const float src0[2], const float src1[2]) override {
    if (nav.ibutton == 0) {
      cad.DragPicked(src1[0], src1[1], src0[0], src0[1]);
      shdr_cad.MakeBuffer(cad);
      // --
      for (unsigned int ip = 0; ip < dmsh.aVec2.size(); ++ip) {
        dmsh.aVec2[ip].p[0] = vtx_xy_when_picked[ip * 2 + 0] + vtx_w[ip] * (src1[0] - picked_pos[0]);
        dmsh.aVec2[ip].p[1] = vtx_xy_when_picked[ip * 2 + 1] + vtx_w[ip] * (src1[1] - picked_pos[1]);
      }
      shdr_dmsh.MakeBuffer(dmsh.aVec2, dmsh.aETri);
    }
  }
 public:
  delfem2::CCad2D cad;
  delfem2::CMeshDynTri2D dmsh;
  delfem2::CMesher_Cad2D mesher;
  std::vector<double> vtx_w;
  std::vector<double> vtx_xy_when_picked;
  std::array<float, 2> picked_pos{0.f, 0.f};
  //
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

  dfm2::CMat4f mMV = viewer.GetModelViewMatrix();
  dfm2::CMat4f mP = viewer.GetProjectionMatrix();
  const dfm2::CMat4f mZ = dfm2::CMat4f::ScaleXYZ(1,1,-1);
  viewer.shdr_cad.Draw((mZ*mP).transpose().data(), mMV.transpose().data(), viewer.cad);
  viewer.shdr_dmsh.Draw((mZ*mP).transpose().data(), mMV.transpose().data());
  glfwSwapBuffers(window);
  glfwPollEvents();
}

int main() {
  dfm2::glfw::InitGLNew();
  
  viewer.OpenWindow();
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

