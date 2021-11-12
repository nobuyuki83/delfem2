/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <stack>
#include <filesystem>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/glfw/viewer2.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/cad2dtriv2.h"
#include "delfem2/openglstb/glyph.h"
#include "delfem2/cad2.h"
#include "delfem2/cad2_mesher.h"
#include "delfem2/cad2_io_svg.h"

#ifndef M_PI
#  define M_PI 3.141592653589793
#endif

namespace dfm2 = delfem2;

// -------------------------------------

int main() {
  delfem2::openglstb::CGlyph glyph(
      std::filesystem::path(PATH_INPUT_DIR) / "myFont.png");
  glyph.ParseGlyphInfo(
      std::filesystem::path(PATH_INPUT_DIR) / "myFont.fnt");
  delfem2::CCad2D cad;
  // --------------------
  delfem2::glfw::CViewer2 viewer;
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  glyph.InitGL();
  delfem2::opengl::setSomeLighting();
  for (int iproblem = 0;;iproblem=(iproblem+1)%6) {
    {
      std::filesystem::path path_svg;
      if (iproblem == 0) { path_svg = std::filesystem::path(PATH_INPUT_DIR) / "shape0.svg"; }
      if (iproblem == 1) { path_svg = std::filesystem::path(PATH_INPUT_DIR) / "shape1.svg"; }
      if (iproblem == 2) { path_svg = std::filesystem::path(PATH_INPUT_DIR) / "shape2.svg"; }
      if (iproblem == 3) { path_svg = std::filesystem::path(PATH_INPUT_DIR) / "tshirt.svg"; }
      if (iproblem == 4) { path_svg = std::filesystem::path(PATH_INPUT_DIR) / "ltshirt.svg"; }
      if (iproblem == 5) { path_svg = std::filesystem::path(PATH_INPUT_DIR) / "lraglan.svg"; }
      dfm2::ReadSVG_Cad2D(
          cad,
          path_svg.string(), 1.0);
      dfm2::CBoundingBox2<double> bb = cad.BB();
      viewer.trans[0] = static_cast<float>(-(bb.x_min + bb.x_max) * 0.5);
      viewer.trans[1] = static_cast<float>(-(bb.y_min + bb.y_max) * 0.5);
      viewer.view_height = static_cast<float>(0.5 * bb.LengthDiagonal());
      viewer.scale = 1.0;
      cad.iedge_picked = -1;
    }
    if (glfwWindowShouldClose(viewer.window)) { break; }
    delfem2::CMeshDynTri2D dmesh;
    {
      delfem2::CMesher_Cad2D mesher;
      mesher.edge_length = -1;
      mesher.Meshing(dmesh, cad);
    }
    for(unsigned int iframe=0;iframe<30;++iframe) {
      viewer.DrawBegin_oldGL();
      delfem2::opengl::Draw_CCad2D(cad);
      delfem2::opengl::DrawMeshDynTri_Edge(dmesh.aETri,dmesh.aVec2);
      ::glColor3d(0.8, 0.8, 0.8);
      delfem2::opengl::DrawMeshDynTri_FaceNorm(dmesh.aETri,dmesh.aVec2);
      {
        ::glTranslated(0, 0, +0.9);
        for (unsigned int ie = 0; ie < cad.topo.edges.size(); ++ie) {
          unsigned int iv0 = cad.topo.edges[ie].iv0;
          unsigned int iv1 = cad.topo.edges[ie].iv1;
          dfm2::CVec2d p = (cad.aVtx[iv0].pos + cad.aVtx[iv1].pos) * 0.5;
          glyph.DrawStringAt(std::to_string(ie), 0.8, p.x, p.y);
        }
        ::glTranslated(0, 0, -0.9);
      }
      viewer.SwapBuffers();
      glfwPollEvents();
    }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}







/*
 void myGlutMotion( int x, int y ){
 nav.glutMotion(x,y);
 if( nav.imodifier != 0){ return; }
 float px0,py0, px1,py1; nav.PosMove2D(px0,py0, px1,py1);
 cad.DragPicked(px1,py1, px0,py0);
 }
 
 void myGlutMouse(int button, int state, int x, int y)
 {
 nav.glutMouse(button,state,x,y);
 if( nav.imodifier == GLUT_ACTIVE_SHIFT || nav.imodifier == GLUT_ACTIVE_ALT ) return;
 if( state == GLUT_DOWN ){
 float px, py; nav.PosMouse2D(px, py);
 cad.Pick(px, py, nav.camera.view_height);
 }
 }
 */


