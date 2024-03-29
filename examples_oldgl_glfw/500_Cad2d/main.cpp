/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/cad2.h"
#include "delfem2/cad2_mesher.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/cad2dtriv2.h"

#ifndef M_PI
#  define M_PI 3.141592653589793
#endif

namespace dfm2 = delfem2;

// -------------------------------------

int main()
{
  delfem2::CCad2D cad;
  // --------------------
  delfem2::glfw::CViewer3 viewer(2);
  delfem2::glfw::InitGLOld();
  viewer.OpenWindow();
  delfem2::opengl::setSomeLighting();
  for(int iproblem=0;;iproblem=(iproblem+1)%5){
    if( iproblem == 0 ){
      cad.Clear();
      cad.AddPolygon(std::vector<double>{
        -1,-1,
        +1,-1,
        +1,+1,
        -1,+1});
    }
    else if( iproblem == 1 ){
      cad.AddVtxFace(0.0, 0.1, 0);
    }
    else if( iproblem == 2 ){
      cad.aEdge[0].SetCubicBezierCurve({-0.5,-1.5}, {+0.5, -1.5});
    }
    else if( iproblem == 3 ){
      cad.AddVtxEdge(-0.0, +0.8, 2);
    }
    else if( iproblem == 4 ){
      double x0 = 2.1, y0 = 0.0;
      cad.AddPolygon(std::vector<double>{
        x0-1, y0-1,
        x0+1, y0-1,
        x0+1, y0+1,
        x0-1, y0+1} );
      cad.AddVtxEdge(x0, -0.2, 5);
    }
    {
      dfm2::CBoundingBox2<double> bb = cad.BB();
      viewer.trans[0] = -(bb.x_min+bb.x_max)*0.5;
      viewer.trans[1] = -(bb.y_min+bb.y_max)*0.5;
      viewer.trans[2] = 0.0;
      viewer.scale = 1.0;
    }
    // --------------------
    delfem2::CMeshDynTri2D dmesh;
    {
      delfem2::CMesher_Cad2D mesher;
      mesher.edge_length = -1;
      mesher.Meshing(dmesh, cad);
    }
    for(unsigned int iframe=0;iframe<100;++iframe) {
      if( glfwWindowShouldClose(viewer.window) ){ break; }
      viewer.DrawBegin_oldGL();
      delfem2::opengl::DrawCad2Vtxs(cad, UINT_MAX);
      delfem2::opengl::DrawCad2Edges(cad, UINT_MAX);
      delfem2::opengl::DrawMeshDynTri_Edge(dmesh.aETri,dmesh.aVec2);
      ::glColor3d(0.8, 0.8, 0.8);
      delfem2::opengl::DrawMeshDynTri_FaceNorm(dmesh.aETri,dmesh.aVec2);
      viewer.SwapBuffers();
      glfwPollEvents();
    }
    if( glfwWindowShouldClose(viewer.window) ){ break; }
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


