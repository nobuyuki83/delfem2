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

#include "delfem2/cad2_dtri2.h"
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
  int iframe = 0;
  while(true){
    {
      const int nframe_interval = 10;
      if( iframe == nframe_interval*0 ){
        cad.Clear();
        cad.AddPolygon(std::vector<double>{-1,-1, +1,-1, +1,+1, -1,+1});
      }
      else if( iframe == nframe_interval*1 ){
        cad.AddVtxFace(0.0, 0.1, 0);
      }
      else if( iframe == nframe_interval*2 ){
        std::vector<double> vparam{0.2, 0.3, 0.8, 0.3};
        cad.SetEdgeType( 0, dfm2::CCad2D_EdgeGeo::BEZIER_CUBIC, vparam );
      }
      else if( iframe == nframe_interval*3 ){
        cad.AddVtxEdge(-0.0, +0.8, 2);
      }
      else if( iframe == nframe_interval*4 ){
        double x0 = 2.1, y0 = 0.0;
        const double poly[8] = {x0-1,y0-1, x0+1,y0-1, x0+1,y0+1, x0-1,y0+1};
        cad.AddPolygon(std::vector<double>(poly,poly+8) );
        cad.AddVtxEdge(x0, -0.2, 5);
      }
      if( iframe % nframe_interval == 0 ){
        dfm2::CBoundingBox2<double> bb = cad.BB();
        viewer.trans[0] = -(bb.x_min+bb.x_max)*0.5;
        viewer.trans[1] = -(bb.y_min+bb.y_max)*0.5;
        viewer.trans[2] = 0.0;
        //viewer.projection.view_height = 0.5*bb.LengthDiagonal();
        viewer.scale = 1.0;
      }
      iframe = (iframe+1)%(nframe_interval*5);
    }
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
    // --------------------
    viewer.DrawBegin_oldGL();
    delfem2::opengl::Draw_CCad2D(cad);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


