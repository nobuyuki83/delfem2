/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#include <stack>
#include "delfem2/funcs.h"
#include "delfem2/cad2d_v2dtri.h"

// --------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_v23.h"
#include "delfem2/opengl/glold_v23dtricad.h"
//
#include "delfem2/opengl/glfw/viewer_glfw.h"

#ifndef M_PI
#  define M_PI 3.141592653589793
#endif

namespace dfm2 = delfem2;

// -------------------------------------

int main(int argc,char* argv[])
{
  delfem2::CCad2D cad;
  // --------------------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  delfem2::opengl::setSomeLighting();
  int iframe = 0;
  const unsigned int nframe_interval = 10;
  while(true){
    if( iframe % nframe_interval == 0 ){
      std::string path_svg;
      if( iframe == nframe_interval*0 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape0.svg"; }
      if( iframe == nframe_interval*1 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape1.svg"; }
      if( iframe == nframe_interval*2 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape2.svg"; }
      if( iframe == nframe_interval*3 ){ path_svg = std::string(PATH_INPUT_DIR)+"/tshirt.svg"; }
      dfm2::ReadSVG_Cad2D(cad,
                          path_svg, 1.0);
//      std::cout << Str_SVGPolygon(cad.XY_VtxCtrl_Face(0),1) << std::endl;
      dfm2::CBoundingBox2D bb = cad.BB();
      viewer.nav.camera.trans[0] = -(bb.x_min+bb.x_max)*0.5;
      viewer.nav.camera.trans[1] = -(bb.y_min+bb.y_max)*0.5;
      viewer.nav.camera.trans[2] = 0.0;
      viewer.nav.camera.view_height = 0.5*bb.LengthDiagonal();
      viewer.nav.camera.scale = 1.0;
      cad.iedge_picked = 22;
    }
    iframe = (iframe+1)%(nframe_interval*4);
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
    // --------------------
    viewer.DrawBegin_oldGL();
    delfem2::opengl::Draw_CCad2D(cad);
    viewer.DrawEnd_oldGL();
  }
EXIT:
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


