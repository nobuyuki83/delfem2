/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/cad2_dtri2.h"
#include <cmath>
#include <stack>
// --------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/v2_glold.h"
#include "delfem2/opengl/cad2dtriv2_glold.h"
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
  while(true){
    {
      const int nframe_interval = 10;
      if( iframe == nframe_interval*0 ){
        cad.Clear();
        const double poly[8] = {-1,-1, +1,-1, +1,+1, -1,+1};
        cad.AddPolygon(std::vector<double>(poly,poly+8));
      }
      else if( iframe == nframe_interval*1 ){
        cad.AddVtxFace(0.0, 0.1, 0);
      }
      else if( iframe == nframe_interval*2 ){
        double param[4] = {0.2, 0.3, 0.8, 0.3};
        std::vector<double> vparam(param,param+4);
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
        dfm2::CBoundingBox2D bb = cad.BB();
        viewer.nav.camera.trans[0] = -(bb.x_min+bb.x_max)*0.5;
        viewer.nav.camera.trans[1] = -(bb.y_min+bb.y_max)*0.5;
        viewer.nav.camera.trans[2] = 0.0;
        viewer.nav.camera.view_height = 0.5*bb.LengthDiagonal();
        viewer.nav.camera.scale = 1.0;
      }
      iframe = (iframe+1)%(nframe_interval*5);
    }
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


