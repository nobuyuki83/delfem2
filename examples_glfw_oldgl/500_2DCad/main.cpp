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
#include "delfem2/cad2d.h"

// --------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_v23.h"
#include "delfem2/opengl/glold_v23dtricad.h"

#ifndef M_PI
  #define M_PI 3.141592653589793
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
  while(!glfwWindowShouldClose(viewer.window)){
    {
      static int iframe = 0;
      const int nframe = 300;
      if( iframe == nframe*0 ){
        cad.Clear();
        const double poly[8] = {-1,-1, +1,-1, +1,+1, -1,+1};
        cad.AddPolygon(std::vector<double>(poly,poly+8));
      }
      else if( iframe == nframe*1 ){
        cad.AddVtxFace(0.0, 0.1, 0);
      }
      else if( iframe == nframe*2 ){
        double param[4] = {0.2, 0.3, -0.2, 0.3};
        std::vector<double> vparam(param,param+4);
        cad.SetEdgeType( 0, 1, vparam );
      }
      else if( iframe == nframe*3 ){
        cad.AddVtxEdge(-0.0, +0.8, 2);
      }
      else if( iframe == nframe*4 ){
        double x0 = 2.1, y0 = 0.0;
        const double poly[8] = {x0-1,y0-1, x0+1,y0-1, x0+1,y0+1, x0-1,y0+1};
        cad.AddPolygon(std::vector<double>(poly,poly+8) );
        cad.AddVtxEdge(x0, -0.2, 5);
      }
      else if( iframe == nframe*5 || iframe == nframe*6 || iframe == nframe*7 )  {
        std::string path_svg;
        if( iframe == nframe*5 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape0.svg"; }
        if( iframe == nframe*6 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape1.svg"; }
        if( iframe == nframe*7 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape2.svg"; }
          //    std::string path_svg = std::string(PATH_INPUT_DIR)+"/shape2.svg";
          //    std::string path_svg = std::string(PATH_INPUT_DIR)+"/shape3.svg";
        std::vector<delfem2::CCad2D_EdgeGeo> aEdge;
        LoopEdgeCCad2D_ReadSVG(aEdge,
                               path_svg);
        Transform_LoopEdgeCad2D(aEdge,false,true,1.0,1.0);
        if( AreaLoop(aEdge) < 0 ){ aEdge = InvertLoop(aEdge); }
        aEdge = RemoveEdgeWithZeroLength(aEdge);
        for(auto & ie : aEdge){ ie.GenMesh(-1); }
        std::cout << aEdge.size() << "  " << AreaLoop(aEdge) << std::endl;
        cad.Clear();
        cad.AddFace(aEdge);
        std::cout << Str_SVGPolygon(cad.XY_VtxCtrl_Face(0),1) << std::endl;
      }
      if( iframe % nframe == 0 ){
        dfm2::CBoundingBox2D bb = cad.BB();
        viewer.nav.camera.trans[0] = -(bb.x_min+bb.x_max)*0.5;
        viewer.nav.camera.trans[1] = -(bb.y_min+bb.y_max)*0.5;
        viewer.nav.camera.trans[2] = 0.0;
        viewer.nav.camera.view_height = 0.5*sqrt( (bb.x_max-bb.x_min)*(bb.x_max-bb.x_min) + (bb.y_max-bb.y_min)*(bb.y_max-bb.y_min) );
        viewer.nav.camera.scale = 1.0;
      }
      iframe = (iframe+1)%(nframe*8);
    }
    // --------------------
    viewer.DrawBegin_oldGL();
    delfem2::opengl::Draw_CCad2D(cad);
    viewer.DrawEnd_oldGL();
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


