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
  int iframe = 0;
  const unsigned int nframe_interval = 10;
  while(true){
    if( iframe % nframe_interval == 0 ){
      std::string path_svg;
      if( iframe == nframe_interval*0 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape0.svg"; }
      if( iframe == nframe_interval*1 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape1.svg"; }
      if( iframe == nframe_interval*2 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape2.svg"; }
      if( iframe == nframe_interval*3 ){ path_svg = std::string(PATH_INPUT_DIR)+"/tshirt.svg"; }
      path_svg = std::string(PATH_INPUT_DIR)+"/tshirt.svg";
      std::vector< std::vector<delfem2::CCad2D_EdgeGeo> > aaEdge;
      std::cout << "########################" << std::endl;
      LoopEdgeCCad2D_ReadSVG(aaEdge,
                             path_svg);
      std::cout << "### " << aaEdge.size() << std::endl;
      cad.Clear();
      for(int iae=0;iae<aaEdge.size();++iae){
        std::vector<delfem2::CCad2D_EdgeGeo> aEdge = aaEdge[iae];
        Transform_LoopEdgeCad2D(aEdge,false,true,1.0,1.0);
        std::cout << aEdge.size() << "  " << AreaLoop(aEdge) << std::endl;
        if( AreaLoop(aEdge) < 0 ){ aEdge = InvertLoop(aEdge); }
        aEdge = RemoveEdgeWithZeroLength(aEdge);
        for(auto & ie : aEdge){ ie.GenMesh(-1); }
        std::cout << aEdge.size() << "  " << AreaLoop(aEdge) << std::endl;
        cad.AddFace(aEdge);
      }
      std::cout << Str_SVGPolygon(cad.XY_VtxCtrl_Face(0),1) << std::endl;
      dfm2::CBoundingBox2D bb = cad.BB();
      viewer.nav.camera.trans[0] = -(bb.x_min+bb.x_max)*0.5;
      viewer.nav.camera.trans[1] = -(bb.y_min+bb.y_max)*0.5;
      viewer.nav.camera.trans[2] = 0.0;
      viewer.nav.camera.view_height = 0.5*sqrt( (bb.x_max-bb.x_min)*(bb.x_max-bb.x_min) + (bb.y_max-bb.y_min)*(bb.y_max-bb.y_min) );
      viewer.nav.camera.scale = 1.0;
    }
    iframe = (iframe+1)%(nframe_interval*4);
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
    // --------------------
    viewer.DrawBegin_oldGL();
    cad.iedge_picked = 5;
    cad.ipicked_elem = 1;
    cad.is_draw_face = true;
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


