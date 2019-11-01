/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <sstream>
#include <math.h>
#include <complex>
#include <set>
#include <stack>
#include "delfem2/funcs.h"
#include "delfem2/dtri_v2.h"
#include "delfem2/cad2d.h"

// ------------------------

#if defined(__APPLE__) && defined(__MACH__)
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl2_v23.h"
#include "delfem2/opengl/gl2_v23dtricad.h"
#include "../glut_cam.h"

#ifndef M_PI
  #define M_PI 3.141592653589793
#endif

CNav3D_GLUT nav;
const double view_height = 2.0;
bool is_animation = false;
int imode_draw = 0;

CCad2D cad;
std::vector<CCad2D_EdgeGeo> aEdge;

// --------------------------------------------------------------

void myGlutDisplay(void)
{
//  	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
	::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
//	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
  
  nav.SetGL_Camera();
  
  Draw_CCad2D(cad);
  for(int ie=0;ie<aEdge.size();++ie){
    Draw_CCad2DEdge(aEdge[ie], false, -1);
  }
    
  ::glColor3d(0,0,0);
  ShowFPS();
  
  ::glutSwapBuffers();
}

void myGlutIdle(){
  if( is_animation ){
  }
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
	::glutPostRedisplay();
}

void myGlutSpecial(int Key, int x, int y)
{
  nav.glutSpecial(Key,x,y);
}

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

void myGlutKeyboard(unsigned char Key, int x, int y)
{
	switch(Key)
	{
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
      break;
    case ' ':
    {
      static int istep = 0;
      std::cout << istep << std::endl;
      if( istep == 0 ){
        cad.Clear();
        const double poly[8] = {-1,-1, +1,-1, +1,+1, -1,+1};
        cad.AddPolygon(std::vector<double>(poly,poly+8));
      }
      else if( istep == 1 ){
        cad.AddVtxFace(0.0, 0.1, 0);
      }
      else if( istep == 2 ){
        double param[4] = {0.2, 0.3, -0.2, 0.3};
        std::vector<double> vparam(param,param+4);
        cad.SetEdgeType( 0, 1, vparam );
      }
      else if( istep == 3 ){
        cad.AddVtxEdge(-0.0, +0.8, 2);
      }
      else if( istep == 4 ){
        double x0 = 2.1, y0 = 0.0;
        const double poly[8] = {x0-1,y0-1, x0+1,y0-1, x0+1,y0+1, x0-1,y0+1};
        cad.AddPolygon(std::vector<double>(poly,poly+8) );
        cad.AddVtxEdge(x0, -0.2, 5);
      }
      else if( istep == 5 || istep == 6 || istep == 7 )  {
        std::string path_svg;
        if( istep == 5 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape0.svg"; }
        if( istep == 6 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape1.svg"; }
        if( istep == 7 ){ path_svg = std::string(PATH_INPUT_DIR)+"/shape2.svg"; }
        //    std::string path_svg = std::string(PATH_INPUT_DIR)+"/shape2.svg";
        //    std::string path_svg = std::string(PATH_INPUT_DIR)+"/shape3.svg";
        std::vector<CCad2D_EdgeGeo> aEdge;
        LoopEdgeCCad2D_ReadSVG(aEdge,
                               path_svg);
        Transform_LoopEdgeCad2D(aEdge,false,true,1.0,1.0);
        if( AreaLoop(aEdge) < 0 ){ aEdge = InvertLoop(aEdge); }
        aEdge = RemoveEdgeWithZeroLength(aEdge);
        for(int ie=0;ie<aEdge.size();++ie){ aEdge[ie].GenMesh(-1); }
        std::cout << aEdge.size() << "  " << AreaLoop(aEdge) << std::endl;
        cad.Clear();
        cad.AddFace(aEdge);
        std::cout << Str_SVGPolygon(cad.XY_VtxCtrl_Face(0),1) << std::endl;
      }
      { // set camera
        CBoundingBox2D bb = cad.BB();
        std::cout << "bb: " << bb.x_min << " " << bb.x_max << " " << bb.y_min << " " << bb.y_max << std::endl;
        nav.camera.trans[0] = -(bb.x_min+bb.x_max)*0.5;
        nav.camera.trans[1] = -(bb.y_min+bb.y_max)*0.5;
        nav.camera.trans[2] = 0.0;
        nav.camera.view_height = 0.5*sqrt( (bb.x_max-bb.x_min)*(bb.x_max-bb.x_min) + (bb.y_max-bb.y_min)*(bb.y_max-bb.y_min) );
        nav.camera.scale = 1.0;
      }
      istep = (istep+1)%8;
      break;
    }
    case 'd':
    {
      imode_draw = (imode_draw+1)%3;
      break;
    }
    case 'w':
    {
      WriteCAD_DXF("hoge.dxf",
                   cad,1.0);
    }
  }
  ::glutPostRedisplay();
}

int main(int argc,char* argv[])
{
  glutInit(&argc, argv);
  
	// Initialize GLUT window 3D
  glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
// 	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_STENCIL);
  glutCreateWindow("3D View");
	glutDisplayFunc(myGlutDisplay);
	glutIdleFunc(myGlutIdle);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
  
  // -----------------------------------
  // setting camera
  nav.camera.view_height = view_height;
//  win.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  nav.camera.camera_rot_mode = CAMERA_ROT_YTOP;
//    win.camera.camera_rot_mode = CAMERA_ROT_ZTOP;
  
  setSomeLighting();
  
  const double poly[8] = {-1,-1, +1,-1, +1,+1, -1,+1};
  cad.AddPolygon(std::vector<double>(poly,poly+8));
  
  glutMainLoop();
	return 0;
}


