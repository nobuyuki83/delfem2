#include <iostream>
#include <sstream>
#include <math.h>

#include <complex>
#include <set>
#include <stack>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/dyntri_v3.h"
#include "delfem2/cad2d.h"

#include "delfem2/funcs_glut.h"
#include "delfem2/funcs_gl.h"
#include "delfem2/v23q_gl.h"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif



///////////////////////////////////////////////////////////////////////////////////////////////


CGlutWindowManager win;
const double view_height = 2.0;
bool is_animation = false;
int imode_draw = 0;

CCad2D cad;

////////////////////////////////////////////////////////////////////////////////////////////////

void myGlutDisplay(void)
{
//  	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
	::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
//	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
  
  win.SetGL_Camera();
  
//  if(      imode_draw == 0 ){ cad.DrawFace_RightSelected(false); }
//  else if( imode_draw == 1 ){ cad.DrawFace_RightSelected(true); }
//  else if( imode_draw == 2 ){ cad.DrawFace_LeftRight(); }
//  cad.DrawVtxEdgeHandler(win.camera.view_height);
  
  cad.Draw();
  
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
  win.glutSpecial(Key,x,y);
}

void myGlutMotion( int x, int y ){
  win.glutMotion(x,y);
  if( win.imodifier != 0){ return; }
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  CVector2 sp0(win.mouse_x-win.dx, win.mouse_y-win.dy);
  CVector2 sp1(win.mouse_x, win.mouse_y);
  const CVector3 src_pick = screenUnProjection(CVector3(sp0.x,sp0.y, 0.0), mMV,mPj);
  const CVector3 dir_pick = screenUnProjectionDirection(CVector3(0.0,  0, -1.0 ), mMV,mPj);
  /////
//  cad.MouseMotion(src_pick,dir_pick, sp0,sp1, mMV, mPj);
}

void myGlutMouse(int button, int state, int x, int y)
{
  win.glutMouse(button,state,x,y);
  if( win.imodifier == GLUT_ACTIVE_SHIFT || win.imodifier == GLUT_ACTIVE_ALT ) return;
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  CVector2 sp0(win.mouse_x, win.mouse_y);
  const CVector3 src_pick = screenUnProjection(CVector3(sp0.x,sp0.y, 0.0), mMV,mPj);
  const CVector3 dir_pick = screenUnProjectionDirection(CVector3(0.0,  0, -1.0 ), mMV,mPj);
  if( state == GLUT_DOWN ){
//    cad.MouseDown(src_pick, dir_pick, sp0, mMV, mPj, view_height);
  }
  if( state == GLUT_UP ){
//    cad.MouseUp(mMV,mPj,view_height);
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
    case 'a':
    {
      is_animation = !is_animation;
      break;
    }
    case 'd':
    {
      imode_draw = (imode_draw+1)%3;
      break;
    }
    case 'b':
    {
      break;
    }
    case 'p':
    {
      break;
    }
    case 'e':
    {
//      MakeItSmooth(cad.aVertex, cad.aEdge, cad.aFace);
      break;
    }
    case 'f':
    {
      break;
    }
    case 's':
    {
//      cad.Initialize_Sphere();
      break;
    }
    case 'c':
    {
//      cad.Initialize_Cube();
      break;
    }
    case 'n':
    {
      break;
    }
    case 't':
    {
      break;
    }
      /*
    case 'w':
    {
      std::ofstream fout;
      fout.open("hoge.txt",std::ios::out);
      cad.WriteFile(fout);
      break;
    }
    case 'r':
    {
      std::ifstream fin;
      fin.open("hoge.txt",std::ios::in);
      cad.ReadFile(fin);
      break;
    }
    case '1':
    {
      cad.imode_edit = CCad3D::EDIT_MOVE;
      break;
    }
    case '2':
    {
      cad.imode_edit = CCad3D::EDIT_ADD_CROSS_SECTION;
      break;
    }
    case '3':
    {
      cad.imode_edit = CCad3D::EDIT_ADD_POINT_EDGE;
      break;
    }
    case '4':
    {
      cad.imode_edit = CCad3D::EDIT_SKETCH;
      break;
    }
    case '5':
    {
      break;
    }
    case '+':
    {
      cad.ChangeEdgeLength(cad.elen*0.9);
      break;
    }
    case '-':
    {
      cad.ChangeEdgeLength(cad.elen/0.9);
      break;
    }
       */
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
  
  ////////////////////////
  win.camera.view_height = view_height;
//  win.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  win.camera.camera_rot_mode = CAMERA_ROT_YTOP;
//    win.camera.camera_rot_mode = CAMERA_ROT_ZTOP;
  
  setSomeLighting();
  
  const double poly[8] = {-1,-1, +1,-1, +1,+1, -1,+1};
  cad.AddPolygon(std::vector<double>(poly,poly+8));
  
  glutMainLoop();
	return 0;
}


