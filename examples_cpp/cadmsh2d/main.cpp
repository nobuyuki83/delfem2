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

#include "delfem2/dyntri_v2.h"
#include "delfem2/cad2d.h"

#include "delfem2/glut_funcs.h"
#include "delfem2/gl_funcs.h"
#include "delfem2/gl_v23q.h"
#include "delfem2/gl_cad_dyntri_v23.h"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif



///////////////////////////////////////////////////////////////////////////////////////////////


CGlutWindowManager win;
const double view_height = 2.0;
bool is_animation = false;
int imode_draw = 0;
std::vector<double> aXY;
std::vector<unsigned int> aTri;
std::vector<double> aW;
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
  
  Draw_CCad2D(cad);
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0.8, 0.8, 0.8);
  DrawMeshTri2D_Face(aTri,aXY);
  ::glLineWidth(1);
  DrawMeshTri2D_Edge(aTri,aXY);
  
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
  const CVector3 src_pick0 = screenUnProjection(CVector3(sp0.x,sp0.y, 0.0), mMV,mPj);
  const CVector3 src_pick1 = screenUnProjection(CVector3(sp1.x,sp1.y, 0.0), mMV,mPj);
  const CVector3 dir_pick = screenUnProjectionDirection(CVector3(0.0,  0, -1.0 ), mMV,mPj);
  /////
  cad.DragPicked(src_pick1[0],src_pick1[1], src_pick0[0],src_pick0[1]);
  ////
  std::vector<double> aXY_bound = cad.XY_VtxCtrl_Face(0);
  int npb = aXY_bound.size()/2;
  int np = aXY.size()/2;
  for(int ip=0;ip<np;++ip){
    aXY[ip*2+0] = 0.0;
    aXY[ip*2+1] = 0.0;
    for(int ipb=0;ipb<npb;++ipb){
      aXY[ip*+2+0] += aW[ip*npb+ipb]*aXY_bound[ipb*2+0];
      aXY[ip*+2+1] += aW[ip*npb+ipb]*aXY_bound[ipb*2+1];
    }
  }
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
    cad.Pick(src_pick[0],src_pick[1],view_height);
  }
  if( state == GLUT_UP ){
    cad.ivtx_picked = -1;
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
  cad.is_draw_face = false;
  {
    std::vector<int> aFlgPnt, aFlgTri;
    CMeshDynTri2D dmsh;
    CMesher_Cad2D mesher;
    mesher.edge_length = 0.5;
    mesher.Meshing(dmsh, cad);
    dmsh.Export_StlVectors(aXY, aTri);
  }
  
  const int nv = 4;
  std::vector<double> aXY_bound = cad.XY_VtxCtrl_Face(0);
  const int nXY = aXY.size()/2;
  aW.resize(nXY*nv);
  for(int ip=0;ip<nXY;++ip){
    MeanValueCoordinate2D(aW.data()+nv*ip,
                          aXY[ip*2+0], aXY[ip*2+1],
                          aXY_bound.data(), aXY_bound.size()/2);
    double sum = 0.0;
    for(int ipb=0;ipb<aXY_bound.size()/2;++ipb){
      sum += aW[nv*ip+ipb];
    }
    assert( fabs(sum-1)<1.0e-10 );
  }
  glutMainLoop();
	return 0;
}


