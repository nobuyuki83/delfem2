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

#include "delfem2/dtri_v2.h"
#include "delfem2/cad2d.h"

#include "delfem2/gl2_funcs.h"
#include "delfem2/gl2_v23.h"
#include "delfem2/gl2_v23dtricad.h"

#include "../glut_cam.h"

#ifndef M_PI
  #define M_PI 3.141592653589793
#endif

// ----------------------------------------------

CNav3D_GLUT nav;
const double view_height = 2.0;
bool is_animation = false;
int imode_draw = 0;
std::vector<double> aXY;
std::vector<unsigned int> aTri;
std::vector<double> aW;
CCad2D cad;

// ----------------------------------------------

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
  nav.glutSpecial(Key,x,y);
}

void myGlutMotion( int x, int y ){
  nav.glutMotion(x,y);
  if( nav.imodifier != 0){ return; }
  float x0,y0, x1,y1; nav.PosMove2D(x0,y0, x1,y1);
  cad.DragPicked(x1,y1, x0,y0);
  // ----
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
  nav.glutMouse(button,state,x,y);
  if( nav.imodifier == GLUT_ACTIVE_SHIFT || nav.imodifier == GLUT_ACTIVE_ALT ) return;
  if( state == GLUT_DOWN ){
    float x0,y0; nav.PosMouse2D(x0, y0);
    cad.Pick(x0,y0,view_height);
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
  
  // -----------------
  nav.camera.view_height = view_height;
//  win.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  nav.camera.camera_rot_mode = CAMERA_ROT_YTOP;
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


