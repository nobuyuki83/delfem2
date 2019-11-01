/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <iostream>
#include <math.h>
#include "delfem2/vec3.h"
#include "delfem2/quat.h"
#include "delfem2/v23m3q.h"

// ---

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl2_color.h"
#include "delfem2/opengl/gl2_v23.h"
#include "../glut_cam.h"

// ----------------------------

class CHandlerRotation{
public:
  CHandlerRotation(){
    size = 1.0;
    quat[0] = 1.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 0.0;
    ielem_picked = -1;
    pos = CVector3(0,0,0);
  }
  void Draw() const{
    DrawHandlerRotation_PosQuat(pos, quat, size, ielem_picked);
  }
  void Pick(bool is_down, double spx, double spy, float mMV[16], float mPj[16], double tol){
    if( !is_down ){
      ielem_picked = -1;
      return;
    }
    CVector3 src = screenUnProjection(CVector3(spx,spy,0), mMV, mPj);
    CVector3 dir = screenUnProjectionDirection(CVector3(0,0,1),mMV, mPj);
    double wh = 1.0/mPj[5];
    ielem_picked = PickHandlerRotation_PosQuat(src,dir, CVector3(0,0,0), quat,size,wh*tol);
  }
  void Drag(double sp0x,double sp0y, double sp1x,double sp1y, float mMV[16], float mPj[16]){
    const CVector2 sp0(sp0x,sp0y);
    const CVector2 sp1(sp1x,sp1y);
    DragHandlerRot_PosQuat(quat,
                   ielem_picked,sp0,sp1,pos,mMV,mPj);
  }
public:
  double size;
  CVector3 pos;
  double quat[4];
  int ielem_picked;
};

// -------------------------------------------

// display data
CHandlerRotation hndlr_rot;
CNav3D_GLUT nav;


void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClearStencil(0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  
  DrawBackground(CColor::Blue());
  
  nav.SetGL_Camera();
  
  ::glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
  {
    ::glPushMatrix();
    double M[16]; Mat4_Quat(M, hndlr_rot.quat);
    glMultMatrixd(M);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,0,0);
    ::glutWireTeapot(1.01);
    ::glEnable(GL_LIGHTING);
    ::glColor3d(1,1,1);
    ::glutSolidTeapot(1.0);
    ::glPopMatrix();
  }
  
  ::glDisable(GL_DEPTH_TEST);
  hndlr_rot.Draw();
  ::glEnable(GL_DEPTH_TEST);
  
  ::glColor3d(0,0,0);
  ShowFPS();
  ::glutSwapBuffers();
}

void myGlutIdle(){
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  ::glViewport(0,0,w,h);
  ::glutPostRedisplay();
}

void myGlutSpecial(int Key, int x, int y)
{
  nav.glutSpecial(Key, x, y);
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y )
{
  nav.glutMotion(x, y);
  double sp1x = nav.mouse_x;
  double sp1y = nav.mouse_y;
  double sp0x = nav.mouse_x-nav.dx;
  double sp0y = nav.mouse_y-nav.dy;
  nav.SetGL_Camera();
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  hndlr_rot.Drag(sp0x,sp0y,sp1x,sp1y,mMV,mPj);
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  nav.glutMouse(button, state, x, y);
  double spx = nav.mouse_x;
  double spy = nav.mouse_y;
  nav.SetGL_Camera();
  float mMV[16]; glGetFloatv(GL_MODELVIEW_MATRIX, mMV);
  float mPj[16]; glGetFloatv(GL_PROJECTION_MATRIX, mPj);
  hndlr_rot.Pick(state==GLUT_DOWN,spx,spy,mMV,mPj,0.05);
  ::glutPostRedisplay();
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
    case 'a':
      break;
    case '1':
      break;
    case '2':
      break;
    case '3':
      break;
    case '4':
      break;
    case 'i': // one iteration
      break;
    case 'd': // change draw mode
      break;
    case 'f': //
      break;
    case 's': //
      break;
    case ' ':
    {
      static int ifile = 0;
      ifile++;
      if( ifile >= 8 ){ ifile=0; }
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
  
  nav.camera.view_height = 2.0;
  nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  setSomeLighting();
  
  glutMainLoop();
  return 0;
}


