/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include <iostream>
#include <math.h>

#if defined(__APPLE__) && defined(__MACH__)
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

#include "delfem2/gl2_funcs.h"
#include "delfem2/gl2_color.h"

#include "../glut_cam.h"

// ---------------------------------------------

bool is_animation;
CNav3D_GLUT nav;

// ---------------------------------------------

void drawObject(){ // for shadow
  ::glutSolidTeapot(1.0);
}

void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClearStencil(0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  nav.SetGL_Camera();
  
  DrawBackground();
  
  ::glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  ::glutWireTeapot(1.01);
  ::glEnable(GL_LIGHTING);
  ::glColor3d(1,1,1);
  ::glutSolidTeapot(1.0);
  
  /*
   {
   glDisable(GL_LIGHTING);
   float plane[4] = {0,1,0,5-0.001};
   float lpos[4] = {0,5,0,1};
   float m_shadow[16]; ShadowMatrix(m_shadow, plane, lpos);
   glPushMatrix();
   glMultMatrixf(m_shadow);
   glColor3d(0,0,0);
   glutSolidTeapot(1.0);
   glPopMatrix();
   
   ::glEnable(GL_LIGHTING);
   ::glBegin(GL_QUADS);
   ::glColor3d(1,1,1);
   ::glNormal3d(0,1,0);
   ::glVertex3d(-20,-5,-20);
   ::glVertex3d(+20,-5,-20);
   ::glVertex3d(+20,-5,+20);
   ::glVertex3d(-20,-5,+20);
   ::glEnd();
   }
   */
  
  ::drawFloorShadow(drawObject, -5, 20);
  
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
  ::glViewport(0,0,w,h);
//  std::cout << glutGet(GLUT_WINDOW_BUFFER_SIZE) << std::endl;
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
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  nav.glutMouse(button, state, x, y);
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
      is_animation = !is_animation;
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
  // ----------------------------------
  
  nav.camera.view_height = 2.0;
  
  setSomeLighting();
  
  glutMainLoop();
  return 0;
}


