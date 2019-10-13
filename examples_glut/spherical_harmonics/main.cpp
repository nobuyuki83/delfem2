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

#include "delfem2/specialfuncs.h"

#include "delfem2/gl2_color.h"
#include "delfem2/gl2_funcs.h"
#include "delfem2/gl2_v23.h"

#include "../glut_cam.h"

static void drawShphere_Heatmap
(double (*value)(double,double,double),
 void (*color)(double) )
{
  const double pi = 3.1415926535;
  int nla = 32; double dl = pi/nla;
  int nlo = 64; double dr = 2.0*pi/nlo;
  ::glBegin(GL_QUADS);
  for(int ila=0;ila<nla-1;ila++){
    int ila0 = ila;
    int ila1 = ila+1;
    double y0 = cos(dl*ila0);
    double y1 = cos(dl*ila1);
    double r0 = sin(dl*ila0);
    double r1 = sin(dl*ila1);
    for(int ilo=0;ilo<nlo;ilo++){
      int ilo0 = ilo;
      int ilo1 = ilo+1;
      double x0 = sin(dr*ilo0);
      double x1 = sin(dr*ilo1);
      double z0 = cos(dr*ilo0);
      double z1 = cos(dr*ilo1);
      CVector3 a(r0*x0,y0,r0*z0);
      CVector3 b(r0*x1,y0,r0*z1);
      CVector3 c(r1*x1,y1,r1*z1);
      CVector3 d(r1*x0,y1,r1*z0);
      {
        double x = (r0+r1)*(x0+x1)*0.25;
        double y = (y0+y1)*0.5;
        double z = (r0+r1)*(z0+z1)*0.25;
        double invl = 1.0/sqrt(x*x+y*y+z*z);
        x *= invl;
        y *= invl;
        z *= invl;
        glNormal3d(x,y,z);
      }
      double va = value(a.x,a.y,a.z);
      double vb = value(b.x,b.y,b.z);
      double vc = value(c.x,c.y,c.z);
      double vd = value(d.x,d.y,d.z);
      color(va+0.5); ::myGlVertex(a);
      color(vb+0.5); ::myGlVertex(b);
      color(vc+0.5); ::myGlVertex(c);
      color(vd+0.5); ::myGlVertex(d);
    }
  }
  ::glEnd();
}

static void drawShphere_Radius
(double (*value)(double,double,double),
 void (*color)(double) )
{
  const double pi = 3.1415926535;
  int nla = 64; double dla = pi/nla;
  int nlo =128; double dlo = 2.0*pi/nlo;
  ::glBegin(GL_TRIANGLES);
  for(int ila=0;ila<nla-1;ila++){
    int ila0 = ila;
    int ila1 = ila+1;
    double y0 = cos(dla*ila0);
    double y1 = cos(dla*ila1);
    double r0 = sin(dla*ila0);
    double r1 = sin(dla*ila1);
    for(int ilo=0;ilo<nlo;ilo++){
      int ilo0 = ilo;
      int ilo1 = ilo+1;
      double x0 = sin(dlo*ilo0);
      double x1 = sin(dlo*ilo1);
      double z0 = cos(dlo*ilo0);
      double z1 = cos(dlo*ilo1);
      CVector3 a(r0*x0,y0,r0*z0);
      CVector3 b(r0*x1,y0,r0*z1);
      CVector3 c(r1*x1,y1,r1*z1);
      CVector3 d(r1*x0,y1,r1*z0);
      double va = value(a.x,a.y,a.z);
      double vb = value(b.x,b.y,b.z);
      double vc = value(c.x,c.y,c.z);
      double vd = value(d.x,d.y,d.z);
      myGlNormal(va*a,vc*c,vb*b);
      color(va+0.5); ::myGlVertex(fabs(va)*a);
      color(vb+0.5); ::myGlVertex(fabs(vb)*b);
      color(vc+0.5); ::myGlVertex(fabs(vc)*c);
      myGlNormal(vd*d,vc*c,va*a);
      color(vd+0.5); ::myGlVertex(fabs(vd)*d);
      color(vc+0.5); ::myGlVertex(fabs(vc)*c);
      color(va+0.5); ::myGlVertex(fabs(va)*a);
    }
  }
  ::glEnd();
}


// -------------------------------

bool is_animation;
CNav3D_GLUT nav;
const int nl = 10;
int l=0;
int m=0;

std::vector<COMPLEX> aDecmpSH;
std::vector<double> aDecmpSHReal;

// -------------------------------

double evaluateSH(double x, double y, double z)
{
  double a[100];
//  makeEvaluationArraySH(l+1,a, x,y,z);
  makeArray_SphericalHarmonics(a, l+1, x,y,z);
  int ish = l*(l+1)+m;
  return a[ish];
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
  
  ::glEnable(GL_LIGHTING);
  drawShphere_Radius(evaluateSH        ,heatmap_glDiffuse);
  
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
    case 'c':
      break;
    case 'd': // change draw mode
      break;
    case 'f': //
      break;
    case 's': //
      break;
    case ' ':
    {
      ++m;
      if( m > l ){ l++; m=-l; }
      if( l >= nl ){ l=0; m=0; }
      std::cout << l << " " << m << " " << (l+1)*l+m << std::endl;
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
  
  // -----------------------------
  
  nav.camera.view_height = 2.0;
  
  setSomeLighting();
  
  glutMainLoop();
	return 0;
}


