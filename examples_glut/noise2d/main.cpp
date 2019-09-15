#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>
#include <assert.h>
#include <vector>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/gl_color.h"
#include "delfem2/noise.h"
#include "delfem2/gl_camera.h"


std::vector<int> aP;
std::vector<double> aGrad;

int nH, nW;
std::vector<double> aV;

////////////////////////////////////////////////////////////////////////////////////

void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
  glutPostRedisplay();
}

void myGlutDisplay(void)
{
  
  //	::glClearColor(0.2, .7, 0.7, 1.0);
  ::glClearColor(0.0, .0, 0.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1, 4.0 );
  
  setGL_Camera2D();
  
  ::glColor3d(1,1,1);
  ::glBegin(GL_LINE_LOOP);
  ::glVertex2d(-0.5, -0.5);
  ::glVertex2d(+0.5, -0.5);
  ::glVertex2d(+0.5, +0.5);
  ::glVertex2d(-0.5, +0.5);
  ::glEnd();
  
  std::vector<std::pair<double, CColor> > colorMap;
  makeHeatMap_BlueCyanGreenYellowRed(colorMap, -0.5, +0.5);
  //  makeHeatMap_BlueGrayRed(colorMap, -0.8, +0.8);
  ::glBegin(GL_QUADS);
  for(int jh=0;jh<nH-1;++jh){
    for(int jw=0;jw<nW-1;++jw){
      int i00 = (jh+0)*nW+(jw+0);
      int i10 = (jh+0)*nW+(jw+1);
      int i11 = (jh+1)*nW+(jw+1);
      int i01 = (jh+1)*nW+(jw+0);
      double v00 = aV[i00];
      double v10 = aV[i10];
      double v11 = aV[i11];
      double v01 = aV[i01];
      double x0 = -0.5+1.0/(nW-1)*(jw+0);
      double x1 = -0.5+1.0/(nW-1)*(jw+1);
      double y0 = -0.5+1.0/(nH-1)*(jh+0);
      double y1 = -0.5+1.0/(nH-1)*(jh+1);
      heatmap(v00, colorMap); ::glVertex2d(x0,y0);
      heatmap(v10, colorMap); ::glVertex2d(x1,y0);
      heatmap(v11, colorMap); ::glVertex2d(x1,y1);
      heatmap(v01, colorMap); ::glVertex2d(x0,y1);
    }
  }
  ::glEnd();
  
  ::glPointSize(5);
  
  // yerrow: input
  ::glLineWidth(1);
  ::glPointSize(5);
  ::glColor3d(1,1,0);
  
  // magenda: last
  ::glLineWidth(1);
  ::glPointSize(5);
  ::glColor3d(1,0,1);
  
  glutSwapBuffers();
}

void myGlutIdle(){
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
}

void ComputePerlin(){
  
  aP.resize(256);
  for(int i=0;i<256;++i){ aP[i]=i; }
  Shuffle(aP);
  aP.resize(512);
  for(int i=0;i<256;++i){ aP[256+i]=i; }
  
  aGrad.clear();
  
  /*
   for(int i=0;i<12;++i){
   double x = (rand()/(RAND_MAX+1.0))-0.5;
   double y = (rand()/(RAND_MAX+1.0))-0.5;
   double len_inv = 1.0/sqrt(x*x+y*y);
   x *= len_inv;
   y *= len_inv;
   aGrad.push_back(x);
   aGrad.push_back(y);
   }
   */
  
  aGrad.push_back(-1); aGrad.push_back(-1);
  aGrad.push_back(-1); aGrad.push_back(+1);
  aGrad.push_back(+1); aGrad.push_back(-1);
  aGrad.push_back(+1); aGrad.push_back(+1);
  
  
  nH = 256;
  nW = 256;
  aV.resize(nH*nW);
  int nrep = 8;
  for(int ih=0;ih<nH;++ih){
    for(int iw=0;iw<nW;++iw){
      double x = (double)iw/(nW)*nrep;
      double y = (double)ih/(nH)*nrep;
      aV[ih*nW+iw] = noise_perlin_2d_oct(x,y,nrep, 4, 0.9, aGrad,aP);
    }
  }
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
    case 'q':
    case 'Q':
    case '\033':  /* '\033' ÇÕ ESC ÇÃ ASCII ÉRÅ[Éh */
      exit(0);
      break;
    case 'a':
      break;
    case ' ':
      ComputePerlin();
      break;
    default:
      break;
  }
}

void myGlutSpecial(int key, int x, int y){
  switch(key){
    case GLUT_KEY_PAGE_UP:
      break;
    case GLUT_KEY_PAGE_DOWN:
      break;
  }
  ::myGlutResize(-1,-1);
  ::glutPostRedisplay();
}

int main(int argc,char* argv[])
{
  // Initailze GLUT
  ::glutInitWindowPosition(200,200);
  ::glutInitWindowSize(400, 300);
  ::glutInit(&argc, argv);
  ::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  ::glutCreateWindow("Cad View");
  
  // Set callback function
  ::glutMotionFunc(myGlutMotion);
  ::glutMouseFunc(myGlutMouse);
  ::glutDisplayFunc(myGlutDisplay);
  ::glutReshapeFunc(myGlutResize);
  ::glutKeyboardFunc(myGlutKeyboard);
  ::glutSpecialFunc(myGlutSpecial);
  ::glutIdleFunc(myGlutIdle);
  
  ComputePerlin();
  
  // Enter main loop
  ::glutMainLoop();
  return 0;
}
