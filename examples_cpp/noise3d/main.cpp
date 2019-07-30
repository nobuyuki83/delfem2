#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>
#include <assert.h>
#include <vector>

#ifdef __APPLE__
#include <GL/glew.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include "delfem2/noise.h"

#include "delfem2/gl_color.h"
#include "delfem2/gl_funcs.h"
#include "delfem2/glut_funcs.h"

////////////////////////////////////////////////////////////////////////////////////


CGlutWindowManager win;

std::vector<int> aP;
std::vector<double> aGrad;

int nH, nW, nD;
std::vector<unsigned char> aV;

////////////////////////////////////////////////////////////////////////////////////

void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
  glutPostRedisplay();
}

void myGlutDisplay(void)
{
  int viewport[4]; glGetIntegerv(GL_VIEWPORT,viewport);
  
  ::glClearColor(0.2, .7, 0.7, 1.0);
  //	::glClearColor(0.0, .0, 0.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1, 4.0 );
  
  win.SetGL_Camera();
  
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_TEXTURE_3D);
  glEnable(GL_TEXTURE_GEN_S);
  glEnable(GL_TEXTURE_GEN_T);
  glEnable(GL_TEXTURE_GEN_R);
  ::glutSolidTeapot(1.0);
  //    glutSolidSphere(1.0, 32, 16);
  //  glutSolidDodecahedron();
  ::glDisable(GL_TEXTURE_3D);
  glDisable(GL_TEXTURE_GEN_S);
  glDisable(GL_TEXTURE_GEN_T);
  glDisable(GL_TEXTURE_GEN_R);
  
  ShowFPS();
  
  glutSwapBuffers();
}

void myGlutIdle(){
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
  win.glutMotion(x,y);
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  win.glutMouse(button, state, x, y);
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
  
  aGrad.push_back(-1); aGrad.push_back(-1); aGrad.push_back(+0);
  aGrad.push_back(-1); aGrad.push_back(+1); aGrad.push_back(+0);
  aGrad.push_back(+1); aGrad.push_back(-1); aGrad.push_back(+0);
  aGrad.push_back(+1); aGrad.push_back(+1); aGrad.push_back(+0);
  aGrad.push_back(+0); aGrad.push_back(-1); aGrad.push_back(-1);
  aGrad.push_back(+0); aGrad.push_back(-1); aGrad.push_back(+1);
  aGrad.push_back(+0); aGrad.push_back(+1); aGrad.push_back(-1);
  aGrad.push_back(+0); aGrad.push_back(+1); aGrad.push_back(+1);
  aGrad.push_back(-1); aGrad.push_back(+0); aGrad.push_back(-1);
  aGrad.push_back(-1); aGrad.push_back(+0); aGrad.push_back(+1);
  aGrad.push_back(+1); aGrad.push_back(+0); aGrad.push_back(-1);
  aGrad.push_back(+1); aGrad.push_back(+0); aGrad.push_back(+1);
  
  
  nH = 128;
  nW = 128;
  nD = 128;
  aV.resize(nH*nW*nD*4);
  int nrep = 4;
  for(int id=0;id<nD;++id){
    for(int ih=0;ih<nH;++ih){
      for(int iw=0;iw<nW;++iw){
        double x = (double)iw/nH*nrep;
        double y = (double)ih/nW*nrep;
        double z = (double)id/nD*nrep;
        double v = noise_perlin_3d_oct(x,y,z,nrep, 4,0.5, aGrad,aP);
        //        double v = noise_perlin_3d(x,y,z, aGrad,aP);
        double v0 = v*128+128;
        if( v0 < 0   ){ v0 =   0; }
        if( v0 > 255 ){ v0 = 255; }
        unsigned char ucv = (unsigned char)v0;
        aV[(id*nW*nH+ih*nW+iw)*4+0] = ucv;
        aV[(id*nW*nH+ih*nW+iw)*4+1] = ucv;
        aV[(id*nW*nH+ih*nW+iw)*4+2] = ucv;
        aV[(id*nW*nH+ih*nW+iw)*4+3] = 255;
      }
    }
  }
  
  /*
   Noise3 noise(5, 5, 5);
   for (int id = 0; id < nD; ++id) {
   double p = (double)id / (double)nD;
   for (int ih = 0; ih < nH; ++ih) {
   double q = (double)ih / (double)nH;
   for (int iw = 0; iw < nW; ++iw) {
   double r = (double)iw / (double)nW;
   GLubyte ub = (GLubyte)(noise.noise(p, q, r) * 127.5 + 127.5);
   //        texture[k][j][i][0] = texture[k][j][i][1] = texture[k][j][i][2] = ub;
   //        texture[k][j][i][3] = 255;
   aV[(id*nW*nH+ih*nW+iw)*4+0] = ub;
   aV[(id*nW*nH+ih*nW+iw)*4+1] = ub;
   aV[(id*nW*nH+ih*nW+iw)*4+2] = ub;
   aV[(id*nW*nH+ih*nW+iw)*4+3] = 255;
   }
   }
   }
   */
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
      glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, nW, nH, nD, 0,
                   GL_RGBA, GL_UNSIGNED_BYTE, aV.data());
      break;
    default:
      break;
  }
}

void myGlutSpecial(int key, int x, int y){
  win.glutSpecial(key, x, y);
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
  
  glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
  glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
  glTexGeni(GL_R, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
  static double genfunc[][4] = {
    { 1.0, 0.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0, 0.0 },
    { 0.0, 0.0, 1.0, 0.0 },
  };
  glTexGendv(GL_S, GL_OBJECT_PLANE, genfunc[0]);
  glTexGendv(GL_T, GL_OBJECT_PLANE, genfunc[1]);
  glTexGendv(GL_R, GL_OBJECT_PLANE, genfunc[2]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, nW, nH, nD, 0,
               GL_RGBA, GL_UNSIGNED_BYTE, aV.data());
  
  win.camera.view_height = 2.0;
  
  setSomeLighting();
  
  // Enter main loop
  ::glutMainLoop();
  return 0;
}
