#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>
#include <vector>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/dyntri_v2.h"

double AreaCGCurve(const std::vector<double>& aCV, double cg[2])
{
  const unsigned int nCV = aCV.size()/2;
  double area = 0;
  cg[0] = 0;
  cg[1] = 0;
  for(unsigned int idiv=0;idiv<nCV;idiv++){
    unsigned int ipo0 = idiv;
    unsigned int ipo1 = idiv+1;
    if( idiv == nCV-1 ){ ipo1 = 0; }
    const double p0[2] = { aCV[ipo0*2+0], aCV[ipo0*2+1] };
    const double p1[2] = { aCV[ipo1*2+0], aCV[ipo1*2+1] };
    const double p2[2] = { 0, 0 };
    double a0 = TriArea2D(p0, p1, p2);
    double cg0[2] = { (p0[0]+p1[0]+p2[0])/3.0, (p0[1]+p1[1]+p2[1])/3.0 };
    cg[0] += cg0[0]*a0;
    cg[1] += cg0[1]*a0;
    area += a0;
  }
  cg[0] /= area;
  cg[1] /= area;
  return area;
}

void MakeRandomCV(unsigned int nCV, std::vector<double>& aCV)
{
  aCV.clear();
  for(unsigned int icv=0;icv<nCV;icv++){
    /*
     {
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     aCV.push_back((double)rand()/(RAND_MAX+1.0));
     }
     */
    {
      double tht = icv*3.1415*2.0/nCV;
      double r = (double)rand()/(RAND_MAX+1.0);
      double px = r*sin(tht);
      double py = r*cos(tht);
      aCV.push_back(px);
      aCV.push_back(py);
    }
  }
  
  {
    double cnt[2];
    double area = AreaCGCurve(aCV,cnt);
    if( area < 0 ){
      std::vector<double> aCV0;
      for(unsigned int idiv=0;idiv<nCV;idiv++){
        unsigned int idiv0 = nCV-1-idiv;
        aCV0.push_back( aCV[idiv0*2+0] );
        aCV0.push_back( aCV[idiv0*2+1] );
      }
      aCV = aCV0;
      area = -area;
    }
    if( area < 1.0e-5 ){
      aCV.clear();
      return;
    }
  }
}

void MakeCurveSpline(const std::vector<double>& aCV, std::vector<double>& aVecCurve)
{
  aVecCurve.resize(0);
  const unsigned int nCV = aCV.size()/2;
  unsigned int ndiv = 5;
  for(unsigned int icv=0;icv<nCV;icv++){
    int icv0=icv;   if( icv0 >= nCV ){ icv0-=nCV; }
    int icv1=icv+1; if( icv1 >= nCV ){ icv1-=nCV; }
    int icv2=icv+2; if( icv2 >= nCV ){ icv2-=nCV; }
    const double p0[2] = { aCV[icv0*2+0], aCV[icv0*2+1] };
    const double p1[2] = { aCV[icv1*2+0], aCV[icv1*2+1] };
    const double p2[2] = { aCV[icv2*2+0], aCV[icv2*2+1] };
    for(unsigned int idiv=0;idiv<ndiv;idiv++){
      const double t = 1.0-(double)idiv/ndiv;
      const double w[3] = {0.5*t*t, -t*t + t + 0.5, 0.5*(1-t)*(1-t) };
      const double px = w[0]*p0[0] + w[1]*p1[0] + w[2]*p2[0];
      const double py = w[0]*p0[1] + w[1]*p1[1] + w[2]*p2[1];
      aVecCurve.push_back(px);
      aVecCurve.push_back(py);
    }
  }
}

void myGlVertex2D(const std::vector<double>& vec, unsigned int i)
{
  ::glVertex3d(vec[i*2],vec[i*2+1],+0.5);
}


void drawCurve
(const std::vector<double>& vec,
 const std::vector<double>& aVecCurve0)
{
  if( aVecCurve0.size() < 2 ){ return; }
  ::glBegin(GL_LINES);
  const unsigned int nvec = vec.size()/2;
  for(unsigned int ivec=0;ivec<nvec;ivec++){
    unsigned int jvec = ivec+1; if( jvec >= nvec ){ jvec -= nvec; }
    myGlVertex2D(vec,ivec);
    myGlVertex2D(vec,jvec);
  }
  ::glEnd();
  
  ::glBegin(GL_POINTS);
  for(unsigned int ivec=0;ivec<nvec;ivec++){
    myGlVertex2D(vec,ivec);
  }
  ::glEnd();
}

void drawMesh
(std::vector<int>& aTri,
 std::vector<double>& aXY)
{
  const unsigned int ntri = aTri.size()/3;
  const unsigned int nxys = aXY.size()/2;
  ::glColor3d(1,1,1);
  ::glBegin(GL_TRIANGLES);
  //  double mag = 20;
  for(unsigned int itri=0;itri<ntri;itri++){
    const int ino0 = aTri[itri*3+0];
    const int ino1 = aTri[itri*3+1];
    const int ino2 = aTri[itri*3+2];
    ::glVertex2d( aXY[ino0*2+0], aXY[ino0*2+1] );
    ::glVertex2d( aXY[ino1*2+0], aXY[ino1*2+1] );
    ::glVertex2d( aXY[ino2*2+0], aXY[ino2*2+1] );
  }
  ::glEnd();
  ////////////////
  ::glColor3d(0,0,0);
  ::glBegin(GL_LINES);
  for(unsigned int itri=0;itri<ntri;itri++){
    const int ino0 = aTri[itri*3+0];
    const int ino1 = aTri[itri*3+1];
    const int ino2 = aTri[itri*3+2];
    ::glVertex2d( aXY[ino0*2+0], aXY[ino0*2+1] );
    ::glVertex2d( aXY[ino1*2+0], aXY[ino1*2+1] );
    ::glVertex2d( aXY[ino1*2+0], aXY[ino1*2+1] );
    ::glVertex2d( aXY[ino2*2+0], aXY[ino2*2+1] );
    ::glVertex2d( aXY[ino2*2+0], aXY[ino2*2+1] );
    ::glVertex2d( aXY[ino0*2+0], aXY[ino0*2+1] );
  }
  ::glEnd();
  ////////////////
  ::glPointSize(2);
  ::glColor3d(0,0,0);
  ::glBegin(GL_POINTS);
  for(unsigned int ino=0;ino<nxys;ino++){
    ::glVertex2d( aXY[ino*2+0], aXY[ino*2+1] );
  }
  ::glEnd();
  
}



////////////////////////////////////////////////////////////////////////////////////

std::vector<CEPo2> aPo2D;
std::vector<CVector2> aVec2;
std::vector<ETri> aETri;
std::vector<int> loopIP_ind, loopIP;

int idp_nearest = -1;

int press_button = -1;
double mov_begin_x, mov_begin_y;
bool is_animation = true;
double mag = 1.0;

//////////////////////////////////

void GenMesh(){
  std::vector<double> aCV0; MakeRandomCV(8,aCV0); // current cv
  std::vector<double> aVecCurve0;  MakeCurveSpline(aCV0,aVecCurve0); // current curve
  ////
  std::vector< std::vector<double> > aaXY;
  aaXY.push_back( aVecCurve0 );
  ////
  Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                 aaXY, 0.03);
}

//////////////////////////////////

void myGlutResize(int w, int h)
{
  if( w < 0 ){
    int view[4]; glGetIntegerv(GL_VIEWPORT,view);
    w = view[2];
    h = view[3];
  }
  glViewport(0, 0, w, h);
  ::glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  ::glOrtho(-w/300.0*mag,w/300.0*mag, -h/300.0*mag,h/300.0*mag, -1,1);
  glutPostRedisplay();
}



void myGlutDisplay(void)
{
  ::glClearColor(1.0, 1.0, 1.0, 1.0);
  //  ::glClearColor(0.0, .0, 0.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1, 4.0 );
  
  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
  
  ::glPointSize(5);
  ::glLineWidth(1);
  ::glPointSize(5);
  ::glColor3d(1,1,0);
  
  DrawMeshDynTri_Edge(aETri, aVec2);
  ::glColor3d(0.8, 0.8, 0.8);
  DrawMeshDynTri_FaceNorm(aETri, aVec2);
  
  
  ::glLineWidth(3);
  ::glColor3d(0,0,0);
  for(int iloop=0;iloop<(int)loopIP_ind.size()-1;iloop++){
    ::glBegin(GL_LINE_LOOP);
    for(int iip=loopIP_ind[iloop];iip<loopIP_ind[iloop+1];iip++){
      const int ip = loopIP[iip];
      ::glVertex3d(aVec2[ip].x, aVec2[ip].y, 0.1);
    }
    ::glEnd();
  }
  
  glutSwapBuffers();
}


void myGlutIdle(){
  if( is_animation ){  
    GenMesh();
  }
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
  GLint viewport[4];
  ::glGetIntegerv(GL_VIEWPORT,viewport);
  const int win_w = viewport[2];
  const int win_h = viewport[3];
  const double mov_end_x = (2.0*x-win_w)/win_w;
  const double mov_end_y = (win_h-2.0*y)/win_h;
  mov_begin_x = mov_end_x;
  mov_begin_y = mov_end_y;
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  GLint viewport[4];
  ::glGetIntegerv(GL_VIEWPORT,viewport);
  const int win_w = viewport[2];
  const int win_h = viewport[3];
  mov_begin_x = (2.0*x-win_w)/win_w;
  mov_begin_y = (win_h-2.0*y)/win_h;
  press_button = button;
  if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ){
  }
  else if( button == GLUT_LEFT_BUTTON && state == GLUT_UP){
    press_button = -1;
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
    {
      is_animation = !is_animation;
      break;
    }
    default:
      break;
  }
}

void myGlutSpecial(int key, int x, int y){
  switch(key){
    case GLUT_KEY_PAGE_UP:
      mag *= 1.0/0.9;
      break;
    case GLUT_KEY_PAGE_DOWN:
      mag *= 0.9;
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
  
  // Enter main loop
  ::glutMainLoop();
  return 0;
}
