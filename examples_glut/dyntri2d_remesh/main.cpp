/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>
#include <vector>
#include <set>

#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

#include "delfem2/dtri_v2.h"
#include "delfem2/gl_cad_dyntri_v23.h"

// --------------------------------------------

std::vector<CEPo2> aPo2D;
std::vector<CVector2> aVec2;
std::vector<ETri> aETri;
std::vector<int> loopIP_ind, loopIP;

int idp_nearest = -1;

int press_button = -1;
double mov_begin_x, mov_begin_y;
bool is_animation = true;
double mag = 1.0;

// --------------------------------------------

void Refine(double px, double py)
{
  CCmdRefineMesh aCmd;
  RefinementPlan_EdgeLongerThan_InsideCircle(aCmd,
                                             0.05, px, py, 0.1,
                                             aPo2D, aVec2, aETri);
  RefineMesh(aPo2D, aETri, aVec2, aCmd);
}

void Coarse(double px, double py)
{
  for(int ip=aPo2D.size()-1;ip>=0;--ip){
    if( aPo2D[ip].e == -1 ){ continue; }
    if( Distance(aVec2[ip],CVector2(px,py)) > 0.1 ){ continue; }
    std::vector< std::pair<int,int> > aTriSuP;
    GetTriArrayAroundPoint(aTriSuP,
                           ip,aPo2D,aETri);
    std::vector<int> aPSuP(aTriSuP.size());
    const int npsup = aPSuP.size();
    for(int iit=0;iit<npsup;++iit){
      const int itri0 = aTriSuP[iit].first;
      const int inotri0 = aTriSuP[iit].second;
      assert( ip == aETri[itri0].v[inotri0] );
      aPSuP[iit] = aETri[itri0].v[ (inotri0+2)%3 ];
    }
    std::map<double,int> mapDistTri;
    for(int iit=0;iit<npsup;++iit){
      int ip1 = aPSuP[iit];
      double d01 = Distance(aVec2[ip],aVec2[ip1]);
      double min_area = 0.0;
      for(int jjt=0;jjt<npsup-2;++jjt){
        const int ip2 = aPSuP[(iit+jjt+1)%npsup];
        const int ip3 = aPSuP[(iit+jjt+2)%npsup];
        double area = TriArea(aVec2[ip1],aVec2[ip2],aVec2[ip3]);
        if( jjt == 0 || area < min_area ){ min_area = area; }
      }
      if( min_area > 1.0e-10 ){
        mapDistTri.insert( std::make_pair(d01,iit) );
      }
    }
    if( mapDistTri.size() > 0 ){
      const int iit0 = mapDistTri.begin()->second;
      assert( iit0>=0 && iit0 < aPSuP.size() );
      double dist0 = mapDistTri.begin()->first;
      if( dist0 < 0.05 ){
        const int itri0 = aTriSuP[iit0].first;
        const int inotri0 = aTriSuP[iit0].second;
        Collapse_ElemEdge(itri0, (inotri0+1)%3, aPo2D, aETri);
        const int ip1 = aETri[itri0].v[ (inotri0+2)%3 ];
        DelaunayAroundPoint(ip1, aPo2D, aETri, aVec2);
      }
    }
  }
}


void GenMesh()
{
  std::vector< std::vector<double> > aaXY;
  aaXY.resize(1);
  {
    double xys[8] = {-0.5,-0.5, +0.5,-0.5, +0.5,+0.5, -0.5,+0.5};
    aaXY[0].assign(xys,xys+8);
  }
  /////
  const double elen = 0.11;
  {
    JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                         aaXY);
    if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
      return;
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind,aVec2);
    if( elen > 10e-10 ){
      ResamplingLoop(loopIP_ind,loopIP,aVec2,
                     elen );
    }
  }
  ////
  Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                 loopIP_ind,loopIP);
  if( elen > 1.0e-10 ){
    CInputTriangulation_Uniform param(1.0);
    std::vector<int> aFlgPnt(aPo2D.size()), aFlgTri(aETri.size());
    MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
                  aVec2.size(), 0, elen, param);
  }
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
  if( is_animation ){}
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
    Refine(mov_begin_x, mov_begin_y);
  }
  if( button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN ){
    Coarse(mov_begin_x, mov_begin_y);
  }
  if( state == GLUT_UP){
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
    case 'b':
    {

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
  
  GenMesh();
  
  // Enter main loop
  ::glutMainLoop();
  return 0;
}
