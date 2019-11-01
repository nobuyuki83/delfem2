#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>
#include <vector>
#include <set>
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mshtopo.h"
#include "delfem2/dtri.h"

#include "delfem2/objfunc_v23.h"
#include "delfem2/objfunc_v23dtri.h"
#include "delfem2/dtri_v2.h"

// --------------

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/opengl/gl2_v23dtricad.h"
#include "delfem2/opengl/gl2_funcs.h"
#include "../glut_cam.h"


// ------------------------------------

std::vector<CEPo2> aPo2D;
std::vector<CVector2> aVec2;
std::vector<ETri> aETri;
std::vector<unsigned int> aLine;
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aXYZt;
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
const double mass_point = 0.01;
const double dt = 0.01;
const double gravity[3] = {0.0, 0.0, -10.0};

CNav3D_GLUT nav;
bool is_animation = false;

// -------------------------------------

void StepTime()
{
  PBD_Pre3D(aXYZt,
            dt, gravity, aXYZ, aUVW, aBCFlag);
  PBD_TriStrain(aXYZt.data(),
                aXYZt.size()/3, aETri, aVec2);
  PBD_Bend(aXYZt.data(),
           aXYZt.size()/3, aETri, aVec2);
  PBD_Seam(aXYZt.data(),
           aXYZt.size()/3, aLine.data(), aLine.size()/2);
  PBD_Post(aXYZ, aUVW,
           dt, aXYZt, aBCFlag);

}

//////////////////////////////////

void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
  ::glutPostRedisplay();
}



void myGlutDisplay(void)
{
  ::glClearColor(1.0, 1.0, 1.0, 1.0);
  //  ::glClearColor(0.0, .0, 0.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1, 4.0 );
  
  nav.SetGL_Camera();
  
  ::glPointSize(5);
  ::glLineWidth(1);
  {
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0.8, 0.8, 0.8);
    /*
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
    ::glEnable(GL_DEPTH_TEST);
     */
//    DrawMeshTri3D_FaceNorm(aXYZ, aTri);
  }
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  DrawMeshDynTri3D_Edge(aXYZ, aETri);

  
  ShowFPS();
  
  glutSwapBuffers();
}


void myGlutIdle(){
  if( is_animation ){
    StepTime();
  }
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y )
{
  nav.glutMotion(x,y);
}

void myGlutMouse(int button, int state, int x, int y)
{
  nav.glutMouse(button, state, x,y);
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
  }
}

void myGlutSpecial(int key, int x, int y){
  nav.glutSpecial(key,x,y);
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
  
  {
    std::vector< std::vector<double> > aaXY;
    aaXY.resize(1);
    double xys[12] = {
      -0.5,-0.5,
      +0.5,-0.5,
      +0.5,+0.5,
      +0.1,+0.6,
      -0.1,+0.6,
      -0.5,+0.5,
    };
    aaXY[0].assign(xys,xys+12);
    GenMesh(aPo2D,aETri,aVec2,
            aaXY, 0.05, 0.05);
  }
  /////////
  const int np = aPo2D.size();
  aXYZ.resize(np*3);
  for(int ip=0;ip<np;++ip){
    aXYZ[ip*3+0] = aVec2[ip].x;
    aXYZ[ip*3+1] = aVec2[ip].y;
    aXYZ[ip*3+2] = 0.0;
  }
  aXYZt = aXYZ;
  aUVW.resize(np*3,0.0);
  aBCFlag.resize(np,0);
  for(int ip=0;ip<np;++ip){
    if( aXYZ[ip*3+1]  > +0.59 ){
      aBCFlag[ip] = 1;
    }
  }
  aLine.clear();
  { // make aLine
    std::map<int,int> mapY2Ip;
    for(int ip=0;ip<np;++ip){
      if( aXYZ[ip*3+0]  > +0.49 ){
        double y0 = aXYZ[ip*3+1];
        int iy = (int)(y0/0.0132);
        mapY2Ip[iy] = ip;
      }
    }
    for(int ip=0;ip<np;++ip){
      if( aXYZ[ip*3+0]  < -0.49 ){
        double y1 = aXYZ[ip*3+1];
        int iy = (int)(y1/0.0132);
        assert( mapY2Ip.find(iy) != mapY2Ip.end() );
        int ip0 = mapY2Ip[iy];
        aLine.push_back(ip);
        aLine.push_back(ip0);
      }
    }
  }
  
  nav.camera.view_height = 1.0;
  nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  setSomeLighting();
  // Enter main loop
  ::glutMainLoop();
  return 0;
}
