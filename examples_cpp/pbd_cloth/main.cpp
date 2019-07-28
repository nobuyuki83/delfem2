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

#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mshtopo.h"
#include "delfem2/dyntri.h"

#include "delfem2/objfunc_v23.h"
#include "delfem2/dyntri_v2.h"

#include "delfem2/gl_cad_dyntri_v23.h"
#include "delfem2/gl_funcs.h"
#include "delfem2/glut_funcs.h"


static void FetchData
(double* val_to,
 int nno, int ndim,
 const int* aIP,
 const double* val_from,
 int nstride)
{
  assert( nstride >= ndim );
  for(int ino=0;ino<nno;++ino){
    int ip = aIP[ino];
    for(int idim=0;idim<ndim;++idim){
      val_to[ino*ndim+idim] = val_from[ip*nstride+idim];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////

std::vector<CEPo2> aPo2D;
std::vector<CVector2> aVec2;
std::vector<ETri> aETri;
std::vector<int> aLine;
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aXYZt;
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
const double mass_point = 0.01;
const double dt = 0.01;
const double gravity[3] = {0.0, 0.0, -10.0};

CGlutWindowManager win;
bool is_animation = false;

//////////////////////////////////

void StepTime()
{
  PBD_Pre3D(aXYZt,
            dt, gravity, aXYZ, aUVW, aBCFlag);
  for(unsigned int it=0;it<aETri.size();++it){
    const int aIP[3] = {aETri[it].v[0],aETri[it].v[1],aETri[it].v[2]};
    const double P[3][2] = {
      {aVec2[aIP[0]].x,aVec2[aIP[0]].y},
      {aVec2[aIP[1]].x,aVec2[aIP[1]].y},
      {aVec2[aIP[2]].x,aVec2[aIP[2]].y} };
    double p[3][3]; FetchData(&p[0][0], 3, 3, aIP, aXYZt.data(), 3);
    double C[3], dCdp[3][9];  PBD_CdC_TriStrain2D3D(C, dCdp, P, p);
    double m[3] = {1,1,1};
    PBD_Update_Const3(aXYZt.data(), 3, 3, m, C, &dCdp[0][0], aIP);
  }
  for(unsigned int it=0;it<aETri.size();++it){
    for(int ie=0;ie<3;++ie){
      const int jt0 = aETri[it].s2[ie];
      if( jt0 == -1 ){ continue; }
      if( jt0 > it ){ continue; }
      const int rt0 = aETri[it].r2[ie];
      const int je0 = (6-rt0-ie)%3;
      assert( aETri[jt0].s2[je0] == it);
      const int aIP[4] = {aETri[it].v[ie],aETri[jt0].v[je0],aETri[it].v[(ie+1)%3],aETri[it].v[(ie+2)%3]};
      const double P[4][3] = {
        {aVec2[aIP[0]].x,aVec2[aIP[0]].y, 0.0},
        {aVec2[aIP[1]].x,aVec2[aIP[1]].y, 0.0},
        {aVec2[aIP[2]].x,aVec2[aIP[2]].y, 0.0},
        {aVec2[aIP[3]].x,aVec2[aIP[3]].y, 0.0} };
      double p[4][3]; FetchData(&p[0][0], 4, 3, aIP, aXYZt.data(), 3);
      double C[3], dCdp[3][12];
      PBD_CdC_QuadBend(C, dCdp,
                       P, p);
      double m[4] = {1,1,1,1};
      PBD_Update_Const3(aXYZt.data(), 4,3, m, C, &dCdp[0][0], aIP);
    }
  }
  for(unsigned int il=0;il<aLine.size()/2;++il){
    int ip0 = aLine[il*2+0];
    int ip1 = aLine[il*2+1];
    const double p[2][3] = {
      {aXYZt[ip0*3+0], aXYZt[ip0*3+1], aXYZt[ip0*3+2]},
      {aXYZt[ip1*3+0], aXYZt[ip1*3+1], aXYZt[ip1*3+2]} };
    double d0 = Distance3D(p[0], p[1]);
    double dLen = 0.01;
    if( d0 > dLen ){
      double n01[3] = {p[1][0]-p[0][0], p[1][1]-p[0][1], p[1][2]-p[0][2]};
      double l01 = Length3D(n01);
      double invl01 = 1.0/l01;
      n01[0] *= invl01;
      n01[1] *= invl01;
      n01[2] *= invl01;
      aXYZt[ip0*3+0] += n01[0]*dLen*0.5;
      aXYZt[ip0*3+1] += n01[1]*dLen*0.5;
      aXYZt[ip0*3+2] += n01[2]*dLen*0.5;
      aXYZt[ip1*3+0] -= n01[0]*dLen*0.5;
      aXYZt[ip1*3+1] -= n01[1]*dLen*0.5;
      aXYZt[ip1*3+2] -= n01[2]*dLen*0.5;
    }
    else{
      aXYZt[ip0*3+0] = (p[0][0]+p[1][0])*0.5;
      aXYZt[ip0*3+1] = (p[0][1]+p[1][1])*0.5;
      aXYZt[ip0*3+2] = (p[0][2]+p[1][2])*0.5;
      aXYZt[ip1*3+0] = (p[0][0]+p[1][0])*0.5;
      aXYZt[ip1*3+1] = (p[0][1]+p[1][1])*0.5;
      aXYZt[ip1*3+2] = (p[0][2]+p[1][2])*0.5;
    }
  }
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
  
  win.SetGL_Camera();
  
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
  win.glutMotion(x,y);
}

void myGlutMouse(int button, int state, int x, int y)
{
  win.glutMouse(button, state, x,y);
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
  win.glutSpecial(key,x,y);
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
  
  win.camera.view_height = 1.0;
  win.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  setSomeLighting();
  // Enter main loop
  ::glutMainLoop();
  return 0;
}
