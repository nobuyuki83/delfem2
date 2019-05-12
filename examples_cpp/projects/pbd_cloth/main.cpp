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

#include "delfem2/v23m3q.h"
#include "delfem2/dyntri_v2.h"
#include "delfem2/mshtopo.h"

#include "delfem2/funcs_gl.h"
#include "delfem2/funcs_glut.h"


static void FetchData
(double* val_to,
 int nno, int ndim,
 const int* aIP,
 const double* val_from,
 int nstride, int noffset)
{
  assert( nstride >= ndim );
  for(int ino=0;ino<nno;++ino){
    int ip = aIP[ino];
    for(int idim=0;idim<ndim;++idim){
      val_to[ino*ndim+idim] = val_from[ip*nstride+noffset+idim];
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////

std::vector<CEPo2> aPo2D;
std::vector<CVector2> aVec2;
std::vector<ETri> aETri;
//std::vector<double> aXY0; // undeformed vertex positions
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aXYZt;
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
double mass_point = 0.01;

CGlutWindowManager win;
bool is_animation = false;

//////////////////////////////////

void StepTime()
{
  const double dt = 0.01;
  const double gravity[3] = {0.0, 0.0, -10.0};
  const int ndof = aXYZ.size();
  const int np = ndof/3;
  for(int ip=0;ip<np;++ip){
    aUVW[ip*3+0] += dt*gravity[0];
    aUVW[ip*3+1] += dt*gravity[1];
    aUVW[ip*3+2] += dt*gravity[2];
  }
  for(int idof=0;idof<ndof;++idof){
    if( aBCFlag[idof] != 0 ){
      aXYZt[idof] = aXYZ[idof];
    }
    else{
      aXYZt[idof] = aXYZ[idof]+dt*aUVW[idof];
    }
  }
  for(int it=0;it<aETri.size();++it){
    const int aIP[3] = {aETri[it].v[0],aETri[it].v[1],aETri[it].v[2]};
    double P[3][2] = {
      {aVec2[aIP[0]].x,aVec2[aIP[0]].y},
      {aVec2[aIP[1]].x,aVec2[aIP[1]].y},
      {aVec2[aIP[2]].x,aVec2[aIP[2]].y} };
    double p[3][3]; FetchData(&p[0][0], 3, 3, aIP, aXYZt.data(), 3, 0);
    double C[3], dCdp[3][9];  ConstraintProjection_CST(C, dCdp, P, p);
//    Check(P,p);
//    break;
    const double m = 1.0;
    const double M[9] = {m,0,0, 0,m,0, 0,0,m};
    double Minv[9]; InverseMat3(Minv, M);
    double MinvC[3][9];
    for(int idim=0;idim<3;++idim){
      const double dC0dpi[3] = {dCdp[0][0*3+idim],dCdp[0][1*3+idim],dCdp[0][2*3+idim]};
      const double dC1dpi[3] = {dCdp[1][0*3+idim],dCdp[1][1*3+idim],dCdp[1][2*3+idim]};
      const double dC2dpi[3] = {dCdp[2][0*3+idim],dCdp[2][1*3+idim],dCdp[2][2*3+idim]};
      double y0[3]; MatVec3(y0, Minv, dC0dpi);
      double y1[3]; MatVec3(y1, Minv, dC1dpi);
      double y2[3]; MatVec3(y2, Minv, dC2dpi);
      MinvC[0][0*3+idim] = y0[0];
      MinvC[0][1*3+idim] = y0[1];
      MinvC[0][2*3+idim] = y0[2];
      MinvC[1][0*3+idim] = y1[0];
      MinvC[1][1*3+idim] = y1[1];
      MinvC[1][2*3+idim] = y1[2];
      MinvC[2][0*3+idim] = y2[0];
      MinvC[2][1*3+idim] = y2[1];
      MinvC[2][2*3+idim] = y2[2];
    }
    double A[9] = {0,0,0, 0,0,0, 0,0,0};
    for(int i=0;i<9;++i){
      A[0*3+0] += MinvC[0][i]*dCdp[0][i];
      A[0*3+1] += MinvC[0][i]*dCdp[1][i];
      A[0*3+2] += MinvC[0][i]*dCdp[2][i];
      A[1*3+0] += MinvC[1][i]*dCdp[0][i];
      A[1*3+1] += MinvC[1][i]*dCdp[1][i];
      A[1*3+2] += MinvC[1][i]*dCdp[2][i];
      A[2*3+0] += MinvC[2][i]*dCdp[0][i];
      A[2*3+1] += MinvC[2][i]*dCdp[1][i];
      A[2*3+2] += MinvC[2][i]*dCdp[2][i];
    }
    double Ainv[9]; InverseMat3(Ainv, A);
    double lmd[3]; MatVec3(lmd, Ainv, C);
    for(int ine=0;ine<3;++ine){
      aXYZt[aIP[ine]*3+0] -= MinvC[0][ine*3+0]*lmd[0] + MinvC[1][ine*3+0]*lmd[1] + MinvC[2][ine*3+0]*lmd[2];
      aXYZt[aIP[ine]*3+1] -= MinvC[0][ine*3+1]*lmd[0] + MinvC[1][ine*3+1]*lmd[1] + MinvC[2][ine*3+1]*lmd[2];
      aXYZt[aIP[ine]*3+2] -= MinvC[0][ine*3+2]*lmd[0] + MinvC[1][ine*3+2]*lmd[1] + MinvC[2][ine*3+2]*lmd[2];
    }
  }
  for(int idof=0;idof<ndof;++idof){
    if( aBCFlag[idof] != 0 ) continue;
    aUVW[idof] = (aXYZt[idof]-aXYZ[idof])/dt;
  }
  for(int idof=0;idof<ndof;++idof){
    if( aBCFlag[idof] != 0 ) continue;
    aXYZ[idof] = aXYZt[idof];
  }

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
    case 'b':
    {

    }
    default:
      break;
  }
}

void myGlutSpecial(int key, int x, int y){
  win.glutSpecial(key,x,y);
}

void GenMesh(const std::vector< std::vector<double> >& aaXY)
{
  std::vector<int> loopIP_ind, loopIP;
  const double elen = 0.07;
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
    MeshingInside(aPo2D,aETri,aVec2, loopIP,
                  elen, param);
  }
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
  
  double lenx = 1.0;
  {
    std::vector< std::vector<double> > aaXY;
    aaXY.resize(1);
    double xys[8] = {-0.5,-0.5, +0.5,-0.5, +0.5,+0.5, -0.5,+0.5};
    aaXY[0].assign(xys,xys+8);
    GenMesh(aaXY);
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
  aBCFlag.resize(np*3,0);
  for(int ip=0;ip<np;++ip){
    if( aXYZ[ip*3+0]  < -0.49*lenx ){
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aBCFlag[ip*3+2] = 1;
    }
  }
  /*
  {
    CMatrix3 m;
    m.SetRotMatrix_Cartesian(0.0, 0.0, 1.0);
    const int np = aXYZ.size()/3;
    for(int ip=0;ip<np;++ip){
      double p0[3] = {aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]};
      double p1[3];  m.MatVec(p0,p1);
      aXYZ[ip*3+0] = p1[0]*1.0;
      aXYZ[ip*3+1] = p1[1]*1.0;
      aXYZ[ip*3+2] = p1[2]*1.0;
    }
  }
   */
  
  win.camera.view_height = 1.0;
  win.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  setSomeLighting();
  // Enter main loop
  ::glutMainLoop();
  return 0;
}
