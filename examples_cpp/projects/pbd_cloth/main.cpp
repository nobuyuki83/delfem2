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

#include "delfem2/objfunc_v23.h"
#include "delfem2/dyntri_v2.h"
#include "delfem2/mshtopo.h"
#include "delfem2/cad_dyntri_v23_gl.h"

#include "delfem2/funcs_gl.h"
#include "delfem2/funcs_glut.h"


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

void GenMesh
(std::vector<CEPo2>& aPo2D,
 std::vector<ETri>& aETri,
 std::vector<CVector2>& aVec2,
 const std::vector< std::vector<double> >& aaXY,
 double resolution_edge,
 double resolution_face)
{
  std::vector<int> loopIP_ind, loopIP;
  {
    JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                         aaXY);
    if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
      return;
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind,aVec2);
    if( resolution_edge > 10e-10 ){
      ResamplingLoop(loopIP_ind,loopIP,aVec2,
                     resolution_edge );
    }
  }
  ////
  Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                 loopIP_ind,loopIP);
  if( resolution_face > 1.0e-10 ){
    CInputTriangulation_Uniform param(1.0);
    std::vector<int> flg_pnt(aVec2.size());
    std::vector<int> flg_tri(aETri.size());
    MeshingInside(aPo2D,aETri,aVec2, flg_pnt,flg_tri,
                  aVec2.size(), 0, resolution_face, param);
  }
}



////////////////////////////////////////////////////////////////////////////////////

std::vector<CEPo2> aPo2D;
std::vector<CVector2> aVec2;
std::vector<ETri> aETri;
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
  PBD_Pre3D(aXYZt,
            dt, gravity, aXYZ, aUVW, aBCFlag);
  for(int it=0;it<aETri.size();++it){
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
  for(int it=0;it<aETri.size();++it){
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
