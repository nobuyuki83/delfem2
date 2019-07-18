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

#include "delfem2/dyntri_v2.h"
#include "delfem2/mshtopo.h"
#include "delfem2/matrix_sparse.h"
#include "delfem2/ilu_sparse.h"
#include "delfem2/fem_ematrix.h"
#include "delfem2/cloth_internal.h"

#include "delfem2/gl_funcs.h"
#include "delfem2/glut_funcs.h"


class CInput_ContactNothing: public CInput_Contact
{
public:
  double penetrationNormal(double& nx, double &ny, double& nz,
                           double px, double py, double pz) const
  {
    return -100;
  }
};


////////////////////////////////////////////////////////////////////////////////////

std::vector<CEPo2> aPo2D;
std::vector<CVector2> aVec2;
std::vector<ETri> aETri;
std::vector<double> aXYZ0; // undeformed vertex positions
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
std::vector<unsigned int> aTri;  // index of triangles
std::vector<unsigned int> aQuad; // index of 4 vertices required for bending
CMatrixSparse mat_A; // coefficient matrix
CPreconditionerILU  ilu_A; // ilu decomposition of the coefficient matrix
double mass_point = 0.01;

int idp_nearest = -1;
int press_button = -1;
double mov_begin_x, mov_begin_y;
bool is_animation = true;
double mag = 1.0;

//////////////////////////////////

void GenMesh(const std::vector< std::vector<double> >& aaXY)
{
  std::vector<int> loopIP_ind, loopIP;
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
                  aVec2.size(),0, elen, param);
  }
}

void StepTime()
{
  const double lambda = 1.0; // Lame's 1st parameter
  const double myu    = 4.0; // Lame's 2nd parameter
  const double stiff_bend = 0.0e-3; // bending stiffness
  const double areal_density = 1.0; // areal density of a cloth
  const double gravity[3] = {0,-1,0}; // gravitatinal accereration
  double time_step_size = 0.03; // size of time step
  const double stiff_contact = 1.0e+3;
  const double contact_clearance = 0.02;
  ////
  CInput_ContactNothing c1;
  StepTime_InternalDynamicsILU(aXYZ, aUVW, mat_A, ilu_A,
                               aXYZ0, aBCFlag,
                               aTri, aQuad,
                               time_step_size,
                               lambda, myu, stiff_bend,
                               gravity, mass_point,
                               stiff_contact,contact_clearance,c1);
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
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1, 4.0 );
  
  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
  
  ::glPointSize(5);
  ::glLineWidth(1);
  ::glPointSize(5);
  ::glColor3d(1,1,0);
    //////////////
  {
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0.8, 0.8, 0.8);
    /*
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
    ::glEnable(GL_DEPTH_TEST);
     */
    DrawMeshTri3D_FaceNorm(aXYZ, aTri);
  }
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  DrawMeshTri3D_Edge(aXYZ, aTri);

  
  ShowFPS();
  
  glutSwapBuffers();
}


void myGlutIdle(){
  if( is_animation ){
    StepTime();
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
//    Refine(mov_begin_x, mov_begin_y);
  }
  if( button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN ){
//    Coarse(mov_begin_x, mov_begin_y);
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
  
  double lenx = 1.0;
  {
    std::vector< std::vector<double> > aaXY;
    aaXY.resize(1);
    double xys[8] = {-0.5,-0.5, +0.5,-0.5, +0.5,+0.5, -0.5,+0.5};
    aaXY[0].assign(xys,xys+8);
    GenMesh(aaXY);
  }
  /////////
  aXYZ0.resize(aPo2D.size()*3);
  for(int ip=0;ip<aPo2D.size();++ip){
    aXYZ0[ip*3+0] = aVec2[ip].x;
    aXYZ0[ip*3+1] = aVec2[ip].y;
    aXYZ0[ip*3+2] = 0.0;
  }
  aXYZ = aXYZ0;
  /////
  aTri.resize(aETri.size()*3);
  for(int it=0;it<aETri.size();++it){
    aTri[it*3+0] = aETri[it].v[0];
    aTri[it*3+1] = aETri[it].v[1];
    aTri[it*3+2] = aETri[it].v[2];
  }
  ElemQuad_DihedralTri(aQuad,aTri.data(),aTri.size()/3,aXYZ0.size()/3);
  {
    const int np = aXYZ0.size()/3;
    mat_A.Initialize(np,3,true);
    std::vector<int> psup_ind,psup;
    JArray_MeshOneRingNeighborhood(psup_ind, psup,
                                   aQuad.data(),aQuad.size()/4, 4, np);
    JArray_Sort(psup_ind, psup);
    mat_A.SetPattern(psup_ind.data(),psup_ind.size(), psup.data(),psup.size());
    ilu_A.Initialize_ILU0(mat_A);
  }
  aUVW.resize(aXYZ.size(),0.0);
  aBCFlag.resize(aXYZ.size(),0);
  for(int ip=0;ip<aXYZ.size()/3;++ip){
    if( aXYZ[ip*3+0]  < -0.49*lenx ){
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aBCFlag[ip*3+2] = 1;
    }
  }
  
  setSomeLighting();
  // Enter main loop
  ::glutMainLoop();
  return 0;
}
