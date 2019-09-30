
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <stack>
#include <algorithm>
#include <cstdlib>
#include <math.h>
#include <time.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/msh.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vec2.h"
#include "delfem2/mat3.h"
#include "delfem2/mats.h"
#include "delfem2/dtri.h"

#include "delfem2/v23m3q.h"
#include "delfem2/ilu_mats.h"
#include "delfem2/fem_emats.h"
#include "delfem2/dtri_v2.h"

#include "delfem2/gl_color.h"
#include "delfem2/gl_v23.h"
#include "delfem2/gl2_funcs.h"

#include "../glut_funcs.h"




void GenMesh
(std::vector<CVector2>& aVec2,
 std::vector<CEPo2>& aPo2D,
 std::vector<ETri>& aETri,
 double elen,
 const std::vector< std::vector<double> >& aaXY)
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
    std::vector<int> aFlgPnt(aVec2.size());
    std::vector<int> aFlgTri(aETri.size());
    MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
                  aVec2.size(), 0, elen, param);
  }
}

void RotationAtMeshPoints
(std::vector<double>& aR,
 const std::vector<double>& aXYZ,
 const std::vector<double>& aDisp,
 const std::vector<int>& psup_ind,
 const std::vector<int>& psup)
{
  const unsigned int np = aXYZ.size()/3;
  aR.resize(np*9);
  for(int ip=0;ip<aXYZ.size()/3;++ip){
    CVector3 Pi(aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]);
    CVector3 pi(aXYZ[ip*3+0]+aDisp[ip*3+0],
                aXYZ[ip*3+1]+aDisp[ip*3+1],
                aXYZ[ip*3+2]+aDisp[ip*3+2]);
    CMatrix3 A;
    A.SetZero();
    for(int jjp=psup_ind[ip];jjp<psup_ind[ip+1];++jjp){
      int jp = psup[jjp];
      CVector3 Pj(aXYZ[jp*3+0],aXYZ[jp*3+1],aXYZ[jp*3+2]);
      CVector3 pj(aXYZ[jp*3+0]+aDisp[jp*3+0],
                  aXYZ[jp*3+1]+aDisp[jp*3+1],
                  aXYZ[jp*3+2]+aDisp[jp*3+2]);
      A += Mat3_OuterProduct(pj-pi,Pj-Pi);
    }
    GetRotPolarDecomp(aR.data()+ip*9,
                      A.mat, 100);
  }
}


/////////////////////////////////////////

CNav3D_GLUT nav;
bool is_animatio = false;
bool is_stiffness_warping = true;

std::vector<unsigned int> aTet;
std::vector<double> aXYZ;
std::vector<double> aDisp;
std::vector<double> aVelo;
std::vector<int> aBCFlag;
double dt = 0.03;
double myu = 200.0;
double lambda = 1.0;
double rho = 1.0;
const double gravity[3] = {0.0, -4.0, 0.0};

CMatrixSparse<double> mat_A;
std::vector<double> vec_b;
CPreconditionerILU<double>  ilu_A;
std::vector<int> psup_ind, psup;
std::vector<double> aR;

/////////////////////////////////////////////////////////////////////////////


void InitializeProblem_ShellEigenPB()
{
  const int np = (int)aXYZ.size()/3;
  JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                              aTet.data(), aTet.size()/4, 4,
                                              (int)aXYZ.size()/3);
  JArray_Sort(psup_ind, psup);
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  ilu_A.Initialize_ILU0(mat_A);
}
  
  ////////////////////////////////////////////


void Solve_Linear()
{
  mat_A.SetZero();
  vec_b.assign(aXYZ.size(),0.0);
  MergeLinSys_SolidLinear_BEuler_MeshTet3D(mat_A, vec_b.data(),
                                           myu, lambda,
                                           rho, gravity,
                                           dt,
                                           aXYZ.data(), aXYZ.size()/3,
                                           aTet.data(), aTet.size()/4,
                                           aDisp.data(),
                                           aVelo.data());
  mat_A.SetBoundaryCondition(aBCFlag.data(),aXYZ.size()/3,3);
  setRHS_Zero(vec_b,aBCFlag,0);
  ////
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  const int nDoF = aXYZ.size();
  std::vector<double> dv(nDoF,0.0);
  std::vector<double> aConv = Solve_PBiCGStab(vec_b.data(), dv.data(),
                                              1.0e-4, 1000, mat_A, ilu_A);
  ////
  XPlusAYBZ(aDisp,nDoF,aBCFlag,
            dt, dv,
            dt, aVelo);
  XPlusAY(aVelo,nDoF,aBCFlag,
          1.0, dv);

  std::cout << "conv; " << aConv.size() <<  std::endl;
}


void Solve_StiffnessWarping()
{
  RotationAtMeshPoints(aR,
                       aXYZ,aDisp,psup_ind,psup);
  /////
  mat_A.SetZero();
  vec_b.assign(aXYZ.size(),0.0);
  MergeLinSys_SolidStiffwarp_BEuler_MeshTet3D(mat_A, vec_b.data(),
                                              myu, lambda,
                                              rho, gravity,
                                              dt,
                                              aXYZ.data(), aXYZ.size()/3,
                                              aTet.data(), aTet.size()/4,
                                              aDisp.data(),
                                              aVelo.data(),
                                              aR);
  ////
  mat_A.SetBoundaryCondition(aBCFlag.data(),aXYZ.size()/3,3);
  setRHS_Zero(vec_b,aBCFlag,0);
  /////////////
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  const int nDoF = aXYZ.size();
  std::vector<double> dv(nDoF,0.0);
  std::vector<double> aConv = Solve_PBiCGStab(vec_b.data(), dv.data(),
                                              1.0e-4, 1000, mat_A, ilu_A);
  ////
  XPlusAYBZ(aDisp,nDoF,aBCFlag,
            dt, dv,
            dt, aVelo);
  XPlusAY(aVelo,nDoF,aBCFlag,
          1.0,dv);
  std::cout << "conv; " << aConv.size() <<  std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void myGlutDisplay(void)
{
  if( is_stiffness_warping ){
    ::glClearColor(0.4f, 0.9f, 0.9f ,1.0f);
  }
  else{
    ::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  }
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.0f, 1.0f );

  nav.SetGL_Camera();

//  glEnable(GL_BLEND);
//  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  {
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
//    glShadeModel(GL_SMOOTH);
    glShadeModel(GL_FLAT);
 }
  ::glDisable(GL_LIGHTING);
  
  { // defomred edge
    ::glColor3d(0,0,0);
    DrawMeshTet3D_EdgeDisp(aXYZ.data(),
                           aTet.data(),aTet.size()/4,
                           aDisp.data(),
                           1.0);
  }
  
  ::glDisable(GL_LIGHTING);
  for(int ip=0;ip<aXYZ.size()/3;++ip){
    CVector3 pi(aXYZ[ip*3+0]+aDisp[ip*3+0],
                aXYZ[ip*3+1]+aDisp[ip*3+1],
                aXYZ[ip*3+2]+aDisp[ip*3+2]);
    CVector3 ex(aR[ip*9+0],aR[ip*9+3],aR[ip*9+6]);
    CVector3 ey(aR[ip*9+1],aR[ip*9+4],aR[ip*9+7]);
    CVector3 ez(aR[ip*9+2],aR[ip*9+5],aR[ip*9+8]);
    ::glBegin(GL_LINES);
    ::glColor3d(1,0,0);
    ::myGlVertex(pi);
    ::myGlVertex(pi+0.04*ex);
    ::glColor3d(0,1,0);
    ::myGlVertex(pi);
    ::myGlVertex(pi+0.04*ey);
    ::glColor3d(0,0,1);
    ::myGlVertex(pi);
    ::myGlVertex(pi+0.04*ez);
    ::glEnd();
  }
    
  {
    ::glEnable(GL_LIGHTING);
    {
      float color[4] = {180.0/256.0, 180.0/256.0, 130.0/256.0,1.0f};
      ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
      ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
      glShadeModel(GL_FLAT);
    }
    DrawMeshTet3D_FaceNorm(aXYZ.data(), 
                           aTet.data(), aTet.size()/4);

  }

  
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }

  ::glColor3d(0,0,0);
  ShowFPS();
  ::glutSwapBuffers();
}

void myGlutResize(int w, int h)
{
  ::glViewport(0, 0, w, h);
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
  switch (Key)
  {
    case 'q':
    case 'Q':
    case '\033':
      exit(0);  /* '\033' ? ESC ? ASCII ??? */
    case 'a':
    {
      is_animatio = !is_animatio;
      break;
    }
    case 'l':
    {
      is_stiffness_warping = !is_stiffness_warping;
      break;
    }
    case 's':
    {
      Solve_StiffnessWarping();
      break;
    }
	::glutPostRedisplay();
  }
}

void myGlutIdle(){
  if( is_animatio ){
    if( is_stiffness_warping ){ Solve_StiffnessWarping(); }
    else{                       Solve_Linear();           }
  }
  ::glutPostRedisplay();
}


void myGlutSpecial(int Key, int x, int y)
{
  nav.glutSpecial(Key, x, y);
  ::glutPostRedisplay();
}


int main(int argc,char* argv[])
{
	// Initialize GLUT
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
 	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("FEM View");

	// Define callback functions
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
	glutIdleFunc(myGlutIdle);
  
  
  {
    std::vector< std::vector<double> > aaXY;
    {
      const double aXY[8] = {
        -1,-0.1,
        +1,-0.1,
        +1,+0.1,
        -1,+0.1 };
      aaXY.push_back( std::vector<double>(aXY,aXY+8) );
    }
    std::vector<CVector2> aVec2;
    std::vector<CEPo2> aPo2D;
    std::vector<ETri> aETri;
    GenMesh(aVec2,aPo2D,aETri,
            0.075, aaXY);
    std::vector<double> aXY;
    std::vector<unsigned int> aTri;
    CMeshTri2D(aXY,aTri,
               aVec2,aETri);
    ExtrudeTri2Tet(3, 0.075,
                   aXYZ,aTet,
                   aXY,aTri);
  }
  aDisp.assign(aXYZ.size(), 0.0);
  aVelo.assign(aXYZ.size(), 0.0);
  aBCFlag.assign(aXYZ.size(),0);
  for(int ip=0;ip<aXYZ.size()/3;++ip){
    double x0 = aXYZ[ip*3+0];
    if( fabs(x0+1)<1.0e-10 ){
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aBCFlag[ip*3+2] = 1;
    }
  }
  InitializeProblem_ShellEigenPB();
  RotationAtMeshPoints(aR,
                       aXYZ,aDisp,psup_ind,psup);

  
  nav.camera.view_height = 2.0;
  nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  setSomeLighting();
  glutMainLoop();
	return 0;
}
