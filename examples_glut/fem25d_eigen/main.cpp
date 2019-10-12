
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
#include "delfem2/mats.h"
#include "delfem2/dtri.h"

#include "delfem2/ilu_mats.h"
#include "delfem2/fem_emats.h"
#include "delfem2/dtri_v2.h"

#include "delfem2/gl2_color.h"
#include "delfem2/gl_v23.h"
#include "delfem2/gl2_funcs.h"

#include "../glut_cam.h"


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

void SetValue_SolidEigen3D_MassLumpedSqrtInv_KernelModes6
(double* aMassLumpedSqrtInv,
 double* aModesKer,
 const double* aXYZ, int nXYZ,
 const unsigned int* aTet, int nTet)
{
  const int nDoF = nXYZ*3;
  std::vector<double> aMassLumpedSqrt(nXYZ);
  MassPoint_Tet3D(aMassLumpedSqrt.data(),
                   1, aXYZ, nXYZ, aTet,nTet);
  
  for(int ip=0;ip<nXYZ;++ip){
    aMassLumpedSqrt[ip] = sqrt(aMassLumpedSqrt[ip]);
  }
  
  {
    for(int i=0;i<nDoF*6;++i){ aModesKer[i] = 0.0; }
    double* p0 = aModesKer+nDoF*0;
    double* p1 = aModesKer+nDoF*1;
    double* p2 = aModesKer+nDoF*2;
    double* p3 = aModesKer+nDoF*3;
    double* p4 = aModesKer+nDoF*4;
    double* p5 = aModesKer+nDoF*5;
    for(int ip=0;ip<nXYZ;++ip){
      const double x0 = aXYZ[ip*3+0];
      const double y0 = aXYZ[ip*3+1];
      const double z0 = aXYZ[ip*3+2];
      const double m0 = aMassLumpedSqrt[ip];
      p0[ip*3+0] = m0;
      p1[ip*3+1] = m0;
      p2[ip*3+2] = m0;
      p3[ip*3+2] = -y0*m0;  p3[ip*3+1] = +z0*m0;
      p4[ip*3+0] = -z0*m0;  p4[ip*3+2] = +x0*m0;
      p5[ip*3+1] = -x0*m0;  p5[ip*3+0] = +y0*m0;
    }
    NormalizeX(p0,nDoF);
    OrthogonalizeToUnitVectorX(p1, p0, nDoF);
    OrthogonalizeToUnitVectorX(p2, p0, nDoF);
    OrthogonalizeToUnitVectorX(p3, p0, nDoF);
    OrthogonalizeToUnitVectorX(p4, p0, nDoF);
    OrthogonalizeToUnitVectorX(p5, p0, nDoF);
    NormalizeX(p1,nDoF);
    OrthogonalizeToUnitVectorX(p2, p1, nDoF);
    OrthogonalizeToUnitVectorX(p3, p1, nDoF);
    OrthogonalizeToUnitVectorX(p4, p1, nDoF);
    OrthogonalizeToUnitVectorX(p5, p1, nDoF);
    NormalizeX(p2,nDoF);
    OrthogonalizeToUnitVectorX(p3, p2, nDoF);
    OrthogonalizeToUnitVectorX(p4, p2, nDoF);
    OrthogonalizeToUnitVectorX(p5, p2, nDoF);
    NormalizeX(p3,nDoF);
    OrthogonalizeToUnitVectorX(p4, p3, nDoF);
    OrthogonalizeToUnitVectorX(p5, p3, nDoF);
    NormalizeX(p4,nDoF);
    OrthogonalizeToUnitVectorX(p5, p4, nDoF);
    NormalizeX(p5,nDoF);
  }
  
  for(int ip=0;ip<nXYZ;++ip){
    aMassLumpedSqrtInv[ip] = 1.0/aMassLumpedSqrt[ip];
  }
}


/////////////////////////////////////////

CNav3D_GLUT nav;

std::vector<unsigned int> aTet;
std::vector<double> aXYZ;
std::vector<double> aMassLumpedSqrtInv;
std::vector<double> aTmp0;
std::vector<double> aTmp1;
std::vector<double> aMode;
std::vector<double> aModesKer;
double lamda0 = 0.1;

CMatrixSparse<double> mat_A;
CPreconditionerILU<double>  ilu_A;

/////////////////////////////////////////////////////////////////////////////

void RemoveKernel()
{
  const int nDoF = aXYZ.size();
  const double* p0 = aModesKer.data()+nDoF*0;
  const double* p1 = aModesKer.data()+nDoF*1;
  const double* p2 = aModesKer.data()+nDoF*2;
  const double* p3 = aModesKer.data()+nDoF*3;
  const double* p4 = aModesKer.data()+nDoF*4;
  const double* p5 = aModesKer.data()+nDoF*5;
  double* p = aTmp0.data();
  OrthogonalizeToUnitVectorX(p, p0, nDoF);
  OrthogonalizeToUnitVectorX(p, p1, nDoF);
  OrthogonalizeToUnitVectorX(p, p2, nDoF);
  OrthogonalizeToUnitVectorX(p, p3, nDoF);
  OrthogonalizeToUnitVectorX(p, p4, nDoF);
  OrthogonalizeToUnitVectorX(p, p5, nDoF);
  NormalizeX(p, nDoF);
}




void InitializeProblem_ShellEigenPB()
{
  const int np = (int)aXYZ.size()/3;
  const int nDoF = np*3;
  aTmp0.assign(nDoF, 0.0);
  //////
  std::vector<int> psup_ind, psup;
  JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                              aTet.data(), aTet.size()/4, 4,
                                              (int)aXYZ.size()/3);
  JArray_Sort(psup_ind, psup);
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  ilu_A.Initialize_ILU0(mat_A);
  
  ///////////////////////////////////////////
  aMassLumpedSqrtInv.resize(np);
  aModesKer.resize(nDoF*6);
  SetValue_SolidEigen3D_MassLumpedSqrtInv_KernelModes6(aMassLumpedSqrtInv.data(),
                                          aModesKer.data(),
                                          aXYZ.data(), aXYZ.size()/3,
                                          aTet.data(), aTet.size()/4);
  
  ////////////////////////////////////////////
  
  double myu = 1.0;
  double lambda = 1.0;
  double rho = 1.0;
  mat_A.SetZero();
  aMode.assign(nDoF, 0.0);
  aTmp0.assign(nDoF, 0.0);
  double gravity[3] = {0,0,0};
  MergeLinSys_SolidLinear_Static_MeshTet3D(mat_A, aMode.data(),
                                           myu, lambda, rho, gravity,
                                           aXYZ.data(), aXYZ.size()/3,
                                           aTet.data(), aTet.size()/4,
                                           aTmp0.data());
  MatSparse_ScaleBlk_LeftRight(mat_A,
                               aMassLumpedSqrtInv.data());
  mat_A.AddDia(0.8);
  
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
}

void Solve(){
  aMode.assign(aTmp0.size(),0.0);
  const double conv_ratio = 1.0e-5;
  const int iteration = 1000;
  std::vector<double> aConv;
  aTmp1 = aTmp0;
  aConv = Solve_PCG(aTmp1.data(), aMode.data(),
                    conv_ratio, iteration, mat_A, ilu_A);
  double lam0 = DotX(aTmp0.data(), aMode.data(), aTmp0.size());
  std::cout << 1.0/lam0 << std::endl;
  aTmp0 = aMode;
  ////
  RemoveKernel();

  ////
  for(int ip=0;ip<aTmp0.size()/3;++ip){
    const double s0 = aMassLumpedSqrtInv[ip];
    for(int idim=0;idim<3;++idim){
      aMode[ip*3+idim] = aTmp0[ip*3+idim]*s0;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void myGlutDisplay(void)
{
//	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
	::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 3.1f, 2.0f );

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
  
  {
    ::glColor3d(0,0,0);
    DrawMeshTet3D_EdgeDisp(aXYZ.data(),
                           aTet.data(),aTet.size()/4,
                           aMode.data(),
                          0.1);
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
    case 's':
    {
      Solve();
      break;
    }
	::glutPostRedisplay();
  }
}



void myGlutIdle(){
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
        -1,-0.3,
        +1,-0.3,
        +1,+0.3,
        -1,+0.3 };
      aaXY.push_back( std::vector<double>(aXY,aXY+8) );
    }

    std::vector<CVector2> aVec2;
    std::vector<CEPo2> aPo2D;
    std::vector<ETri> aETri;
    GenMesh(aVec2,aPo2D,aETri,
            0.05, aaXY);
    std::cout << aVec2.size() << " " << aPo2D.size() << " " << aETri.size() << std::endl;
    std::vector<double> aXY;
    std::vector<unsigned int> aTri;
    CMeshTri2D(aXY,aTri,
               aVec2,aETri);
    ExtrudeTri2Tet(1, 0.1,
                   aXYZ,aTet,
                   aXY,aTri);
  }
  
  
  aTmp0.assign(aXYZ.size(),0.0);
  aMode.assign(aXYZ.size(),0.0);
  
  InitializeProblem_ShellEigenPB();
  
  for(int i=0;i<aXYZ.size();++i){
    aTmp0[i] = (rand()+1.0)/(RAND_MAX+1.0);
  }
  RemoveKernel();
  
  
  nav.camera.view_height = 2.0;
  nav.camera.camera_rot_mode = CAMERA_ROT_TBALL;
  
  setSomeLighting();
  glutMainLoop();
	return 0;
}
