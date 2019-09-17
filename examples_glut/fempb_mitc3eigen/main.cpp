#include <iostream>
#include <math.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/emat.h"
#include "delfem2/mshtopo.h"
#include "delfem2/msh.h"
#include "delfem2/mats.h"
#include "delfem2/vec2.h"

#include "delfem2/ilu_mats.h"
#include "delfem2/dtri_v2.h"
#include "delfem2/fem_emats.h"

#include "delfem2/gl_funcs.h"
#include "delfem2/gl_color.h"

#include "../glut_funcs.h"


void SetValue_ShellPBMITC3Eigen_MassLumpedSqrtInv_KernelModes3
(double* aMassLumpedSqrtInv,
 double* aModesKer,
 double rho, double thickness,
 const double* aXY, int nXY,
 const unsigned int* aTri, int nTri)
{
  const unsigned int nDoF = nXY*3;
  std::vector<double> aMassLumpedSqrt(nDoF);
  MassLumped_ShellPlateBendingMITC3(aMassLumpedSqrt.data(),
                                          rho, thickness,
                                          aXY, nXY,
                                          aTri, nTri);
  for(int i=0;i<nDoF;++i){ aMassLumpedSqrt[i] = sqrt(aMassLumpedSqrt[i]); }
  {
    for(int i=0;i<nDoF*3;++i){ aModesKer[i] = 0.0; }
    double* p0 = aModesKer+nDoF*0;
    double* p1 = aModesKer+nDoF*1;
    double* p2 = aModesKer+nDoF*2;
    for(int ip=0;ip<nXY;++ip){
      const double x0 = aXY[ip*2+0];
      const double y0 = aXY[ip*2+1];
      const double m0 = aMassLumpedSqrt[ip*3+0];
      const double m1 = aMassLumpedSqrt[ip*3+1];
      const double m2 = aMassLumpedSqrt[ip*3+2];
      p0[ip*3+0] = m0;
      p1[ip*3+0] = +y0*m0;  p1[ip*3+1] = m1;
      p2[ip*3+0] = -x0*m0;  p2[ip*3+2] = m2;
    }
    NormalizeX(p0,nDoF);
    OrthogonalizeToUnitVectorX(p1, p0, nDoF);
    OrthogonalizeToUnitVectorX(p2, p0, nDoF);
    NormalizeX(p1,nDoF);
    OrthogonalizeToUnitVectorX(p2, p1, nDoF);
    NormalizeX(p2,nDoF);
  }
  for(int i=0;i<nDoF;++i){ aMassLumpedSqrtInv[i] = 1.0/aMassLumpedSqrt[i]; }
}

void RemoveKernel(std::vector<double>& aTmp0,
                  const std::vector<double>& aModesKer)
{
  const unsigned int nDoF = aTmp0.size();
  const double* p0 = aModesKer.data()+nDoF*0;
  const double* p1 = aModesKer.data()+nDoF*1;
  const double* p2 = aModesKer.data()+nDoF*2;
  double* p = aTmp0.data();
  OrthogonalizeToUnitVectorX(p, p0, nDoF);
  OrthogonalizeToUnitVectorX(p, p1, nDoF);
  OrthogonalizeToUnitVectorX(p, p2, nDoF);
  NormalizeX(p, nDoF);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////

// display data
bool is_animation;

CGlutWindowManager window;

std::vector<unsigned int> aTri;
std::vector<double> aXY0;
std::vector<double> aMassLumpedSqrtInv;
std::vector<double> aTmp0;
std::vector<double> aTmp1;
std::vector<double> aMode;
std::vector<double> aModesKer;

CMatrixSparse<double> mat_A;
CPreconditionerILU<double> ilu_A;

const double lenx = 0.5;
const double leny = 0.03;
const double thickness = 0.004;
double EYoung = 68.3*1.0e+9;
double Poisson = 0.34;
//double Poisson = 0.0;
double myu = EYoung/(2*(1.0+Poisson));
double lambda = Poisson*EYoung/(1+Poisson)/(1-2*Poisson);
const double rho = 2700;
const double offset_dia = 0.2;
double freq_theo = 0.0;

///////////////////////////////////////////////////////////////

void MakeMesh(){
  std::vector< std::vector<double> > aaXY;
  {
    aaXY.resize(1);
    aaXY[0].push_back(-lenx*0.5); aaXY[0].push_back(-leny*0.5);
    aaXY[0].push_back(+lenx*0.5); aaXY[0].push_back(-leny*0.5);
    aaXY[0].push_back(+lenx*0.5); aaXY[0].push_back(+leny*0.5);
    aaXY[0].push_back(-lenx*0.5); aaXY[0].push_back(+leny*0.5);
  }
  //////////////////////////////
  std::vector<CEPo2> aPo2D;
  std::vector<ETri> aETri;
  std::vector<CVector2> aVec2;
  GenMesh(aPo2D, aETri, aVec2,
          aaXY, 0.01, 0.01);
  MeshTri2D_Export(aXY0,aTri,
                   aVec2,aETri);
  std::cout<<"  ntri;"<<aTri.size()/3<<"  nXY:"<<aXY0.size()/2<<std::endl;
}

void InitializeProblem_ShellEigenPB()
{
  const int np = (int)aXY0.size()/2;
  const int nDoF = np*3;
  aTmp0.assign(nDoF, 0.0);
  //////
  std::vector<int> psup_ind, psup;
  JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                              aTri.data(), aTri.size()/3, 3,
                                              (int)aXY0.size()/2);
  JArray_Sort(psup_ind, psup);
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  ilu_A.Initialize_ILU0(mat_A);
  
  ///////////////////////////////////////////
  aMassLumpedSqrtInv.resize(nDoF);
  aModesKer.resize(nDoF*3);
  SetValue_ShellPBMITC3Eigen_MassLumpedSqrtInv_KernelModes3(aMassLumpedSqrtInv.data(),
                                                            aModesKer.data(),
                                                            rho, thickness,
                                                            aXY0.data(), aXY0.size()/2,
                                                            aTri.data(), aTri.size()/3);
  
  ////////////////////////////////////////////
  

  mat_A.SetZero();
  aMode.assign(nDoF, 0.0);
  aTmp0.assign(nDoF, 0.0);
  MergeLinSys_ShellStaticPlateBendingMITC3_MeshTri2D(mat_A, aMode.data(),
                                                     thickness,lambda, myu, 0.0, 0.0,
                                                     aXY0.data(), aXY0.size()/2,
                                                     aTri.data(), aTri.size()/3,
                                                     aTmp0.data());
  MatSparse_ScaleBlkLen_LeftRight(mat_A,
                                  aMassLumpedSqrtInv.data());
  mat_A.AddDia(offset_dia);
  
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
}

void Solve(){
  aMode.assign(aTmp0.size(),0.0);
  const double conv_ratio = 1.0e-5;
  const int iteration = 10000;
  std::vector<double> aConv;
  aTmp1 = aTmp0;
  aConv = Solve_PCG(aTmp1.data(), aMode.data(),
                    conv_ratio, iteration, mat_A, ilu_A);
  {
    double lam0 = DotX(aTmp0.data(), aMode.data(), aTmp0.size());
    double freq_sim = sqrt(1.0/lam0-offset_dia)/(2*M_PI);
    std::cout << "freq theo" << freq_theo << "   freq_sim:" << freq_sim << "   " << freq_theo/freq_sim  << "     " << aConv.size() << std::endl;
  }
  aTmp0 = aMode;
  ////
  RemoveKernel(aTmp0,
               aModesKer);
  ////
  for(unsigned int i=0;i<aTmp0.size();++i){
    aMode[i] = aTmp0[i]*aMassLumpedSqrtInv[i];
  }
}

//////////////////////////////////////////////////////////////



void myGlutDisplay(void)
{
  //	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
  ::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
  ::glClearStencil(0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1f, 4.0f );
  window.SetGL_Camera();
  
  DrawBackground();
  
  ::glColor3d(0,0,0);
  DrawMeshTri2D_Edge(aTri,aXY0);
  {
    double scale = (aXY0.size()/2)*1.0e-4;
    assert( aMode.size()/3 == aXY0.size()/2 );
    ::glColor3d(1,0,0);
    ::glBegin(GL_LINES);
    for(int itri=0;itri<aTri.size()/3;itri++){
      const unsigned int i0 = aTri[itri*3+0];
      const unsigned int i1 = aTri[itri*3+1];
      const unsigned int i2 = aTri[itri*3+2];
      const double p0[3] = { aXY0[i0*2+0], aXY0[i0*2+1], aMode[i0*3+0]*scale };
      const double p1[3] = { aXY0[i1*2+0], aXY0[i1*2+1], aMode[i1*3+0]*scale };
      const double p2[3] = { aXY0[i2*2+0], aXY0[i2*2+1], aMode[i2*3+0]*scale };
      ::glVertex3dv( p0 );
      ::glVertex3dv( p1 );
      ::glVertex3dv( p1 );
      ::glVertex3dv( p2 );
      ::glVertex3dv( p2 );
      ::glVertex3dv( p0 );
    }
    for(int ip=0;ip<aXY0.size()/2;++ip){
      const double p0[3] = { aXY0[ip*2+0], aXY0[ip*2+1], aMode[ip*3+0]*scale };
      double rx = aMode[ip*3+1]*scale;
      double ry = aMode[ip*3+2]*scale;
      const double v[3] = {
        thickness*0.5*(+ry),
        thickness*0.5*(-rx),
        thickness*0.5};
      ::glVertex3d( p0[0]+v[0], p0[1]+v[1], p0[2]+v[2] );
      ::glVertex3d( p0[0]-v[0], p0[1]-v[1], p0[2]-v[2] );
    }
    ::glEnd();
  }

  ShowFPS();
  ::glutSwapBuffers();
}

void myGlutIdle(){
  
  ::glutPostRedisplay();
}


void myGlutResize(int w, int h)
{
  ::glViewport(0,0,w,h);
  ::glutPostRedisplay();
}

void myGlutSpecial(int Key, int x, int y)
{
  window.glutSpecial(Key, x, y);
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y )
{
  window.glutMotion(x, y);
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  window.glutMouse(button, state, x, y);
  ::glutPostRedisplay();
}

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
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
    case 'r':
    {
      for(int ip=0;ip<aXY0.size()/2;++ip){
        aTmp0[ip*3+0] = (rand()+1.0)/(RAND_MAX+1.0);
        aTmp0[ip*3+1] = (rand()+1.0)/(RAND_MAX+1.0);
        aTmp0[ip*3+2] = (rand()+1.0)/(RAND_MAX+1.0);
      }
      break;
    }
    case 'c':
    {
      for(int itr=0;itr<200;++itr){
        double C[3][2];
        for(int i=0;i<6;++i){
          (&C[0][0])[i] = 10.0*(rand()/(RAND_MAX+1.0)-0.5);
        }
        double a0 = TriArea2D(C[0], C[1], C[2]);
        if( a0 < 0.1 ) continue;
        double u[3][3];
        for(int i=0;i<9;++i){
          (&u[0][0])[i] = 1.0*(rand()/(RAND_MAX+1.0)-0.5);
        }
        double thickness = (rand()+1.0)/(RAND_MAX+1.0);
        double lambda = (rand()+1.0)/(RAND_MAX+1.0);
        double myu = (rand()+1.0)/(RAND_MAX+1.0);
        double diff = Check_WdWddW_PlateBendingMITC3(C, u,
                                                     thickness,lambda,myu, 1.0e-5);
        std::cout << itr << " " << diff << std::endl;
      }
    }
  }
  ::glutPostRedisplay();
}


int main(int argc,char* argv[])
{
  glutInit(&argc, argv);
  
  // Initialize GLUT window 3D
  glutInitWindowPosition(200,200);
  glutInitWindowSize(400, 300);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_STENCIL);
  glutCreateWindow("3D View");
  glutDisplayFunc(myGlutDisplay);
  glutIdleFunc(myGlutIdle);
  glutReshapeFunc(myGlutResize);
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);
  glutKeyboardFunc(myGlutKeyboard);
  glutSpecialFunc(myGlutSpecial);
  
  ////////////////////////
  
  window.camera.view_height = 1.0;
  window.camera.camera_rot_mode = CAMERA_ROT_ZTOP;
  setSomeLighting();
  
  MakeMesh();
  InitializeProblem_ShellEigenPB();
  for(int ip=0;ip<aXY0.size()/2;++ip){
    aTmp0[ip*3+0] = (rand()+1.0)/(RAND_MAX+1.0);
    aTmp0[ip*3+1] = (rand()+1.0)/(RAND_MAX+1.0);
    aTmp0[ip*3+2] = (rand()+1.0)/(RAND_MAX+1.0);
  }
  
  {
    const double I = thickness*thickness*thickness*leny/12.0;
    const double EI = EYoung*I;
    const double rhoA = rho*thickness*leny;
    // https://www.mapleprimes.com/DocumentFiles/206657_question/Transverse_vibration_of_beams.pdf
    // https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory
    freq_theo = (22.3733)/(lenx*lenx)*sqrt(EI/rhoA)/(2*M_PI);
  }
  
  glutMainLoop();
  return 0;
}


