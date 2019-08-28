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
#include "delfem2/ilu_mats.h"
#include "delfem2/vec2.h"

#include "delfem2/dyntri_v2.h"
#include "delfem2/fem_emats.h"

#include "delfem2/gl_funcs.h"
#include "delfem2/gl_color.h"
#include "delfem2/glut_funcs.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////

// display data
bool is_animation;

CGlutWindowManager window;


std::vector<unsigned int> aTri;
std::vector<double> aXY0;

std::vector<double> aVal;
std::vector<double> aVelo;
std::vector<double> aAcc;
std::vector<int> aBCFlag; // boundary condition flag
std::vector<int> aMSFlag; // master slave flag

CMatrixSparse<double> mat_A;
std::vector<double> vec_b;
CPreconditionerILU<double> ilu_A;

const double lenx = 1.0;
const double leny = 0.1;

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
  std::vector<CVector2> aVec2;
  std::vector<int> loopIP_ind, loopIP; // vtx on loop
  const double elen = 0.05;
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
  {
    std::vector<CEPo2> aPo2D;
    std::vector<ETri> aETri;
    Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                   loopIP_ind,loopIP);
    if( elen > 1.0e-10 ){
      CInputTriangulation_Uniform param(1.0);
      std::vector<int> aFlgPnt(aPo2D.size()), aFlgTri(aETri.size());
      MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
                    aVec2.size(),0,elen, param);
    }
    MeshTri2D_Export(aXY0,aTri,
                     aVec2,aETri);
  }
  std::cout<<"  ntri;"<<aTri.size()/3<<"  nXY:"<<aXY0.size()/2<<std::endl;
}

void InitializeProblem_PlateBendingMITC3()
{
  const int np = (int)aXY0.size()/2;
  aBCFlag.assign(np*3, 0);
  for(int ip=0;ip<np;++ip){
    const double px = aXY0[ip*2+0];
//    const double py = aXY0[ip*2+1];
    if( fabs(px-(-lenx*0.5)) < 0.0001 ){
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aBCFlag[ip*3+2] = 1;
    }
  }
  //////
  std::vector<int> psup_ind, psup;
  JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                              aTri.data(), aTri.size()/3, 3,
                                              (int)aXY0.size()/2);
  JArray_Sort(psup_ind, psup);
  ////
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(),psup.size());
//  ilu_A.Initialize_ILU0(mat_A);
    ilu_A.Initialize_ILUk(mat_A,0);
}

void SolveProblem_PlateBendingMITC3()
{
  const int np = (int)aXY0.size()/2;
  const int nDoF = np*3;
  //////////////////////////
  double thickness = 0.1;
  double myu = 10000.0;
  double lambda = 0.0;
  double rho = 1.0;
  double gravity_z = -10.0;
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  MergeLinSys_ShellStaticPlateBendingMITC3_MeshTri2D(mat_A,vec_b.data(),
                                                     thickness,lambda,myu,
                                                     rho,gravity_z,
                                                     aXY0.data(), aXY0.size()/2,
                                                     aTri.data(), aTri.size()/3,
                                                     aVal.data());
  std::cout << Dot(vec_b, vec_b) << std::endl;
  mat_A.SetBoundaryCondition(aBCFlag.data(),aBCFlag.size()/3,3);
  setRHS_Zero(vec_b, aBCFlag,0);
  //////////////////////////
  std::vector<double> vec_x;
  {
    ilu_A.SetValueILU(mat_A);
    ilu_A.DoILUDecomp();
    vec_x.resize(vec_b.size());
    std::vector<double> conv = Solve_PCG(vec_b.data(), vec_x.data(), 1.0e-5, 100,
                                         mat_A, ilu_A);
    std::cout << "convergence   nitr:" << conv.size() << "    res:" << conv[conv.size()-1] << std::endl;
  }
  //////////////////////////
  XPlusAY(aVal,nDoF,aBCFlag,
          1.0,vec_x);
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
    assert( aVal.size()/3 == aXY0.size()/2 );
    ::glColor3d(1,0,0);
    ::glBegin(GL_LINES);
    for(int itri=0;itri<aTri.size()/3;itri++){
      const unsigned int i0 = aTri[itri*3+0];
      const unsigned int i1 = aTri[itri*3+1];
      const unsigned int i2 = aTri[itri*3+2];
      const double p0[3] = { aXY0[i0*2+0], aXY0[i0*2+1], aVal[i0*3+0] };
      const double p1[3] = { aXY0[i1*2+0], aXY0[i1*2+1], aVal[i1*3+0] };
      const double p2[3] = { aXY0[i2*2+0], aXY0[i2*2+1], aVal[i2*3+0] };
      ::glVertex3dv( p0 );
      ::glVertex3dv( p1 );
      ::glVertex3dv( p1 );
      ::glVertex3dv( p2 );
      ::glVertex3dv( p2 );
      ::glVertex3dv( p0 );
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
  aVal.assign(aXY0.size()/2*3, 0.0);
  InitializeProblem_PlateBendingMITC3();
  SolveProblem_PlateBendingMITC3();
  
  glutMainLoop();
  return 0;
}


