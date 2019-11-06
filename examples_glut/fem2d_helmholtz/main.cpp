#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits>
#include <vector>
#include <complex>
#include "delfem2/mshtopo.h"
#include "delfem2/msh.h"
#include "delfem2/dtri.h"
#include "delfem2/mats.h"

#include "delfem2/dtri_v2.h"
#include "delfem2/ilu_mats.h"
#include "delfem2/fem_emats.h"
#include "delfem2/cad2d.h"

// ---------------

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem2/opengl/gl2_color.h"
#include "delfem2/opengl/gl2_funcs.h"
#include "../glut_cam.h"

namespace dfm2 = delfem2;

// --------------------------------------------

CNav3D_GLUT nav;
std::vector<unsigned int> aTri1;
std::vector<double> aXY1;
std::vector<int> aBCFlag; // boundary condition flag
std::vector< std::vector<unsigned int> > aaIP;
int ipCenter = -1;

std::vector<std::complex<double> > aCVal;
std::vector<double> aVal;
double time_cur  = 0.0;

CMatrixSparse<std::complex<double> > mat_A;
std::vector<std::complex<double> > vec_b;
CPreconditionerILU<std::complex<double> > ilu_A;

bool is_animation = false;


//////////////////////////////////////////////////////////////////////////////////////

void MakeMesh(){
  CCad2D cad2d;
  {
    double xy[8] = {-1,-1, +1,-1, +1,+1, -1,+1};
    std::vector<double> aXY(xy,xy+8);
    cad2d.AddPolygon(aXY);
  }
  cad2d.AddVtxFace(0.0, 0.0, 0);
  CMesher_Cad2D mshr;
  mshr.edge_length = 0.05;
  CMeshDynTri2D dmsh;
  mshr.Meshing(dmsh, cad2d);
  MeshTri2D_Export(aXY1, aTri1,
                   dmsh.aVec2, dmsh.aETri);
  
  aaIP.resize(4);
  aaIP[0] = mshr.IndPoint_IndEdge(0, true, cad2d);
  aaIP[1] = mshr.IndPoint_IndEdge(1, true, cad2d);
  aaIP[2] = mshr.IndPoint_IndEdge(2, true, cad2d);
  aaIP[3] = mshr.IndPoint_IndEdge(3, true, cad2d);
  ipCenter = 4;
}

//////////////////////////////////////////////////////////////////////////////////////////
// iproblem: 0, 1
void InitializeProblem_Scalar()
{
  const int np = (int)aXY1.size()/2;
  aCVal.assign(np, std::complex<double>(0.0));
  aBCFlag.resize(np,0);
  aBCFlag[ipCenter] = 1;
  aCVal[ipCenter] = 1;
  //////
  std::vector<int> psup_ind, psup;
  JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                              aTri1.data(), aTri1.size()/3, 3, (int)aXY1.size()/2);
  JArray_Sort(psup_ind, psup);
  ////
  mat_A.Initialize(np, 1, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(),psup.size());
  ilu_A.Initialize_ILU0(mat_A);
//  ilu_A.Initialize_ILUk(mat_A, 2);
}

// --------------------------------------
// iproblem: 0
void SolveProblem_Poisson()
{
  const int np = (int)aXY1.size()/2;
  const int nDoF = np;
  // ------------
  const double wave_length = 0.4;
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_Helmholtz_MeshTri2D(mat_A,vec_b.data(),
                                        wave_length,
                                        aXY1.data(),aXY1.size()/2,
                                        aTri1.data(),aTri1.size()/3,
                                        aCVal.data());

  for(int ipl=0;ipl<aaIP.size();++ipl){
    dfm2::MergeLinSys_SommerfeltRadiationBC_Polyline2D(mat_A,vec_b.data(),
                                                       wave_length,
                                                       aXY1.data(),aXY1.size()/2,
                                                       aaIP[ipl].data(),aaIP[ipl].size(),
                                                       aCVal.data());
  }


  mat_A.SetBoundaryCondition(aBCFlag.data(),aBCFlag.size(),1);
  setRHS_Zero(vec_b, aBCFlag,0);
  ///////////////////////////
  std::vector<std::complex<double> > vec_x;
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  vec_x.resize(vec_b.size());
  /*
  std::vector<double> aConv = Solve_PBiCGStab(vec_b.data(), vec_x.data(),
                                              1.0e-4, 400, mat_A, ilu_A);
   */
  /*
   std::vector<double> aConv = Solve_BiCGSTAB(vec_b, vec_x,
   1.0e-4,400, mat_A);
   */
  std::vector<double> aConv = Solve_PCOCG(vec_b.data(), vec_x.data(),
                                          1.0e-4, 400, mat_A, ilu_A);
  std::cout << aConv.size() << " " << aConv[ aConv.size()-1 ] << std::endl;

  for(int ic=0;ic<aConv.size();++ic){
    std::cout << ic << " " << aConv[ic] << std::endl;
  }

//  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  ///////////////////////////
  XPlusAY(aCVal,
          nDoF,aBCFlag, std::complex<double>(1.0),vec_x);
}


//////////////////////////////////////////////////////////////////////////////////////

void myGlutResize(int w, int h)
{
  glViewport(0, 0, w, h);
  glutPostRedisplay();
}

void myGlutDisplay(void)
{
  //  ::glClearColor(0.2, .7, 0.7, 1.0);
  ::glClearColor(1.0, 1.0, 1.0, 1.0);
  ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  ::glEnable(GL_DEPTH_TEST);
  
  ::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 1.1, 4.0 );
  
  nav.SetGL_Camera();
  
  delfem2::opengl::DrawMeshTri2D_Edge(aTri1,aXY1);
  ::glPointSize(2);
  ::glColor3d(0,0,0);
  delfem2::opengl::DrawPoints2D_Points(aXY1);
  
  std::vector< std::pair<double,delfem2::CColor> > colorMap;
//  makeHeatMap_BlueGrayRed(colorMap, -0.2, +0.2);
  delfem2::ColorMap_BlueCyanGreenYellowRed(colorMap, -0.2, +0.2);
  delfem2::opengl::DrawMeshTri2D_ScalarP1(aXY1.data(),aXY1.size()/2,
                                 aTri1.data(),aTri1.size()/3,
                                 aVal.data(),1,colorMap);
  
  ::glColor3d(0,0,0);
  ShowFPS();
  
  glutSwapBuffers();
}

void myGlutIdle(){
  if( is_animation ){
    time_cur += 0.3;
    std::complex<double> rot(cos(time_cur),sin(time_cur));
    for(int ip=0;ip<aVal.size();++ip){
      aVal[ip] = (rot*aCVal[ip]).real();
    }
  }
  ::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
  nav.glutMotion(x,y);
}

void myGlutMouse(int button, int state, int x, int y){
    nav.glutMouse(button,state,x,y);
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
      is_animation = !is_animation; 
      break;
    default:
      break;
  }
}

void myGlutSpecial(int key, int x, int y){
  nav.glutSpecial(key,x,y);
  switch(key){
  }
  ::myGlutResize(-1,-1);
  ::glutPostRedisplay();
}

int main(int argc,char* argv[])
{
  // Initailze GLUT
  ::glutInitWindowPosition(200,200);
  ::glutInitWindowSize(1000, 500);
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
  
  MakeMesh();
  InitializeProblem_Scalar();
  SolveProblem_Poisson();
  
  aVal.resize(aCVal.size(),0.1);
  time_cur = 0.0;
  for(int ip=0;ip<aVal.size();++ip){ aVal[ip] = aCVal[ip].real(); }
  
  // Enter main loop
  nav.camera.view_height = 2.0;
  ::glutMainLoop();
  return 0;
}
