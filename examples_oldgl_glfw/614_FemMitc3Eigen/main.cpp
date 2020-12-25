/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/lsilu_mats.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/lsmats.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/fem_emats.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"
#include "delfem2/vec2.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/jagarray.h"
#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>
#include <random>


#ifndef M_PI
#  define M_PI 3.14159265359
#endif

namespace dfm2 = delfem2;

// ----------------------------------


void SetValue_ShellPBMITC3Eigen_MassLumpedSqrtInv_KernelModes3
(double* aMassLumpedSqrtInv,
 double* aModesKer,
 double rho, double thickness,
 const double* aXY, int nXY,
 const unsigned int* aTri, int nTri)
{
  const unsigned int nDoF = nXY*3;
  std::vector<double> aMassLumpedSqrt(nDoF);
  dfm2::MassLumped_ShellPlateBendingMITC3(aMassLumpedSqrt.data(),
                                          rho, thickness,
                                          aXY, nXY,
                                          aTri, nTri);
  for(unsigned int i=0;i<nDoF;++i){ aMassLumpedSqrt[i] = sqrt(aMassLumpedSqrt[i]); }
  {
    for(unsigned int i=0;i<nDoF*3;++i){ aModesKer[i] = 0.0; }
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
    dfm2::NormalizeX(p0,nDoF);
    dfm2::OrthogonalizeToUnitVectorX(p1, p0, nDoF);
    dfm2::OrthogonalizeToUnitVectorX(p2, p0, nDoF);
    dfm2::NormalizeX(p1,nDoF);
    dfm2::OrthogonalizeToUnitVectorX(p2, p1, nDoF);
    dfm2::NormalizeX(p2,nDoF);
  }
  for(unsigned int i=0;i<nDoF;++i){ aMassLumpedSqrtInv[i] = 1.0/aMassLumpedSqrt[i]; }
}

void RemoveKernel(std::vector<double>& aTmp0,
                  const std::vector<double>& aModesKer)
{
  const unsigned int nDoF = aTmp0.size();
  const double* p0 = aModesKer.data()+nDoF*0;
  const double* p1 = aModesKer.data()+nDoF*1;
  const double* p2 = aModesKer.data()+nDoF*2;
  double* p = aTmp0.data();
  dfm2::OrthogonalizeToUnitVectorX(p, p0, nDoF);
  dfm2::OrthogonalizeToUnitVectorX(p, p1, nDoF);
  dfm2::OrthogonalizeToUnitVectorX(p, p2, nDoF);
  dfm2::NormalizeX(p, nDoF);
}



// ----------------------

// display data
bool is_animation;

std::vector<unsigned int> aTri;
std::vector<double> aXY0;
std::vector<double> aMassLumpedSqrtInv;
std::vector<double> aTmp0;
std::vector<double> aTmp1;
std::vector<double> aMode;
std::vector<double> aModesKer;

dfm2::CMatrixSparse<double> mat_A;
dfm2::CPreconditionerILU<double> ilu_A;

const double lenx = 0.3;
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

// --------------------------------------

void MakeMesh(){
  std::vector< std::vector<double> > aaXY;
  {
    aaXY.resize(1);
    aaXY[0].push_back(-lenx*0.5); aaXY[0].push_back(-leny*0.5);
    aaXY[0].push_back(+lenx*0.5); aaXY[0].push_back(-leny*0.5);
    aaXY[0].push_back(+lenx*0.5); aaXY[0].push_back(+leny*0.5);
    aaXY[0].push_back(-lenx*0.5); aaXY[0].push_back(+leny*0.5);
  }
  // ----------------------------
  std::vector<dfm2::CDynPntSur> aPo2D;
  std::vector<dfm2::CDynTri> aETri;
  std::vector<dfm2::CVec2d> aVec2;
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
  // --------------------------------------
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                             aTri.data(), aTri.size()/3, 3,
                             (int)aXY0.size()/2);
  dfm2::JArray_Sort(psup_ind, psup);
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  ilu_A.Initialize_ILU0(mat_A);
  // ------------------------------
  aMassLumpedSqrtInv.resize(nDoF);
  aModesKer.resize(nDoF*3);
  SetValue_ShellPBMITC3Eigen_MassLumpedSqrtInv_KernelModes3(
      aMassLumpedSqrtInv.data(),
      aModesKer.data(),
      rho, thickness,
      aXY0.data(), aXY0.size()/2,
      aTri.data(), aTri.size()/3);
  // -------------------------------
  mat_A.SetZero();
  aMode.assign(nDoF, 0.0);
  aTmp0.assign(nDoF, 0.0);
  dfm2::MergeLinSys_ShellStaticPlateBendingMITC3_MeshTri2D(
      mat_A, aMode.data(),
      thickness,lambda, myu, 0.0, 0.0,
      aXY0.data(), aXY0.size()/2,
      aTri.data(), aTri.size()/3,
      aTmp0.data());
  MatSparse_ScaleBlkLen_LeftRight(
      mat_A,
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
  {
    const std::size_t n = aTmp1.size();
    std::vector<double> tmp0(n), tmp1(n);
    auto vr = dfm2::CVecXd(aTmp1);
    auto vu = dfm2::CVecXd(aMode);
    auto vs = dfm2::CVecXd(tmp0);
    auto vt = dfm2::CVecXd(tmp1);
    aConv = dfm2::Solve_PCG(
        vr,vu,vs,vt,
        conv_ratio, iteration, mat_A, ilu_A);
  }
  {
    double lam0 = dfm2::DotX(aTmp0.data(), aMode.data(), aTmp0.size());
    double freq_sim = sqrt(1.0/lam0-offset_dia)/(2*M_PI);
    std::cout << "freq theo" << freq_theo << "   freq_sim:" << freq_sim << "   " << freq_theo/freq_sim  << "     " << aConv.size() << std::endl;
  }
  aTmp0 = aMode;
  //
  RemoveKernel(aTmp0,
               aModesKer);
  //
  for(unsigned int i=0;i<aTmp0.size();++i){
    aMode[i] = aTmp0[i]*aMassLumpedSqrtInv[i];
  }
}

// ---------------------------------------------

void myGlutDisplay(void)
{
  delfem2::opengl::DrawBackground();
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  delfem2::opengl::DrawMeshTri2D_Edge(aTri,aXY0);
  {
    double scale = (aXY0.size()/2)*1.0e-4;
    assert( aMode.size()/3 == aXY0.size()/2 );
    ::glColor3d(1,0,0);
    ::glBegin(GL_LINES);
    for(unsigned int itri=0;itri<aTri.size()/3;itri++){
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
    for(unsigned int ip=0;ip<aXY0.size()/2;++ip){
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
      for(unsigned int ip=0;ip<aXY0.size()/2;++ip){
        aTmp0[ip*3+0] = (rand()+1.0)/(RAND_MAX+1.0);
        aTmp0[ip*3+1] = (rand()+1.0)/(RAND_MAX+1.0);
        aTmp0[ip*3+2] = (rand()+1.0)/(RAND_MAX+1.0);
      }
      break;
    }
  }
}


int main(int argc,char* argv[])
{
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  // --------------------------------
  viewer.camera.view_height = 0.2;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::ZTOP;
  viewer.camera.Rot_Camera(0.5, 0.5);
  delfem2::opengl::setSomeLighting();
  
  MakeMesh();
  InitializeProblem_ShellEigenPB();
  {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> score(0.0, 1.0);
    for (std::size_t ip = 0; ip < aXY0.size() / 2; ++ip) {
      aTmp0[ip * 3 + 0] = score(mt);
      aTmp0[ip * 3 + 1] = score(mt);
      aTmp0[ip * 3 + 2] = score(mt);
    }
  }
  
  {
    const double I = thickness*thickness*thickness*leny/12.0;
    const double EI = EYoung*I;
    const double rhoA = rho*thickness*leny;
    // https://www.mapleprimes.com/DocumentFiles/206657_question/Transverse_vibration_of_beams.pdf
    // https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory
    freq_theo = (22.3733)/(lenx*lenx)*sqrt(EI/rhoA)/(2*M_PI);
  }
  
  while (true)
  {
    Solve();
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ) goto EXIT;
  }
  
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


