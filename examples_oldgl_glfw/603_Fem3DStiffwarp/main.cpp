/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/opengl/glfw/viewer_glfw.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/ilu_mats.h"
#include "delfem2/fem_emats.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/mshmisc.h"
#include "delfem2/mshuni.h"
#include "delfem2/mats.h"
#include "delfem2/dtri.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/jagarray.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>

namespace dfm2 = delfem2;

// ----------------------------------

void GenMesh
(std::vector<dfm2::CVec2d>& aVec2,
 std::vector<dfm2::CDynPntSur>& aPo2D,
 std::vector<dfm2::CDynTri>& aETri,
 double elen,
 const std::vector< std::vector<double> >& aaXY)
{
  std::vector<int> loopIP_ind, loopIP;
  {
    dfm2::JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                               aaXY);
    if( !dfm2::CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
      return;
    }
    dfm2::FixLoopOrientation(loopIP,
                             loopIP_ind,aVec2);
    if( elen > 10e-10 ){
      dfm2::ResamplingLoop(loopIP_ind,loopIP,aVec2,
                           elen );
    }
  }
  ////
  Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                 loopIP_ind,loopIP);
  if( elen > 1.0e-10 ){
    dfm2::CInputTriangulation_Uniform param(1.0);
    std::vector<int> aFlgPnt(aVec2.size());
    std::vector<unsigned int> aFlgTri(aETri.size(),0);
    MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
                  aVec2.size(), 0, elen, param);
  }
}

void RotationAtMeshPoints
(std::vector<double>& aR,
 const std::vector<double>& aXYZ,
 const std::vector<double>& aDisp,
 const std::vector<unsigned int> &psup_ind,
 const std::vector<unsigned int> &psup)
{
  const unsigned int np = aXYZ.size()/3;
  aR.resize(np*9);
  for(std::size_t ip=0;ip<aXYZ.size()/3;++ip){
    dfm2::CVec3d Pi(aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]);
    dfm2::CVec3d pi(aXYZ[ip*3+0]+aDisp[ip*3+0],
                aXYZ[ip*3+1]+aDisp[ip*3+1],
                aXYZ[ip*3+2]+aDisp[ip*3+2]);
    dfm2::CMat3d A;
    A.SetZero();
    for(unsigned int jjp=psup_ind[ip];jjp<psup_ind[ip+1];++jjp){
      int jp = psup[jjp];
      dfm2::CVec3d Pj(aXYZ[jp*3+0],aXYZ[jp*3+1],aXYZ[jp*3+2]);
      dfm2::CVec3d pj(aXYZ[jp*3+0]+aDisp[jp*3+0],
                  aXYZ[jp*3+1]+aDisp[jp*3+1],
                  aXYZ[jp*3+2]+aDisp[jp*3+2]);
      A += dfm2::Mat3_OuterProduct(pj-pi,Pj-Pi);
    }
    dfm2::GetRotPolarDecomp(aR.data()+ip*9,
                            A.mat, 100);
  }
}


// ---------------------------

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

dfm2::CMatrixSparse<double> mat_A;
std::vector<double> vec_b;
dfm2::CPreconditionerILU<double>  ilu_A;
std::vector<unsigned int> psup_ind, psup;
std::vector<double> aR;

// --------------------------------------------


void InitializeProblem_ShellEigenPB()
{
  const unsigned int np = (int)aXYZ.size()/3;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind, psup,
      aTet.data(), aTet.size()/4, 4,
      (int)aXYZ.size()/3);
  dfm2::JArray_Sort(psup_ind, psup);
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  ilu_A.Initialize_ILU0(mat_A);
}
  
// ------------------------------------------------------

void Solve_Linear()
{
  mat_A.SetZero();
  vec_b.assign(aXYZ.size(),0.0);
  dfm2::MergeLinSys_SolidLinear_BEuler_MeshTet3D(mat_A, vec_b.data(),
                                                 myu, lambda,
                                                 rho, gravity,
                                                 dt,
                                                 aXYZ.data(), aXYZ.size()/3,
                                                 aTet.data(), aTet.size()/4,
                                                 aDisp.data(),
                                                 aVelo.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b,aBCFlag,0);
  //
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  const int nDoF = aXYZ.size();
  std::vector<double> dv(nDoF,0.0);
  std::vector<double> aConv = Solve_PBiCGStab(vec_b.data(), dv.data(),
                                              1.0e-4, 1000, mat_A, ilu_A);
  //
  dfm2::XPlusAYBZ(aDisp,nDoF,aBCFlag,
                  dt, dv,
                  dt, aVelo);
  dfm2::XPlusAY(aVelo,nDoF,aBCFlag,
                1.0, dv);
  std::cout << "conv; " << aConv.size() <<  std::endl;
}


void Solve_StiffnessWarping()
{
  RotationAtMeshPoints(aR,
                       aXYZ,aDisp,psup_ind,psup);
  // ----------------------
  mat_A.SetZero();
  vec_b.assign(aXYZ.size(),0.0);
  dfm2::MergeLinSys_SolidStiffwarp_BEuler_MeshTet3D(mat_A, vec_b.data(),
                                                    myu, lambda,
                                                    rho, gravity,
                                                    dt,
                                                    aXYZ.data(), aXYZ.size()/3,
                                                    aTet.data(), aTet.size()/4,
                                                    aDisp.data(),
                                                    aVelo.data(),
                                                    aR);
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b,aBCFlag,0);
  // -------------------
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  const int nDoF = aXYZ.size();
  std::vector<double> dv(nDoF,0.0);
  std::vector<double> aConv = Solve_PBiCGStab(vec_b.data(), dv.data(),
                                              1.0e-4, 1000, mat_A, ilu_A);
  dfm2::XPlusAYBZ(aDisp,nDoF,aBCFlag,
                  dt, dv,
                  dt, aVelo);
  dfm2::XPlusAY(aVelo,nDoF,aBCFlag,
                1.0,dv);
  std::cout << "conv; " << aConv.size() <<  std::endl;
}

// --------------------------------------------------------------

void myGlutDisplay()
{
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
    delfem2::opengl::DrawMeshTet3D_EdgeDisp(aXYZ.data(),
                                            aTet.data(),aTet.size()/4,
                                            aDisp.data(),
                                            1.0);
  }
  
  ::glDisable(GL_LIGHTING);
  for(std::size_t ip=0;ip<aXYZ.size()/3;++ip){
    dfm2::CVec3d pi(aXYZ[ip*3+0]+aDisp[ip*3+0],
                      aXYZ[ip*3+1]+aDisp[ip*3+1],
                      aXYZ[ip*3+2]+aDisp[ip*3+2]);
    dfm2::CVec3d ex(aR[ip*9+0],aR[ip*9+3],aR[ip*9+6]);
    dfm2::CVec3d ey(aR[ip*9+1],aR[ip*9+4],aR[ip*9+7]);
    dfm2::CVec3d ez(aR[ip*9+2],aR[ip*9+5],aR[ip*9+8]);
    ::glBegin(GL_LINES);
    ::glColor3d(1,0,0);
    delfem2::opengl::myGlVertex(pi);
    delfem2::opengl::myGlVertex(pi+0.04*ex);
    ::glColor3d(0,1,0);
    delfem2::opengl::myGlVertex(pi);
    delfem2::opengl::myGlVertex(pi+0.04*ey);
    ::glColor3d(0,0,1);
    delfem2::opengl::myGlVertex(pi);
    delfem2::opengl::myGlVertex(pi+0.04*ez);
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
    delfem2::opengl::DrawMeshTet3D_FaceNorm(aXYZ.data(),
                                   aTet.data(), aTet.size()/4);

  }

}

int main(int argc,char* argv[])
{
  {
    std::vector< std::vector<double> > aaXY;
    {
      const double aXY[8] = {
        -1,-0.1,
        +1,-0.1,
        +1,+0.1,
        -1,+0.1 };
      aaXY.emplace_back(aXY,aXY+8 );
    }
    std::vector<dfm2::CVec2d> aVec2;
    std::vector<dfm2::CDynPntSur> aPo2D;
    std::vector<dfm2::CDynTri> aETri;
    GenMesh(aVec2,aPo2D,aETri,
            0.075, aaXY);
    std::vector<double> aXY;
    std::vector<unsigned int> aTri;
    CMeshTri2D(aXY,aTri,
               aVec2,aETri);
    dfm2::ExtrudeTri2Tet(3, 0.075,
        aXYZ,aTet,
        aXY,aTri);
  }
  aDisp.assign(aXYZ.size(), 0.0);
  aVelo.assign(aXYZ.size(), 0.0);
  aBCFlag.assign(aXYZ.size(),0);
  for(std::size_t ip=0;ip<aXYZ.size()/3;++ip){
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

  delfem2::opengl::CViewer_GLFW viewer;
  viewer.camera.view_height = 2.0;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  viewer.Init_oldGL();
  delfem2::opengl::setSomeLighting();

  while(!glfwWindowShouldClose(viewer.window)){
    if( is_stiffness_warping ){ Solve_StiffnessWarping(); }
    else{                       Solve_Linear();           }
    // -----
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);

}
