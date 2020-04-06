/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <set>
#include "delfem2/mshtopo.h"
#include "delfem2/mats.h"
#include "delfem2/emat.h"

#include "delfem2/ilu_mats.h"
#include "delfem2/cloth_internal.h"
#include "delfem2/dtri_v2.h"

// ------------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// -------------------------------------------


class CInput_ContactNothing: public dfm2::CInput_Contact
{
public:
  double penetrationNormal(double& nx, double &ny, double& nz,
                           double px, double py, double pz) const
  {
    return -100;
  }
};


// ---------------------------

std::vector<dfm2::CDynPntSur> aPo2D;
std::vector<dfm2::CDynTri> aETri;
std::vector<dfm2::CVec2d> aVec2;
std::vector<double> aXYZ0; // undeformed vertex positions
std::vector<double> aXYZ; // deformed vertex positions
std::vector<double> aUVW; // deformed vertex velocity
std::vector<int> aBCFlag;  // boundary condition flag (0:free 1:fixed)
std::vector<unsigned int> aTri;  // index of triangles
std::vector<unsigned int> aQuad; // index of 4 vertices required for bending
dfm2::CMatrixSparse<double> mat_A; // coefficient matrix
dfm2::CPreconditionerILU<double>  ilu_A; // ilu decomposition of the coefficient matrix
double mass_point = 0.01;

int idp_nearest = -1;
int press_button = -1;
double mov_begin_x, mov_begin_y;
bool is_animation = true;
double mag = 1.0;

// ----------------------------------

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
    dfm2::CInputTriangulation_Uniform param(1.0);
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

// -----------------------------------------------

void myGlutDisplay()
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
  {
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0.8, 0.8, 0.8);
    /*
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
    ::glEnable(GL_DEPTH_TEST);
     */
    delfem2::opengl::DrawMeshTri3D_FaceNorm(aXYZ, aTri);
  }
  
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  delfem2::opengl::DrawMeshTri3D_Edge(aXYZ, aTri);

}

int main(int argc,char* argv[])
{
  double lenx = 1.0;
  {
    std::vector< std::vector<double> > aaXY;
    aaXY.resize(1);
    double xys[8] = {-0.5,-0.5, +0.5,-0.5, +0.5,+0.5, -0.5,+0.5};
    aaXY[0].assign(xys,xys+8);
    GenMesh(aaXY);
  }
  // -------------------------------
  aXYZ0.resize(aPo2D.size()*3);
  for(std::size_t ip=0;ip<aPo2D.size();++ip){
    aXYZ0[ip*3+0] = aVec2[ip].x();
    aXYZ0[ip*3+1] = aVec2[ip].y();
    aXYZ0[ip*3+2] = 0.0;
  }
  aXYZ = aXYZ0;
  /////
  aTri.resize(aETri.size()*3);
  for(std::size_t it=0;it<aETri.size();++it){
    aTri[it*3+0] = aETri[it].v[0];
    aTri[it*3+1] = aETri[it].v[1];
    aTri[it*3+2] = aETri[it].v[2];
  }
  dfm2::ElemQuad_DihedralTri(aQuad,
                             aTri.data(),aTri.size()/3,aXYZ0.size()/3);
  {
    const std::size_t np = aXYZ0.size()/3;
    mat_A.Initialize(np,3,true);
    std::vector<unsigned int> psup_ind,psup;
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                                                      aQuad.data(),aQuad.size()/4, 4, np);
    dfm2::JArray_Sort(psup_ind, psup);
    mat_A.SetPattern(psup_ind.data(),psup_ind.size(), psup.data(),psup.size());
    ilu_A.Initialize_ILU0(mat_A);
  }
  aUVW.resize(aXYZ.size(),0.0);
  aBCFlag.resize(aXYZ.size(),0);
  for(std::size_t ip=0;ip<aXYZ.size()/3;++ip){
    if( aXYZ[ip*3+0]  < -0.49*lenx ){
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aBCFlag[ip*3+2] = 1;
    }
  }

  // --------------------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  delfem2::opengl::setSomeLighting();
  while(!glfwWindowShouldClose(viewer.window)) {
    StepTime();
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
