/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#include "delfem2/mat3.h"
#include "delfem2/quat.h"
#include "delfem2/mats.h"
#include "delfem2/emat.h"
#include "delfem2/mshtopo.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/primitive.h"
#include "delfem2/vecxitrsol.h"
//
#include "delfem2/fem_emats.h"
#include "delfem2/ilu_mats.h"

// ----------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_color.h"
#include "delfem2/opengl/glold_v23.h"

namespace dfm2 = delfem2;

// ------------------------------------

void SetFixedBoundaryCondition(
    std::vector<double>& aRhs,
    unsigned int iframe,
    const std::vector<double>& aXYZ0,
    const std::vector<int>& aBCFlag)
{
  double A[16];
  {
    dfm2::AffMat3_Identity(A);
    const double trans0[3] = {0, -0.8, 0};
    dfm2::Translate_AffMat3(A,
                            trans0);
    const double axis0[3] = {0, +2.0*sin(0.03*iframe), 1.0*sin(0.07*iframe)};
    dfm2::Rotate_AffMat3_Rodriguez(A,
                                   axis0);
    const double trans1[3] = {0.2*sin(0.03*iframe), +0.6+0.2*cos(0.05*iframe), 0};
    dfm2::Translate_AffMat3(A,
                            trans1);
  }
  const unsigned int np = aRhs.size()/3;
  for(unsigned int ip=0;ip<np;++ip){
    if( aBCFlag[ip*3+0] == 0 ){ continue; }
    if( aBCFlag[ip*3+0] == 1 ){
      aRhs[ip*3+0] = aXYZ0[ip*3+0];
      aRhs[ip*3+1] = aXYZ0[ip*3+1];
      aRhs[ip*3+2] = aXYZ0[ip*3+2];
    }
    if( aBCFlag[ip*3+0] == 2 ) {
      dfm2::Vec3_AffMat3Vec3Projection(aRhs.data()+ip*3, A, aXYZ0.data()+ip*3);
    }
  }
}


// ------------------------------------

std::vector<unsigned int> aTri;
std::vector<double> aXYZ0;
std::vector<double> aXYZ1;
std::vector<int> aBCFlag;

std::vector<double> aQuat;

unsigned int imode = 0;
std::vector<unsigned int> psup_ind, psup;

// ------------------------------------------

void SetProblem()
{
  dfm2::MeshTri3D_CylinderClosed(aXYZ0, aTri,
                                 0.2, 1.6,
                                 16, 16);
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                             aTri.data(), aTri.size()/3, 3,
                             (int)aXYZ0.size()/3);
  dfm2::JArray_Sort(psup_ind, psup);

  {
    const unsigned int np = aXYZ0.size() / 3;
    aBCFlag.assign(np * 3, 0);
    for(unsigned int ip=0;ip<np;++ip) {
      double y0 = aXYZ0[ip*3+1];
      if( y0 < -0.65 ){
        aBCFlag[ip*3+0] = 1;
        aBCFlag[ip*3+1] = 1;
        aBCFlag[ip*3+2] = 1;
      }
      if( y0 > +0.65 ){
        aBCFlag[ip*3+0] = 2;
        aBCFlag[ip*3+1] = 2;
        aBCFlag[ip*3+2] = 2;
      }
    }
  }
  aXYZ1 = aXYZ0; // intial guess  
}

// -------------------------------------------

void myGlutDisplay()
{
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1,aTri);
  ::glColor3d(0.8,0.8,0.8);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ0.data(),aXYZ0.size()/3,
                                   aTri.data(),aTri.size()/3);
  ::glColor3d(0,0,0);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ1.data(),aXYZ1.size()/3,
                                   aTri.data(),aTri.size()/3);
  
  { // draw bc as a point
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    for(unsigned int ip=0;ip<aXYZ0.size()/3;++ip){
      if( aBCFlag[ip*3+0] == 0 ){ continue; }
      if( aBCFlag[ip*3+0] == 1 ){ ::glColor3d(0,0,1); }
      if( aBCFlag[ip*3+0] == 2 ){ ::glColor3d(0,1,0); }
      ::glVertex3dv(aXYZ1.data()+ip*3);
    }
    ::glEnd();
  }
  if( imode == 2 ){ // draw bc as a point
    ::glLineWidth(2);
    ::glBegin(GL_LINES);
    double l = 0.04;
    for(unsigned int ip=0;ip<aXYZ1.size()/3;++ip){
      const double* p = aXYZ1.data()+ip*3;
      {
        double ex0[3] = {1,0,0};
        double ex[3]; dfm2::QuatVec(ex, aQuat.data()+ip*4, ex0);
        ::glColor3d(1,0,0);
        ::glVertex3dv(p);
        ::glVertex3d(p[0]+l*ex[0], p[1]+l*ex[1], p[2]+l*ex[2]);
      }
      {
        double ey0[3] = {0,1,0};
        double ey[3]; dfm2::QuatVec(ey, aQuat.data()+ip*4, ey0);
        ::glColor3d(0,1,0);
        ::glVertex3dv(p);
        ::glVertex3d(p[0]+l*ey[0], p[1]+l*ey[1], p[2]+l*ey[2]);
      }
      {
        double ez0[3] = {0,0,1};
        double ez[3]; dfm2::QuatVec(ez, aQuat.data()+ip*4, ez0);
        ::glColor3d(0,0,1);
        ::glVertex3dv(p);
        ::glVertex3d(p[0]+l*ez[0], p[1]+l*ez[1], p[2]+l*ez[2]);
      }
    }
    ::glEnd();
  }
}

void ARAP_LinearDisponly(int iframe)
{
  const double weight_bc = 100.0;
  class CAtA {
  public:
    CAtA(const std::vector<unsigned int>& psup_ind0,
         const std::vector<unsigned int>& psup0,
         double weight_bc0,
         const std::vector<int>& aBCFlag0) :
    psup_ind(psup_ind0),
    psup(psup0),
    weight_bc(weight_bc0),
    aBCFlag(aBCFlag0)
    {
      const unsigned int np = aBCFlag.size()/3;
      const unsigned int ne = psup.size();
      // -----
      aMatEdge.resize(ne*9*2);
      assert(psup_ind.size()==np+1);
      for(unsigned int ip=0;ip<np;++ip){
        for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
//          unsigned int jp = psup[ipsup];
          dfm2::Mat3_Identity(aMatEdge.data()+ipsup*18,   +1.0);
          dfm2::Mat3_Identity(aMatEdge.data()+ipsup*18+9, -1.0);
        }
      }
      // ---
      vec_tmp.resize(ne*3);
    }
    void JacobiTVecTmp(double*y ,
                       double alpha, double beta) const {
      const unsigned int np = aBCFlag.size()/3;
      for(int i=0;i<np*3;++i){ y[i] *= beta; }
      for(unsigned int ip=0;ip<np;++ip){
        for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
          unsigned int jp0 = psup[ipsup];
          dfm2::MatTVec3_ScaleAdd(y+ip*3,
                                  aMatEdge.data()+ipsup*18,
                                  vec_tmp.data()+ipsup*3,
                                  alpha, 1.0);
          dfm2::MatTVec3_ScaleAdd(y+jp0*3,
                                  aMatEdge.data()+ipsup*18+9,
                                  vec_tmp.data()+ipsup*3,
                                  alpha, 1.0);
        }
      }
    }
    void MatVec(double* y,
                double alpha, const double* vec,  double beta) const {
      const unsigned int np = aBCFlag.size()/3;
      std::fill(vec_tmp.begin(),vec_tmp.end(), 0.0);
      for(unsigned int ip=0;ip<np;++ip){
        for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
          unsigned int jp0 = psup[ipsup];
          dfm2::MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
                                 aMatEdge.data()+ipsup*18,
                                 vec+ip*3,
                                 1.0, 1.0);
          dfm2::MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
                                 aMatEdge.data()+ipsup*18+9,
                                 vec+jp0*3,
                                 1.0, 1.0);
        }
      }
      this->JacobiTVecTmp(y,
                          alpha, beta);
      // add diagonal for fixed boundary condition
      for(int i=0;i<aBCFlag.size();++i){
        if( aBCFlag[i] == 0 ){ continue; }
        y[i] += weight_bc*vec[i];
      }
    }
    void MakeRHS(double* aRhs,
                 const double* aXYZ0,
                 const double* aXYZ1,
                 const double* aGoal) const
    {
      const unsigned int np = aBCFlag.size()/3;
      const unsigned int ne = psup.size();
      vec_tmp.assign(ne*3,0);
      for(unsigned int ip=0;ip<np;++ip){
        for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
          const unsigned int jp0 = psup[ipsup];
          const double d0[3] = { aXYZ0[jp0*3+0]-aXYZ0[ip*3+0], aXYZ0[jp0*3+1]-aXYZ0[ip*3+1], aXYZ0[jp0*3+2]-aXYZ0[ip*3+2] };
          const double d1[3] = { aXYZ1[jp0*3+0]-aXYZ1[ip*3+0], aXYZ1[jp0*3+1]-aXYZ1[ip*3+1], aXYZ1[jp0*3+2]-aXYZ1[ip*3+2] };
          vec_tmp[ipsup*3+0] += +(d0[0] - d1[0]);
          vec_tmp[ipsup*3+1] += +(d0[1] - d1[1]);
          vec_tmp[ipsup*3+2] += +(d0[2] - d1[2]);
        }
      }
      this->JacobiTVecTmp(aRhs,
                          -1.0, 0.0);
     // making RHS vector for fixed boundary condition
      for(int i=0;i<np*3;++i){
        if( aBCFlag[i] == 0 ){ continue; }
        aRhs[i] += (aGoal[i]-aXYZ1[i])*weight_bc;
      }

    }
  public:
    const std::vector<unsigned int>& psup_ind;
    const std::vector<unsigned int>& psup;
    const double weight_bc;
    const std::vector<int>& aBCFlag;
    // -------------
    std::vector<double> aMatEdge;
    mutable std::vector<double> vec_tmp;
  } MatAtA(psup_ind, psup, weight_bc, aBCFlag);
  // -----
  for(int i=0;i<aXYZ0.size();++i){ // adding noise for debuggng purpose
    aXYZ1[i] += 0.02*(double)rand()/(RAND_MAX+1.0)-0.01;
  }
  std::vector<double> aGoal(aXYZ0.size());
  SetFixedBoundaryCondition(aGoal,
                            iframe,aXYZ0,aBCFlag);
// making RHS vector for elastic deformation
  std::vector<double> aRhs(aXYZ0.size(),0.0);
  MatAtA.MakeRHS(aRhs.data(),
                 aXYZ0.data(), aXYZ1.data(), aGoal.data());
  std::vector<double> aUpd(aXYZ0.size(),0.0);
  std::vector<double> aRes = dfm2::Solve_CG(aRhs.data(), aUpd.data(),
                                            aRhs.size(), 1.0e-7, 300, MatAtA);
  std::cout << aRes.size() << std::endl;
  for(int i=0;i<aBCFlag.size();++i){ aXYZ1[i] += aUpd[i]; }
}

// --------------------------------------------------

int main(int argc,char* argv[])
{
  SetProblem();
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();
  int iframe = 0;
  while (!glfwWindowShouldClose(viewer.window))
  {
    ARAP_LinearDisponly(iframe);
    // ------
    iframe++;
    viewer.DrawBegin_oldGL();
    myGlutDisplay();
    viewer.DrawEnd_oldGL();
  }
  
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


