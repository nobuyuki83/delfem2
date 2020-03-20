/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#include "delfem2/mat3.h"
#include "delfem2/mat4.h"
#include "delfem2/mats.h"
#include "delfem2/mshtopo.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"
#include "delfem2/primitive.h"
#include "delfem2/vecxitrsol.h"

// ----------------
#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_color.h"
#include "delfem2/opengl/glold_v23.h"

namespace dfm2 = delfem2;

// ------------------------------------

void SetLinSys_LaplaceGraph_MeshTri3(
    dfm2::CMatrixSparse<double>& mat_A)
{
  mat_A.SetZero();
  for(unsigned int ip=0;ip<mat_A.nblk_col;++ip){
    unsigned int n1r = mat_A.colInd[ip+1] - mat_A.colInd[ip];
    for(int icrs=mat_A.colInd[ip];icrs<mat_A.colInd[ip+1];++icrs){
      mat_A.valCrs[icrs*9+0*3+0] = -1;
      mat_A.valCrs[icrs*9+1*3+1] = -1;
      mat_A.valCrs[icrs*9+2*3+2] = -1;
    }
    mat_A.valDia[ip*9+0*3+0] = n1r;
    mat_A.valDia[ip*9+1*3+1] = n1r;
    mat_A.valDia[ip*9+2*3+2] = n1r;
  }
}

void SetDisplacementAtFixedBoundary(
    std::vector<double>& aRhs,
    unsigned int iframe,
    const std::vector<double>& aXYZ0,
    const std::vector<int>& aBCFlag)
{
  double A[16];
  {
    dfm2::Mat4_Identity(A);
    const double trans0[3] = {0, -0.8, 0};
    dfm2::Translate_Mat4Affine(A,
                            trans0);
    const double axis0[3] = {0, +2.0*sin(0.03*iframe), 1.0*sin(0.07*iframe)};
    dfm2::Rotate_Mat4AffineRodriguez(A,
                                   axis0);
    const double trans1[3] = {0.2*sin(0.03*iframe), +0.6+0.2*cos(0.05*iframe), 0};
    dfm2::Translate_Mat4Affine(A,
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
      dfm2::Vec3_Mat4Vec3_AffineProjection(aRhs.data()+ip*3, A, aXYZ0.data()+ip*3);
    }
  }
}

void SetProblem(
    std::vector<unsigned int>& aTri,
    std::vector<double>& aXYZ0,
    std::vector<int>& aBCFlag)
{
  dfm2::MeshTri3D_CylinderClosed(aXYZ0, aTri,
                                 0.2, 1.6,
                                 16, 16);
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
}


// ------------------------------------
// ------------------------------------------



// -------------------------------------------

void myGlutDisplay(
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aXYZ1,
    const std::vector<unsigned int>& aTri,
    const std::vector<int>& aBCFlag)
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
}

void LapDef_LinearDirectDisponly(
  std::vector<double>& aXYZ1,
  dfm2::CMatrixSparse<double>& mat_A,
  int iframe,
  const std::vector<double>& aXYZ0,
  const std::vector<int>& aBCFlag)
{
  SetLinSys_LaplaceGraph_MeshTri3(mat_A);
  std::vector<double> aDelta(aXYZ0.size());
  mat_A.MatVec(aDelta.data(),
               1.0, aXYZ0.data(), 0.0);
  // ----------
  std::vector<double> aRhs = aDelta;
  SetDisplacementAtFixedBoundary(aRhs,
                            iframe,aXYZ0,aBCFlag);
  mat_A.SetFixedBC_Dia(aBCFlag.data(), 1.0);
  mat_A.SetFixedBC_Row(aBCFlag.data());
  aXYZ1 = aXYZ0;
  std::vector<double> aRes = Solve_BiCGStab(aRhs, aXYZ1,
                                            1.0e-5, 100, mat_A);
  std::cout << aRes.size() << std::endl;
}

// --------------------------------------------------


class CDeformer_Laplace{
public:
  CDeformer_Laplace(const dfm2::CMatrixSparse<double>& A0,
       const std::vector<int>& aFlg0,
       double weight_bc0): A(A0), aBCFlag(aFlg0), weight_bc(weight_bc0){
    const unsigned int np = aFlg0.size()/3;
    vec_tmp.resize(np*3);
    
    // make jacobi preconditioner
    aDiaInv.assign(np*9,0.0);
    for(unsigned int ip=0;ip<np;++ip){
      for(unsigned int icrs=A.colInd[ip];icrs<A.colInd[ip+1];++icrs){
        unsigned int jp0 = A.rowPtr[icrs];
        dfm2::MatTMat3_ScaleAdd(aDiaInv.data()+jp0*9,
                                A.valCrs.data()+icrs*9,
                                A.valCrs.data()+icrs*9,
                                1.0, 0.0); // del. prev. value and set new vaue
      }
      {
        dfm2::MatTMat3_ScaleAdd(aDiaInv.data()+ip*9,
                                A.valDia.data()+ip*9,
                                A.valDia.data()+ip*9,
                                1.0, 0.0); // del. prev. value and set new vaue
      }
    }
    for(int ip=0;ip<np;++ip){
      for(int i=0;i<3;++i){
        if( aBCFlag[ip*3+i] == 0 ){ continue; }
        aDiaInv[ip*9+i*3+i] += weight_bc;
      }
    }
    for(unsigned int ip=0;ip<np;++ip){
      dfm2::Inverse_Mat3(aDiaInv.data()+ip*9);
    }
  }
  void MatVec(double* y,
              double alpha, const double* vec,  double beta) const {
    A.MatVec(vec_tmp.data(),
             1, vec, 0.0);
    A.MatTVec(y,
              alpha, vec_tmp.data(), beta);
    // add diagonal for fixed boundary condition
    for(int i=0;i<aBCFlag.size();++i){
      if( aBCFlag[i] == 0 ){ continue; }
      y[i] += weight_bc*vec[i];
    }
  }
  void Solve(double* v) const {
    const unsigned int np = aBCFlag.size()/3;
    for(int ip=0;ip<np;++ip){
      double tmp[3];
      dfm2::MatVec3(tmp, aDiaInv.data()+ip*9, v+ip*3);
      v[ip*3+0] = tmp[0];
      v[ip*3+1] = tmp[1];
      v[ip*3+2] = tmp[2];
    }
  }
public:
  const dfm2::CMatrixSparse<double>& A;
  const std::vector<int>& aBCFlag;
  const double weight_bc;
  std::vector<double> aDiaInv;
  mutable std::vector<double> vec_tmp;
};

void LapDef_LinearEnergyDisponly(
    std::vector<double>& aXYZ1,
    dfm2::CMatrixSparse<double>& mat_A,
    int iframe,
    bool is_preconditioner,
    const std::vector<double>& aXYZ0,
    const std::vector<unsigned int>& aTri,
    const std::vector<int>& aBCFlag)
{
  SetLinSys_LaplaceGraph_MeshTri3(mat_A);
  const double weight_bc = 100.0;
  CDeformer_Laplace mat_AtA(mat_A, aBCFlag, weight_bc);
  
  // ----------
  std::vector<double> aRhs(aXYZ0.size(),0.0);
  for(int i=0;i<aXYZ0.size();++i){ // adding noise for debuggng purpose
    aXYZ1[i] += 0.02*(double)rand()/(RAND_MAX+1.0)-0.01;
  }
  { // making RHS vector for elastic deformation
    const unsigned int np = aXYZ0.size()/3;
    std::vector<double> aTmp(np*3,0.0);
    for(unsigned int ip=0;ip<np;++ip){
      double* tmp = aTmp.data()+ip*3;
      for(unsigned int icrs=mat_A.colInd[ip];icrs<mat_A.colInd[ip+1];++icrs){
        unsigned int jp0 = mat_A.rowPtr[icrs];
        const double d0[3] = { aXYZ0[jp0*3+0]-aXYZ0[ip*3+0], aXYZ0[jp0*3+1]-aXYZ0[ip*3+1], aXYZ0[jp0*3+2]-aXYZ0[ip*3+2] };
        const double d1[3] = { aXYZ1[jp0*3+0]-aXYZ1[ip*3+0], aXYZ1[jp0*3+1]-aXYZ1[ip*3+1], aXYZ1[jp0*3+2]-aXYZ1[ip*3+2] };
        tmp[0] += +(d0[0] - d1[0]);
        tmp[1] += +(d0[1] - d1[1]);
        tmp[2] += +(d0[2] - d1[2]);
      }
    }
    mat_A.MatTVec(aRhs.data(),
                  -1.0, aTmp.data(), 0.0);
  }
  { // making RHS vector for fixed boundary condition
    const unsigned int np = aXYZ0.size()/3;
    std::vector<double> aGoal(np*3);
    SetDisplacementAtFixedBoundary(aGoal,
                              iframe,aXYZ0,aBCFlag);
    for(int i=0;i<np*3;++i){
      if( aBCFlag[i] == 0 ){ continue; }
      aRhs[i] += (aGoal[i]-aXYZ1[i])*weight_bc;
    }
  }
  std::vector<double> aUpd(aXYZ0.size(),0.0);
  if( is_preconditioner ){
    std::vector<double> aRes = dfm2::Solve_PCG(aRhs.data(), aUpd.data(),
                                               aRhs.size(), 1.0e-7, 300, mat_AtA, mat_AtA);
    std::cout << aRes.size() << std::endl;
  }
  else{
    std::vector<double> aRes = dfm2::Solve_CG(aRhs.data(), aUpd.data(),
                                              aRhs.size(), 1.0e-7, 300, mat_AtA);
    std::cout << aRes.size() << std::endl;
  }
  for(int i=0;i<aBCFlag.size();++i){ aXYZ1[i] += aUpd[i]; }
}

int main(int argc,char* argv[])
{
  std::vector<int> aBCFlag;
  std::vector<unsigned int> aTri;
  std::vector<double> aXYZ0;
  SetProblem(aTri, aXYZ0, aBCFlag);
  dfm2::CMatrixSparse<double> mat_A;
  {
    std::vector<unsigned int> psup_ind, psup;
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                               aTri.data(), aTri.size()/3, 3,
                               (int)aXYZ0.size()/3);
    dfm2::JArray_Sort(psup_ind, psup);
    mat_A.Initialize(aXYZ0.size()/3, 3, true);
    mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                     psup.data(),     psup.size());
  }
  std::vector<double> aXYZ1 = aXYZ0;
  
  // ------------------------
  
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.view_height = 1.0;
  viewer.nav.camera.camera_rot_mode = delfem2::CAMERA_ROT_TBALL;
  delfem2::opengl::setSomeLighting();
  
  int iframe = 0;
  while (!glfwWindowShouldClose(viewer.window))
  {
    glfwSetWindowTitle(viewer.window, "Direct Constraint");
    for(;iframe<100;++iframe){
      LapDef_LinearDirectDisponly(aXYZ1,mat_A,
                                  iframe, aXYZ0,aBCFlag);
      viewer.DrawBegin_oldGL();
      myGlutDisplay(aXYZ0,aXYZ1,aTri,aBCFlag);
      viewer.DrawEnd_oldGL();
      if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
    }
    glfwSetWindowTitle(viewer.window, "Energy-based Constraint without Preconditioner");
    for(;iframe<200;++iframe){
      LapDef_LinearEnergyDisponly(aXYZ1,mat_A,
                                  iframe,false, aXYZ0,aTri,aBCFlag);
      viewer.DrawBegin_oldGL();
      myGlutDisplay(aXYZ0,aXYZ1,aTri,aBCFlag);
      viewer.DrawEnd_oldGL();
      if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
    }
    glfwSetWindowTitle(viewer.window, "Energy-based Constraint with Preconditioner");
    for(;iframe<300;++iframe){
      LapDef_LinearEnergyDisponly(aXYZ1,mat_A,
                                  iframe,true, aXYZ0,aTri,aBCFlag);
      viewer.DrawBegin_oldGL();
      myGlutDisplay(aXYZ0,aXYZ1,aTri,aBCFlag);
      viewer.DrawEnd_oldGL();
      if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
    }
    // ------
    iframe = 0;
  }
  
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


