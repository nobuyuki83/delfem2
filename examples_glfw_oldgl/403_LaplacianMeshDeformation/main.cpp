/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <cmath>
#include "delfem2/mshtopo.h"
#include "delfem2/mat3.h"
#include "delfem2/mats.h"
#include "delfem2/emat.h"
#include "delfem2/fem_emats.h"
#include "delfem2/ilu_mats.h"
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




// ------------------------------------

std::vector<unsigned int> aTri;
std::vector<double> aXYZ0;
std::vector<double> aXYZ1;
std::vector<double> aDelta;

dfm2::CMatrixSparse<double> mat_A;
std::vector<int> aBCFlag;

// ------------------------------------------

void SetProblem()
{
  dfm2::MeshTri3D_CylinderClosed(aXYZ0, aTri,
                                 0.2, 1.6,
                                 16, 16);
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
  // ---------------------------------
  SetLinSys_LaplaceGraph_MeshTri3(mat_A);
  aDelta.resize(aXYZ0.size());
  mat_A.MatVec(aDelta.data(),
               1.0, aXYZ0.data(), 0.0);
  mat_A.SetFixedBC_Dia(aBCFlag.data());
  mat_A.SetFixedBC_Row(aBCFlag.data());
}

void SetFixedBoundaryCondition(std::vector<double>& aRhs,
                               unsigned int iframe)
{
  double A[16];
  {
    dfm2::AffMat3_Identity(A);
    const double trans0[3] = {0, -0.8, 0};
    dfm2::Translate_AffMat3(A,
                            trans0);
    const double axis0[3] = {0, +0.0, 0.5*sin(0.1*iframe)};
    dfm2::Rotate_AffMat3_Rodriguez(A,
                                   axis0);
    const double trans1[3] = {0.2*sin(0.03*iframe), +0.6-0.2*cos(0.05*iframe), 0};
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


// -------------------------------------------
// -------------------------------------------

void myGlutDisplay()
{
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1,aTri);
  ::glColor3d(0,0,0);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ0.data(),aXYZ0.size()/3,
                                   aTri.data(),aTri.size()/3);
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
    std::vector<double> aRhs = aDelta;
    SetFixedBoundaryCondition(aRhs,iframe);
    aXYZ1 = aXYZ0;
    std::vector<double> aRes = Solve_BiCGSTAB(aRhs,
                                              aXYZ1, 1.0e-5, 100, mat_A);
    std::cout << aRes.size() << std::endl;

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


