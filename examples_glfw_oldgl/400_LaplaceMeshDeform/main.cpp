#include <iostream>
#include <cmath>
#include <ctime>
#include "delfem2/primitive.h"
#include "delfem2/mshtopo.h"
#include "delfem2/mats.h"
#include "delfem2/emat.h"
#include "delfem2/fem_emats.h"
#include "delfem2/ilu_mats.h"
#include "delfem2/mshio.h"
#include "delfem2/mshmisc.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.hpp"
#include "delfem2/opengl/gl2_funcs.h"
#include "delfem2/opengl/gl2_color.h"
#include "delfem2/opengl/gl2_v23.h"

namespace dfm2 = delfem2;

// ------------------------------------

void FetchData(
    double* val_to,
    int nno, int ndim,
    const unsigned int* aIP,
    const double* val_from,
    int nstride=-1)
{
  if( nstride == -1 ){ nstride = ndim; }
  assert( nstride >= ndim );
  for(int ino=0;ino<nno;++ino){
    unsigned int ip = aIP[ino];
    for(int idim=0;idim<ndim;++idim){
      val_to[ino*ndim+idim] = val_from[ip*nstride+idim];
    }
  }
}

void SetLinSys_LaplaceGraph_MeshTri3D(
    dfm2::CMatrixSparse<double>& mat_A,
    std::vector<double>& vec_b,
    const double* aXYZ,
    unsigned int np,
    const unsigned int* aTri1,
    unsigned int nTri)
{
  const unsigned int nDoF = np;
  std::vector<int> tmp_buffer(nDoF, -1);
  mat_A.SetZero();
  for (unsigned int itri = 0; itri<nTri; ++itri){
    const unsigned int i0 = aTri1[itri*3+0];
    const unsigned int i1 = aTri1[itri*3+1];
    const unsigned int i2 = aTri1[itri*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][3]; FetchData(&coords[0][0],3,3,aIP, aXYZ);
    // ----------------------
    double emat[3][3][3][3];
    for(int i=0;i<81;++i){ (&emat[0][0][0][0])[i] = 0.0; }
    emat[0][1][0][0] = 1;
    emat[0][2][0][0] = 1;
    emat[1][2][0][0] = 1;
    emat[1][0][0][0] = 1;
    emat[2][0][0][0] = 1;
    emat[2][1][0][0] = 1;
    mat_A.Mearge(3, aIP, 3, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
  for(unsigned int ip=0;ip<np;++ip){
    double sum = 0.0;
    for(int icrs=mat_A.colInd[ip];icrs<mat_A.colInd[ip+1];++icrs){
      sum += mat_A.valCrs[icrs*9];
    }
    if( sum == 0.0 ) continue;
    double sum_inv = 1.0/sum;
    for(int icrs=mat_A.colInd[ip];icrs<mat_A.colInd[ip+1];++icrs){
      double v0 = -mat_A.valCrs[icrs*9] * sum_inv;
      mat_A.valCrs[icrs*9+0*3+0] = v0;
      mat_A.valCrs[icrs*9+1*3+1] = v0;
      mat_A.valCrs[icrs*9+2*3+2] = v0;
    }
    mat_A.valDia[ip*9+0*3+0] = 1;
    mat_A.valDia[ip*9+1*3+1] = 1;
    mat_A.valDia[ip*9+2*3+2] = 1;
  }
  // ----------------
  mat_A.MatVec(1.0, aXYZ, 0.0, vec_b.data());
  //std::cout << dfm2::CheckSymmetry(mat_A) << std::endl;
}

// ------------------------------------

std::vector<unsigned int> aTri;
std::vector<double> aXYZ;
std::vector<double> aXYZ1;
std::vector<double> aDelta;

dfm2::CMatrixSparse<double> mat_A;
dfm2::CPreconditionerILU<double>  ilu_A;
std::vector<int> aBCFlag;

// ------------------------------------------

void SetProblem()
{
//  dfm2::MeshTri3D_Cube(aXYZ, aTri, 10);
/*
  dfm2::MeshTri3D_CylinderOpen(aXYZ,aTri,
      0.2, 1.5,
      8,8);
      */

  /*
  const unsigned int nl = 8;
  dfm2::MeshTri3D_CylinderClosed(aXYZ, aTri,
                                 0.2, 1.5,
                                 8, nl);
                                 */
  {
    dfm2::Read_Obj(std::string(PATH_INPUT_DIR)+"/bunny_1k.obj",
        aXYZ,aTri);
    dfm2::Normalize_Points3D(aXYZ,1.0);
    dfm2::Rotate_Points3D(aXYZ, -M_PI*0.5, 0.0, 0.0);
  }

  {
    std::vector<unsigned int> psup_ind, psup;
    JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                                aTri.data(), aTri.size()/3, 3,
                                                (int)aXYZ.size()/3);
    dfm2::JArray_Sort(psup_ind, psup);
    mat_A.Initialize(aXYZ.size()/3, 3, true);
    mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                     psup.data(),     psup.size());
  }
  aDelta.resize(aXYZ.size());
  SetLinSys_LaplaceGraph_MeshTri3D(mat_A, aDelta,
                                   aXYZ.data(), aXYZ.size() / 3,
                                   aTri.data(), aTri.size() / 3);
  {
    const unsigned int np = aXYZ.size() / 3;
    aBCFlag.assign(np * 3, 0);
    aBCFlag[0] = 1;
    aBCFlag[1] = 1;
    aBCFlag[2] = 1;
    aBCFlag[(np - 1) * 3 + 0] = 1;
    aBCFlag[(np - 1) * 3 + 1] = 1;
    aBCFlag[(np - 1) * 3 + 2] = 1;
  }
}

// -------------------------------------------

// -------------------------------------------


void myGlutDisplay()
{
  ::glDisable(GL_LIGHTING);
  ::glColor3d(1,0,0);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ,aTri);
  ::glColor3d(0,0,0);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ.data(),aXYZ.size()/3, aTri.data(),aTri.size()/3);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ1.data(),aXYZ1.size()/3,aTri.data(),aTri.size()/3);
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
    {
      std::vector<double> aRhs = aDelta;
      const unsigned int np = aRhs.size()/3;
      aRhs[0*3+1] = -0.75;
      aRhs[(np-1)*3+1] = +1.0+cos(iframe*0.01);
      aRhs[(np-1)*3+0] = +sin(iframe*0.01);
      mat_A.SetFixedBC_Dia(aBCFlag.data());
      mat_A.SetFixedBC_Row(aBCFlag.data());
      aXYZ1 = aXYZ;
      std::vector<double> aRes = Solve_BiCGSTAB(aRhs,
                                                aXYZ1, 1.0e-5, 100, mat_A);
      std::cout << aRes.size() << std::endl;
    }
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


