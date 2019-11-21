/**
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file fem_emats.cpp
 * @brief implementation of the merge functions for various PDE
 * @author Nobuyuki Umetani
 * @date 2018
 * @details this file only depends on and "emat.h" and "mats.h"
 */

#include <cstdio>
#include <complex>

typedef std::complex<double> COMPLEX;

#include "delfem2/emat.h"
#include "delfem2/mats.h"

#include "delfem2/fem_emats.h"

namespace dfm2 = delfem2;

static void FetchData
(double* val_to,
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


// area of a triangle
static double TriArea2D(const double p0[], const double p1[], const double p2[]){
  return 0.5*((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]));
}

static double TetVolume3D
(const double v1[3],
 const double v2[3],
 const double v3[3],
 const double v4[3])
{
  return
  ((v2[0]-v1[0])*((v3[1]-v1[1])*(v4[2]-v1[2])-(v4[1]-v1[1])*(v3[2]-v1[2]))
   -(v2[1]-v1[1])*((v3[0]-v1[0])*(v4[2]-v1[2])-(v4[0]-v1[0])*(v3[2]-v1[2]))
   +(v2[2]-v1[2])*((v3[0]-v1[0])*(v4[1]-v1[1])-(v4[0]-v1[0])*(v3[1]-v1[1]))
   ) * 0.16666666666666666666666666666667;
}


// caluculate Derivative of Area Coord
static inline void TetDlDx(double dldx[][3], double a[],
                           const double p0[], const double p1[], const double p2[], const double p3[])
{
  const double vol = TetVolume3D(p0, p1, p2, p3);
  const double dtmp1 = 1.0/(vol * 6.0);
  
  a[0] = +dtmp1*(p1[0]*(p2[1]*p3[2]-p3[1]*p2[2])-p1[1]*(p2[0]*p3[2]-p3[0]*p2[2])+p1[2]*(p2[0]*p3[1]-p3[0]*p2[1]));
  a[1] = -dtmp1*(p2[0]*(p3[1]*p0[2]-p0[1]*p3[2])-p2[1]*(p3[0]*p0[2]-p0[0]*p3[2])+p2[2]*(p3[0]*p0[1]-p0[0]*p3[1]));
  a[2] = +dtmp1*(p3[0]*(p0[1]*p1[2]-p1[1]*p0[2])-p3[1]*(p0[0]*p1[2]-p1[0]*p0[2])+p3[2]*(p0[0]*p1[1]-p1[0]*p0[1]));
  a[3] = -dtmp1*(p0[0]*(p1[1]*p2[2]-p2[1]*p1[2])-p0[1]*(p1[0]*p2[2]-p2[0]*p1[2])+p0[2]*(p1[0]*p2[1]-p2[0]*p1[1]));
  
  dldx[0][0] = -dtmp1*((p2[1]-p1[1])*(p3[2]-p1[2])-(p3[1]-p1[1])*(p2[2]-p1[2]));
  dldx[0][1] = +dtmp1*((p2[0]-p1[0])*(p3[2]-p1[2])-(p3[0]-p1[0])*(p2[2]-p1[2]));
  dldx[0][2] = -dtmp1*((p2[0]-p1[0])*(p3[1]-p1[1])-(p3[0]-p1[0])*(p2[1]-p1[1]));
  
  dldx[1][0] = +dtmp1*((p3[1]-p2[1])*(p0[2]-p2[2])-(p0[1]-p2[1])*(p3[2]-p2[2]));
  dldx[1][1] = -dtmp1*((p3[0]-p2[0])*(p0[2]-p2[2])-(p0[0]-p2[0])*(p3[2]-p2[2]));
  dldx[1][2] = +dtmp1*((p3[0]-p2[0])*(p0[1]-p2[1])-(p0[0]-p2[0])*(p3[1]-p2[1]));
  
  dldx[2][0] = -dtmp1*((p0[1]-p3[1])*(p1[2]-p3[2])-(p1[1]-p3[1])*(p0[2]-p3[2]));
  dldx[2][1] = +dtmp1*((p0[0]-p3[0])*(p1[2]-p3[2])-(p1[0]-p3[0])*(p0[2]-p3[2]));
  dldx[2][2] = -dtmp1*((p0[0]-p3[0])*(p1[1]-p3[1])-(p1[0]-p3[0])*(p0[1]-p3[1]));
  
  dldx[3][0] = +dtmp1*((p1[1]-p0[1])*(p2[2]-p0[2])-(p2[1]-p0[1])*(p1[2]-p0[2]));
  dldx[3][1] = -dtmp1*((p1[0]-p0[0])*(p2[2]-p0[2])-(p2[0]-p0[0])*(p1[2]-p0[2]));
  dldx[3][2] = +dtmp1*((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]));
  
  //  std::cout << dldx[0][0]+dldx[1][0]+dldx[2][0]+dldx[3][0] << std::endl;
  //  std::cout << dldx[0][1]+dldx[1][1]+dldx[2][1]+dldx[3][1] << std::endl;
  //  std::cout << dldx[0][2]+dldx[1][2]+dldx[2][2]+dldx[3][2] << std::endl;
  
  //  std::cout << a[0]+dldx[0][0]*p0[0]+dldx[0][1]*p0[1]+dldx[0][2]*p0[2] << std::endl;
  //  std::cout << a[1]+dldx[1][0]*p1[0]+dldx[1][1]*p1[1]+dldx[1][2]*p1[2] << std::endl;
  //  std::cout << a[2]+dldx[2][0]*p2[0]+dldx[2][1]*p2[1]+dldx[2][2]*p2[2] << std::endl;
  //  std::cout << a[3]+dldx[3][0]*p3[0]+dldx[3][1]*p3[1]+dldx[3][2]*p3[2] << std::endl;
}


static void MatVec3
(double y[3],
 const double m[9], const double x[3]){
  y[0] = m[0]*x[0] + m[1]*x[1] + m[2]*x[2];
  y[1] = m[3]*x[0] + m[4]*x[1] + m[5]*x[2];
  y[2] = m[6]*x[0] + m[7]*x[1] + m[8]*x[2];
}

static void MatMat3
(double* C,
 const double* A, const double* B)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      C[i*3+j] = A[i*3+0]*B[0*3+j] + A[i*3+1]*B[1*3+j] + A[i*3+2]*B[2*3+j];
    }
  }
}

static void MatMatTrans3
(double* C,
 const double* A, const double* B)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      C[i*3+j] = A[i*3+0]*B[j*3+0] + A[i*3+1]*B[j*3+1] + A[i*3+2]*B[j*3+2];
    }
  }
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void dfm2::MergeLinSys_Poission_MeshTri2D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double alpha,
 const double source,
 const double* aXY1, int np,
 const unsigned int* aTri1, int nTri,
 const double* aVal)
{
  const int nDoF = np;
  /////
  std::vector<int> tmp_buffer(nDoF, -1);
  for (int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1);
    const double value[3] = { aVal[i0], aVal[i1], aVal[i2] };
    ////
    double eres[3];
    double emat[3][3];
    EMat_Poisson_Tri2D
    (eres,emat,
     alpha, source,
     coords, value);
    for (int ino = 0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    mat_A.Mearge(3, aIP, 3, aIP, 1, &emat[0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_Helmholtz_MeshTri2D
(CMatrixSparse<COMPLEX>& mat_A,
 COMPLEX* vec_b,
 const double wave_length,
 const double* aXY1, int np,
 const unsigned int* aTri1, int nTri,
 const COMPLEX* aVal)
{
  const int nDoF = np;
  std::vector<int> tmp_buffer(nDoF, -1);
  for (int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1);
    const COMPLEX value[3] = { aVal[i0], aVal[i1], aVal[i2] };
    ////
    std::complex<double> eres[3];
    std::complex<double> emat[3][3];
    EMat_Helmholtz_Tri2D(eres,emat,
                         wave_length,
                         coords, value);
    for(int ino=0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    mat_A.Mearge(3, aIP, 3, aIP, 1, &emat[0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_SommerfeltRadiationBC_Polyline2D
(CMatrixSparse<COMPLEX>& mat_A,
 COMPLEX* vec_b,
 const double wave_length,
 const double* aXY1, int np,
 const unsigned int* aIP, int nIP,
 const COMPLEX* aVal)
{
  const int nDoF = np;
  std::vector<int> tmp_buffer(nDoF, -1);
  for(int iel=0; iel<nIP-1; ++iel){
    const unsigned int i0 = aIP[iel+0];
    const unsigned int i1 = aIP[iel+1];
    const unsigned int aIP[2] = {i0,i1};
    double P[2][2]; FetchData(&P[0][0],2,2,aIP, aXY1);
    const COMPLEX val[2] = { aVal[i0], aVal[i1] };
    ////
    COMPLEX eres[2], emat[2][2];
    EMat_SommerfeltRadiationBC_Line2D(eres,emat,
                                         wave_length,P,val);
    for(int ino=0;ino<2;ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    mat_A.Mearge(2, aIP, 2, aIP, 1, &emat[0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_Poission_MeshTet3D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double alpha,
 const double source,
 const double* aXYZ, int nXYZ,
 const unsigned int* aTet, int nTet,
 const double* aVal)
{
  const int np = nXYZ;
  std::vector<int> tmp_buffer(np, -1);
  for (int itet = 0; itet<nTet; ++itet){
    const unsigned int i0 = aTet[itet*4+0];
    const unsigned int i1 = aTet[itet*4+1];
    const unsigned int i2 = aTet[itet*4+2];
    const unsigned int i3 = aTet[itet*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; FetchData(&coords[0][0],4,3,aIP, aXYZ);
    const double value[4] = { aVal[i0], aVal[i1], aVal[i2], aVal[i3] };
    ////
    double eres[4], emat[4][4];
    EMat_Poisson_Tet3D(eres,emat,
                       alpha, source,
                       coords, value);
    for (int ino = 0; ino<4; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    mat_A.Mearge(4, aIP, 4, aIP, 1, &emat[0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_Diffusion_MeshTri2D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double alpha,
 const double rho,
 const double source,
 const double dt_timestep,
 const double gamma_newmark,
 const double* aXY1, int nXY,
 const unsigned int* aTri1, int nTri,
 const double* aVal,
 const double* aVelo)
{
//  const int nDoF = nXY;
  ////
//  mat_A.SetZero();
//  for(int idof=0;idof<nDoF;++idof){ vec_b[idof] = 0.0; }
  std::vector<int> tmp_buffer(nXY, -1);
  for (int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1);
    const double value[3] = { aVal[ i0], aVal[ i1], aVal[ i2] };
    const double velo[ 3] = { aVelo[i0], aVelo[i1], aVelo[i2] };
    ////
    double eres[3];
    double emat[3][3];
    EMat_Diffusion_Tri2D
    (eres,emat,
     alpha, source,
     dt_timestep, gamma_newmark, rho,
     coords, value, velo);
    for (int ino = 0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    mat_A.Mearge(3, aIP, 3, aIP, 1, &emat[0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_Diffusion_MeshTet3D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double alpha,
 const double rho,
 const double source,
 const double dt_timestep,
 const double gamma_newmark,
 const double* aXYZ, int nXYZ,
 const unsigned int* aTet, int nTet,
 const double* aVal,
 const double* aVelo)
{
  const int np = nXYZ;
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTet; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; FetchData(&coords[0][0],4,3,aIP, aXYZ);
    const double value[4] = { aVal[ i0], aVal[ i1], aVal[ i2], aVal[ i3] };
    const double velo[ 4] = { aVelo[i0], aVelo[i1], aVelo[i2], aVelo[i3] };
    ////
    double eres[4];
    double emat[4][4];
    EMat_Diffusion_Newmark_Tet3D(eres,emat,
                                 alpha, source,
                                 dt_timestep, gamma_newmark, rho,
                                 coords, value, velo);
    for (int ino = 0; ino<4; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    mat_A.Mearge(4, aIP, 4, aIP, 1, &emat[0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_SolidLinear_Static_MeshTri2D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double myu,
 const double lambda,
 const double rho,
 const double g_x,
 const double g_y,
 const double* aXY1, int nXY,
 const unsigned int* aTri1, int nTri,
 const double* aVal)
{
  const int np = nXY;
  std::vector<int> tmp_buffer(np, -1);
  for(int iel=0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1);
    double disps[3][2]; FetchData(&disps[0][0],3,2,aIP, aVal);
    ////
    double eres[3][2];
    double emat[3][3][2][2];
    EMat_SolidStaticLinear_Tri2D(eres,emat,
                                 myu, lambda, rho, g_x, g_y,
                                 disps, coords);
    for (int ino = 0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*2+0] += eres[ino][0];
      vec_b[ip*2+1] += eres[ino][1];
    }
    mat_A.Mearge(3, aIP, 3, aIP, 4, &emat[0][0][0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_SolidLinear_NewmarkBeta_MeshTri2D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double myu,
 const double lambda,
 const double rho,
 const double g_x,
 const double g_y,
 const double dt_timestep,
 const double gamma_newmark,
 const double beta_newmark,
 const double* aXY1, int nXY,
 const unsigned int* aTri1, int nTri,
 const double* aVal,
 const double* aVelo,
 const double* aAcc)
{
  const int np = nXY;
//  const int nDoF = np*2;
  ////
//  mat_A.SetZero();
//  for(int idof=0;idof<nDoF;idof++){ vec_b[idof] = 0.0; }
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1);
    double disps[3][2]; FetchData(&disps[0][0],3,2,aIP, aVal);
    double velos[3][2]; FetchData(&velos[0][0],3,2,aIP, aVelo);
    double accs[3][2]; FetchData(&accs[0][0],3,2,aIP, aAcc);
    ////
    double eres[3][2];
    double emat[3][3][2][2];
    EMat_SolidDynamicLinear_Tri2D
    (eres,emat,
     myu, lambda,
     rho, g_x, g_y,
     dt_timestep, gamma_newmark, beta_newmark,
     disps, velos, accs, coords, 
     true);
    for (int ino = 0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*2+0] += eres[ino][0];
      vec_b[ip*2+1] += eres[ino][1];
    }
    // marge dde
    mat_A.Mearge(3, aIP, 3, aIP, 4, &emat[0][0][0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_StokesStatic2D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double myu,
 const double g_x,
 const double g_y,
 const double* aXY1, int nXY,
 const unsigned int* aTri1, int nTri,
 const double* aVal)
{
  const int np = nXY;
//  const int nDoF = np*3;
  ////
//  mat_A.SetZero();
//  for(int idof=0;idof<nDoF;++idof){ vec_b[idof] = 0.0; }
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1);
    double velo_press[3][3]; FetchData(&velo_press[0][0],3,3,aIP, aVal);
    ////
    double eres[3][3];
    double emat[3][3][3][3];
    ////
    EMat_Stokes2D_Static_P1(myu, g_x, g_y, coords, velo_press, emat, eres);
    for (int ino = 0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0];
      vec_b[ip*3+1] += eres[ino][1];
      vec_b[ip*3+2] += eres[ino][2];
    }
    // marge dde
    mat_A.Mearge(3, aIP, 3, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_StokesDynamic2D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double myu,
 const double rho,
 const double g_x,
 const double g_y,
 const double dt_timestep,
 const double gamma_newmark,
 const double* aXY1, int nXY,
 const unsigned int* aTri1, int nTri,
 const double* aVal,
 const double* aVelo)
{
  const int np = nXY;
//  const int nDoF = np*3;
  ////
//  mat_A.SetZero();
//  for(int i=0;i<nDoF;++i){ vec_b[i] = 0.0; }
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1);
    double velo_press[3][3]; FetchData(&velo_press[0][0],3,3,aIP, aVal);
    double acc_apress[3][3]; FetchData(&acc_apress[0][0],3,3,aIP, aVelo);
    ////
    double eres[3][3];
    double emat[3][3][3][3];
    ////
    EMat_Stokes2D_Dynamic_P1(myu, rho,  g_x, g_y,
                                dt_timestep, gamma_newmark,
                                coords, velo_press, acc_apress,
                                emat, eres);
    for (int ino = 0; ino<3; ino++){
      const int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0];
      vec_b[ip*3+1] += eres[ino][1];
      vec_b[ip*3+2] += eres[ino][2];
    }
    mat_A.Mearge(3, aIP, 3, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_NavierStokes2D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double myu,
 const double rho,
 const double g_x,
 const double g_y,
 const double dt_timestep,
 const double gamma_newmark,
 const double* aXY1, int nXY,
 const unsigned int* aTri1, int nTri,
 const double* aVal, // vx,vy,press
 const double* aDtVal) // ax,ay,apress
{
  const int np = nXY;
//  const int nDoF = np*3;
  ////
//  mat_A.SetZero();
//  for(int i=0;i<nDoF;++i){ vec_b[i] = 0.0; }
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1);
    double velo[3][3]; FetchData(&velo[0][0],3,3,aIP, aVal);
    double acc[3][3]; FetchData(&acc[0][0],3,3,aIP, aDtVal);
    ////
    double eres[3][3], emat[3][3][3][3];
    EMat_NavierStokes2D_Dynamic_P1(myu, rho,  g_x, g_y,
                                      dt_timestep, gamma_newmark,
                                      coords, velo, acc,
                                      emat, eres);
    for (int ino = 0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0];
      vec_b[ip*3+1] += eres[ino][1];
      vec_b[ip*3+2] += eres[ino][2];
    }
    mat_A.Mearge(3, aIP, 3, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
}


// compute total energy and its first and second derivatives
double dfm2::MergeLinSys_Cloth
(CMatrixSparse<double>& ddW, // (out) second derivative of energy
 double* dW, // (out) first derivative of energy
 ////
 double lambda, // (in) Lame's 1st parameter
 double myu,  // (in) Lame's 2nd parameter
 double stiff_bend, // (in) bending stiffness
 const double* aPosIni, int np, int ndim,
 const unsigned int* aTri, int nTri, // (in) triangle index
 const unsigned int* aQuad, int nQuad, // (in) index of 4 vertices required for bending
 const double* aXYZ
 )
{
  double W = 0;
  std::vector<int> tmp_buffer(np,-1);
  
  // marge element in-plane strain energy
  for(int itri=0;itri<nTri;itri++){
    const unsigned int aIP[3] = { aTri[itri*3+0], aTri[itri*3+1], aTri[itri*3+2] };
    double C[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double c[3][3];
    for(int ino=0;ino<3;ino++){
      const int ip = aIP[ino];
      for(int i=0;i<ndim;i++){ C[ino][i] = aPosIni[ip*ndim+i]; }
      for(int i=0;i<3;i++){ c[ino][i] = aXYZ[ip*3+i]; }
    }
    double e, de[3][3], dde[3][3][3][3];
    WdWddW_CST( e,de,dde, C,c, lambda,myu );
    W += e;  // marge energy
    // marge de
    for(int ino=0;ino<3;ino++){
      const unsigned int ip = aIP[ino];
      for(int i =0;i<3;i++){ dW[ip*3+i] += de[ino][i]; }
    }
    // marge dde
    ddW.Mearge(3, aIP, 3, aIP, 9, &dde[0][0][0][0], tmp_buffer);
  }
//  std::cout << "cst:" << W << std::endl;
  // marge element bending energy
  for(int iq=0;iq<nQuad;iq++){
    const unsigned int aIP[4] = { aQuad[iq*4+0], aQuad[iq*4+1], aQuad[iq*4+2], aQuad[iq*4+3] };
    double C[4][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
    double c[4][3];
    for(int ino=0;ino<4;ino++){
      const unsigned int ip = aIP[ino];
      for(int i=0;i<ndim;i++){ C[ino][i] = aPosIni[ip*ndim+i]; }
      for(int i=0;i<3;i++){ c[ino][i] = aXYZ [ip*3+i]; }
    }
    double e, de[4][3], dde[4][4][3][3];
    WdWddW_Bend( e,de,dde, C,c, stiff_bend );
    W += e;  // marge energy
    // marge de
    for(int ino=0;ino<4;ino++){
      const int ip = aIP[ino];
      for(int i =0;i<3;i++){ dW[ip*3+i] += de[ino][i]; }
    }
    // marge dde
    ddW.Mearge(4, aIP, 4, aIP, 9, &dde[0][0][0][0], tmp_buffer);
  }
  return W;
}




double dfm2::MergeLinSys_Contact
(CMatrixSparse<double>& ddW,
 double* dW,
 ////
 double stiff_contact,
 double contact_clearance,
 const CInput_Contact& input,
 const double* aXYZ,
 int nXYZ)
{
  const unsigned int np = nXYZ;
  std::vector<int> tmp_buffer(np,-1);
  double W = 0;
  for(unsigned int ip=0;ip<np;ip++){
    double c[3] = { aXYZ[ip*3+0], aXYZ[ip*3+1], aXYZ[ip*3+2] };
    double e, de[3], dde[3][3];
    WdWddW_Contact( e,de,dde, c, stiff_contact,contact_clearance, input );
    W += e;  // marge energy
    // marge de
    for(int i =0;i<3;i++){ dW[ip*3+i] += de[i]; }
    // marge dde
    ddW.Mearge(1, &ip, 1, &ip, 9, &dde[0][0], tmp_buffer);
  }
  return W;
}

/*
 void Solve_LinearSolid_TetP1()
 {
 
 unsigned int ndof = aXYZ.size();
 MatrixXd Mat(ndof,ndof);
 //  Eigen::SparseMatrix<double> Mat(ndof,ndof);
 Mat.setZero();
 
 VectorXd Lhs(ndof);
 Lhs.setZero();
 
 for(unsigned int itet=0;itet<aTet.size()/4;itet++){
 const unsigned int aPo[4] = { aTet[itet*4+0], aTet[itet*4+1], aTet[itet*4+2], aTet[itet*4+3] };
 unsigned int i0=aPo[0], i1=aPo[1], i2=aPo[2], i3=aPo[3];
 double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
 double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
 double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]};
 double p3[3] = {aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2]};
 ////
 const double vol = TetVolume3D(p0,p1,p2,p3);
 double dldx[4][3];
 double zero_order_term[4];
 TetDlDx(dldx, zero_order_term,   p0,p1,p2,p3);
 ////
 double disp[4][3] = { {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0} };
 double emat[4][4][3][3], eres[4][3];
 matRes_LinearSolid_TetP1(emat,eres,
 vol,lambda,myu,
 0,1,0, 1,
 dldx,disp);
 
 for(unsigned int ino=0;ino<4;ino++){
 for(unsigned int jno=0;jno<4;jno++){
 unsigned int ipo = aPo[ino];
 unsigned int jpo = aPo[jno];
 for(unsigned int idim=0;idim<3;idim++){
 for(unsigned int jdim=0;jdim<3;jdim++){
 Mat(ipo*3+idim,jpo*3+jdim) += emat[ino][jno][idim][jdim];
 //            Mat.coeffRef(ipo*3+idim,jpo*3+jdim) += emat[ino][jno][idim][jdim];
 }
 }
 }
 }
 for(unsigned int ino=0;ino<4;ino++){
 unsigned int ipo = aPo[ino];
 for(unsigned int idim=0;idim<3;idim++){
 Lhs(ipo*3+idim) += eres[ino][idim];
 }
 }
 }
 
 std::vector<unsigned int> aFix;
 {
 for(unsigned int ipo=0;ipo<aXYZ.size()/3;ipo++){
 double p[3] = {aXYZ[ipo*3+0], aXYZ[ipo*3+1], aXYZ[ipo*3+2]};
 if( p[0] < 0 ){
 aFix.push_back(ipo);
 }
 }
 }
 for(unsigned int ifix=0;ifix<aFix.size();ifix++){
 unsigned int ifix0 = aFix[ifix];
 for(int i=0;i<3;i++){
 Lhs(ifix0*3+i) = 0;
 }
 for(unsigned int jpo=0;jpo<aXYZ.size()/3;jpo++){
 for(int i=0;i<3;i++){
 for(int j=0;j<3;j++){
 Mat(ifix0*3+i,jpo*3+j) = 0;
 //          Mat.coeffRef(ifix0*3+i,jpo*3+j) = 0;
 }
 }
 }
 for(unsigned int jpo=0;jpo<aXYZ.size()/3;jpo++){
 for(int i=0;i<3;i++){
 for(int j=0;j<3;j++){
 Mat(jpo*3+i,ifix0*3+j) = 0;
 //          Mat.coeffRef(jpo*3+i,ifix0*3+j) = 0;
 }
 }
 }
 {
 for(int i=0;i<3;i++){
 Mat(ifix0*3+i,ifix0*3+i) = 1;
 //        Mat.coeffRef(ifix0*3+i,ifix0*3+i) = 1;
 }
 }
 }
 
 
 //  sMat = Mat;
 //  Mat.makeCompressed();
 
 VectorXd Rhs;
 Rhs = Mat.partialPivLu().solve(Lhs);
 
 //  std::cout << ndof << std::endl;
 //  ConjugateGradient<SparseMatrix<double>> cg;
 //  cg.setMaxIterations(10000);
 //  cg.compute(Mat);
 //  Rhs = cg.solve(Lhs);
 
 for(unsigned int ip=0;ip<aXYZ.size()/3;ip++){
 for(unsigned int i=0;i<3;i++){
 aDisp[ip*3+i] += Rhs(ip*3+i);
 }
 }
 }
 */


void dfm2::MergeLinSys_SolidLinear_Static_MeshTet3D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double myu,
 const double lambda,
 const double rho,
 const double g[3],
 const double* aXYZ, int nXYZ,
 const unsigned int* aTet, int nTet,
 const double* aDisp)
{
  const int np = nXYZ;
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTet; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = { i0, i1, i2, i3 };
    double P[4][3]; FetchData(&P[0][0], 4, 3, aIP, aXYZ);
    double disps[4][3]; FetchData(&disps[0][0], 4, 3, aIP, aDisp);
    ////
    double emat[4][4][3][3];
    for(int i=0;i<144;++i){ (&emat[0][0][0][0])[i] = 0.0; } // zero-clear
    double eres[4][3];
    {
      const double vol = TetVolume3D(P[0],P[1],P[2],P[3]);
      for(int ino=0;ino<4;++ino){
        eres[ino][0] = vol*rho*g[0]*0.25;
        eres[ino][1] = vol*rho*g[1]*0.25;
        eres[ino][2] = vol*rho*g[2]*0.25;
      }
    }
    EMat_SolidLinear_Static_Tet(emat,eres,
                                myu, lambda,
                                P, disps,
                                true); // additive
    for (int ino = 0; ino<4; ino++){
      const int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0];
      vec_b[ip*3+1] += eres[ino][1];
      vec_b[ip*3+2] += eres[ino][2];
    }
    mat_A.Mearge(4, aIP, 4, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_LinearSolid3D_Static_Q1
(CMatrixSparse<double>& mat_A,
 std::vector<double>& vec_b,
 const double myu,
 const double lambda,
 const double rho,
 const double g_x,
 const double g_y,
 const double g_z,
 const std::vector<double>& aXYZ,
 const std::vector<int>& aHex,
 const std::vector<double>& aVal)
{
  const int np = (int)aXYZ.size()/3;
  const int nDoF = np*3;
  ////
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<(int)aHex.size()/8; ++iel){
    const unsigned int i0 = aHex[iel*8+0];
    const unsigned int i1 = aHex[iel*8+1];
    const unsigned int i2 = aHex[iel*8+2];
    const unsigned int i3 = aHex[iel*8+3];
    const unsigned int i4 = aHex[iel*8+4];
    const unsigned int i5 = aHex[iel*8+5];
    const unsigned int i6 = aHex[iel*8+6];
    const unsigned int i7 = aHex[iel*8+7];
    const unsigned int aIP[8] = { i0, i1, i2, i3, i4, i5, i6, i7 };
    double coords[8][3]; FetchData(&coords[0][0], 8, 3, aIP, aXYZ.data());
    double disps[8][3]; FetchData(&disps[0][0], 8, 3, aIP, aVal.data());
    ////
    double eres[8][3];
    double emat[8][8][3][3];
    MakeMat_LinearSolid3D_Static_Q1(myu, lambda,
                                    rho, g_x, g_y, g_z,
                                    coords, disps,
                                    emat,eres);
    for (int ino = 0; ino<8; ino++){
      const int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0];
      vec_b[ip*3+1] += eres[ino][1];
      vec_b[ip*3+2] += eres[ino][2];
    }
    // marge dde
    mat_A.Mearge(8, aIP, 8, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_SolidLinear_NewmarkBeta_MeshTet3D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double myu,
 const double lambda,
 const double rho,
 const double g[3],
 const double dt_timestep,
 const double gamma_newmark,
 const double beta_newmark,
 const double* aXYZ, int nXYZ,
 const unsigned int* aTet, int nTet,
 const double* aVal,
 const double* aVelo,
 const double* aAcc)
{
  const int np = nXYZ;
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTet; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double P[4][3]; FetchData(&P[0][0],4,3,aIP, aXYZ);
    double disps[4][3];  FetchData(&disps[0][0], 4,3,aIP, aVal);
    double velos[4][3];  FetchData(&velos[0][0], 4,3,aIP, aVelo);
    double accs[4][3];   FetchData(&accs[0][0],  4,3,aIP, aAcc);
    ////
    double eres[4][3], emat[4][4][3][3];
    EMat_SolidLinear_NewmarkBeta_MeshTet3D(eres,emat,
                                                  myu, lambda,
                                                  rho, g[0], g[1], g[2],
                                                  dt_timestep, gamma_newmark, beta_newmark,
                                                  disps, velos, accs, P,
                                                  true);
    for (int ino = 0; ino<4; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0];
      vec_b[ip*3+1] += eres[ino][1];
      vec_b[ip*3+2] += eres[ino][2];
    }
    mat_A.Mearge(4, aIP, 4, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_SolidLinear_BEuler_MeshTet3D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double myu,
 const double lambda,
 const double rho,
 const double g[3],
 const double dt,
 const double* aXYZ, int nXYZ,
 const unsigned int* aTet, int nTet,
 const double* aDisp,
 const double* aVelo)
{
  const int np = nXYZ;
  std::vector<int> tmp_buffer(np, -1);
  for(int iel=0; iel<nTet; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = { i0, i1, i2, i3 };
    double P[4][3]; FetchData(&P[0][0], 4, 3, aIP, aXYZ);
    double emat[4][4][3][3];
    double eres[4][3];
    const double vol = TetVolume3D(P[0], P[1], P[2], P[3]);
    {
      double dldx[4][3], const_term[4];
      TetDlDx(dldx, const_term, P[0], P[1], P[2], P[3]);
      ddW_SolidLinear_Tet3D(&emat[0][0][0][0],
                            lambda, myu, vol, dldx, false, 3);
    }
    {
      double u[4][3]; FetchData(&u[0][0], 4, 3, aIP, aDisp);
      double v[4][3]; FetchData(&v[0][0], 4, 3, aIP, aVelo);
      for(int ino=0;ino<4;++ino){
        for(int idim=0;idim<3;++idim){
          eres[ino][idim] = vol*rho*g[idim]*0.25;
          for(int jno=0;jno<4;++jno){
            eres[ino][idim] -= emat[ino][jno][idim][0]*(u[jno][0]+dt*v[jno][0]);
            eres[ino][idim] -= emat[ino][jno][idim][1]*(u[jno][1]+dt*v[jno][1]);
            eres[ino][idim] -= emat[ino][jno][idim][2]*(u[jno][2]+dt*v[jno][2]);
          }
        }
      }
    }
    {
      for(int ino=0;ino<4;++ino){
        emat[ino][ino][0][0] += rho*vol*0.25/(dt*dt);
        emat[ino][ino][1][1] += rho*vol*0.25/(dt*dt);
        emat[ino][ino][2][2] += rho*vol*0.25/(dt*dt);
      }
    }
    for (int ino = 0; ino<4; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0]/dt;
      vec_b[ip*3+1] += eres[ino][1]/dt;
      vec_b[ip*3+2] += eres[ino][2]/dt;
    }
    mat_A.Mearge(4, aIP, 4, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
}


void dfm2::MergeLinSys_SolidStiffwarp_BEuler_MeshTet3D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double myu,
 const double lambda,
 const double rho,
 const double g[3],
 const double dt,
 const double* aXYZ, int nXYZ,
 const unsigned int* aTet, int nTet,
 const double* aDisp,
 const double* aVelo,
 const std::vector<double>& aR)
{
  const int np = nXYZ;
  assert((int)aR.size()==np*9);
  //////
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTet; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = { i0, i1, i2, i3 };
    double P[4][3]; FetchData(&P[0][0], 4, 3, aIP, aXYZ);
    const double vol = TetVolume3D(P[0], P[1], P[2], P[3]);
    ////
    double emat[4][4][3][3];
    { // make stifness matrix with stiffness warping
      double dldx[4][3], const_term[4];
      TetDlDx(dldx, const_term, P[0], P[1], P[2], P[3]);
      double emat0[4][4][3][3];
      ddW_SolidLinear_Tet3D(&emat0[0][0][0][0],
                            lambda, myu, vol, dldx, false, 3);
      double mtmp[9];
      for(int ino=0;ino<4;++ino){
        const double* Mi = aR.data()+aIP[ino]*9;
        for(int jno=0;jno<4;++jno){
          MatMatTrans3(mtmp, &emat0[ino][jno][0][0], Mi);
          MatMat3(&emat[ino][jno][0][0], Mi,mtmp);
        }
      }
    }
    double eres[4][3];
    {
      for(int ino=0;ino<4;++ino){
        eres[ino][0] = vol*rho*g[0]*0.25;
        eres[ino][1] = vol*rho*g[1]*0.25;
        eres[ino][2] = vol*rho*g[2]*0.25;
      }
      double u0[4][3]; FetchData(&u0[0][0], 4, 3, aIP, aDisp);
      double v0[4][3]; FetchData(&v0[0][0], 4, 3, aIP, aVelo);
      for(int ino=0;ino<4;++ino){
        const double* Mi = aR.data()+aIP[ino]*9;
        for(int idim=0;idim<3;++idim){
          for(int jno=0;jno<4;++jno){
            double Pj1[3]; MatVec3(Pj1, Mi,P[jno]);
            double uj1[3] = {
              P[jno][0]+u0[jno][0]+dt*v0[jno][0]-Pj1[0],
              P[jno][1]+u0[jno][1]+dt*v0[jno][1]-Pj1[1],
              P[jno][2]+u0[jno][2]+dt*v0[jno][2]-Pj1[2] };
            eres[ino][idim] -= emat[ino][jno][idim][0]*uj1[0];
            eres[ino][idim] -= emat[ino][jno][idim][1]*uj1[1];
            eres[ino][idim] -= emat[ino][jno][idim][2]*uj1[2];
          }
        }
      }
    }
    for(int ino=0;ino<4;++ino){
      emat[ino][ino][0][0] += rho*vol*0.25/(dt*dt);
      emat[ino][ino][1][1] += rho*vol*0.25/(dt*dt);
      emat[ino][ino][2][2] += rho*vol*0.25/(dt*dt);
    }
    ////////////////
    for (int ino = 0; ino<4; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0]/dt;
      vec_b[ip*3+1] += eres[ino][1]/dt;
      vec_b[ip*3+2] += eres[ino][2]/dt;
    }
    mat_A.Mearge(4, aIP, 4, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_Stokes3D_Static
(CMatrixSparse<double>& mat_A,
 std::vector<double>& vec_b,
 const double myu,
 const double rho,
 const double g_x,
 const double g_y,
 const double g_z,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTet,
 const std::vector<double>& aVal,
 const std::vector<double>& aVelo)
{
  const int np = (int)aXYZ.size()/3;
  const int nDoF = np*4;
  ////
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  std::vector<int> tmp_buffer(np, -1);
  for (int itet = 0; itet<(int)aTet.size()/4; ++itet){
    const unsigned int i0 = aTet[itet*4+0];
    const unsigned int i1 = aTet[itet*4+1];
    const unsigned int i2 = aTet[itet*4+2];
    const unsigned int i3 = aTet[itet*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; FetchData(&coords[0][0],4,3,aIP, aXYZ.data());
    double velo_press[4][4]; FetchData(&velo_press[0][0],4,4,aIP, aVal.data());
    ////
    double eres[4][4];
    double emat[4][4][4][4];
    ////
    MakeMat_Stokes3D_Static_P1(myu, g_x, g_y, g_z,
                               coords, velo_press,
                               emat, eres);
    for (int ino = 0; ino<4; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*4+0] += eres[ino][0];
      vec_b[ip*4+1] += eres[ino][1];
      vec_b[ip*4+2] += eres[ino][2];
      vec_b[ip*4+3] += eres[ino][3];
    }
    // marge dde
    mat_A.Mearge(4, aIP, 4, aIP,16, &emat[0][0][0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_Stokes3D_Dynamic
(CMatrixSparse<double>& mat_A,
 std::vector<double>& vec_b,
 const double myu,
 const double rho,
 const double g_x,
 const double g_y,
 const double g_z,
 const double dt_timestep,
 const double gamma_newmark,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTet,
 const std::vector<double>& aVal,
 const std::vector<double>& aVelo)
{
  const int np = (int)aXYZ.size()/3;
  const int nDoF = np*4;
  ////
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<(int)aTet.size()/4; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; FetchData(&coords[0][0],4,3,aIP, aXYZ.data());
    double velo_press[4][4]; FetchData(&velo_press[0][0],4,4,aIP, aVal.data());
    double acc_apress[4][4]; FetchData(&acc_apress[0][0],4,4,aIP, aVelo.data());
    ////
    double eres[4][4];
    double emat[4][4][4][4];
    MakeMat_Stokes3D_Dynamic_P1(myu, rho,  g_x, g_y,g_z,
                                dt_timestep, gamma_newmark,
                                coords, velo_press, acc_apress,
                                emat, eres);
    for (int ino = 0; ino<4; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*4+0] += eres[ino][0];
      vec_b[ip*4+1] += eres[ino][1];
      vec_b[ip*4+2] += eres[ino][2];
      vec_b[ip*4+3] += eres[ino][3];
    }
    mat_A.Mearge(4, aIP, 4, aIP,16, &emat[0][0][0][0], tmp_buffer);
  }
}

void dfm2::MergeLinSys_NavierStokes3D_Dynamic
(CMatrixSparse<double>& mat_A,
 std::vector<double>& vec_b,
 const double myu,
 const double rho,
 const double g_x,
 const double g_y,
 const double g_z,
 const double dt_timestep,
 const double gamma_newmark,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTet,
 const std::vector<double>& aVal,
 const std::vector<double>& aVelo)
{
  const int np = (int)aXYZ.size()/3;
  const int nDoF = np*4;
  ////
  mat_A.SetZero();
  vec_b.assign(nDoF, 0.0);
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<(int)aTet.size()/4; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; FetchData(&coords[0][0],4,3,aIP, aXYZ.data());
    double velo_press[4][4]; FetchData(&velo_press[0][0],4,4,aIP, aVal.data());
    double acc_apress[4][4]; FetchData(&acc_apress[0][0],4,4,aIP, aVelo.data());
    ////
    double eres[4][4], emat[4][4][4][4];
    MakeMat_NavierStokes3D_Dynamic_P1(myu, rho,  g_x, g_y,g_z,
                                      dt_timestep, gamma_newmark,
                                      coords, velo_press, acc_apress,
                                      emat, eres);
    for (int ino = 0; ino<4; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*4+0] += eres[ino][0];
      vec_b[ip*4+1] += eres[ino][1];
      vec_b[ip*4+2] += eres[ino][2];
      vec_b[ip*4+3] += eres[ino][3];
    }
    mat_A.Mearge(4, aIP, 4, aIP,16, &emat[0][0][0][0], tmp_buffer);
  }
}

/* ------------------------------------------------------------------------- */



/*
 void Solve_LinearSolid_TetP2()
 {
 
 unsigned int ncorner = aXYZ.size()/3;
 unsigned int nedge = aEdge.size()/2;
 unsigned int npoint = ncorner + nedge;
 unsigned int ndof = npoint*3;
 std::cout << "ndof: " << ndof << std::endl;
 MatrixXd Mat(ndof,ndof);
 //  Eigen::SparseMatrix<double> Mat(ndof,ndof);
 //  Mat.reserve(VectorXi::Constant(ndof,100));
 Mat.setZero();
 
 VectorXd Lhs(ndof);
 Lhs.setZero();
 
 double lambda = 0;
 double myu = 1;
 
 for(unsigned int itet=0;itet<aTet.size()/4;itet++){
 const unsigned int aPo[10] = {
 aTet[itet*4+0],
 aTet[itet*4+1],
 aTet[itet*4+2],
 aTet[itet*4+3],
 aTetEdge[itet*6+0] + ncorner,
 aTetEdge[itet*6+1] + ncorner,
 aTetEdge[itet*6+2] + ncorner,
 aTetEdge[itet*6+3] + ncorner,
 aTetEdge[itet*6+4] + ncorner,
 aTetEdge[itet*6+5] + ncorner };
 unsigned int i0=aPo[0], i1=aPo[1], i2=aPo[2], i3=aPo[3];
 double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
 double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
 double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]};
 double p3[3] = {aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2]};
 ////
 const double vol = TetVolume3D(p0,p1,p2,p3);
 double dldx[4][3];
 double zero_order_term[4];
 TetDlDx(dldx, zero_order_term,   p0,p1,p2,p3);
 ////
 double disp[10][3];
 for(unsigned int ino=0;ino<10;ino++){
 disp[ino][0] = 0;
 disp[ino][1] = 0;
 disp[ino][2] = 0;
 }
 double emat[10][10][3][3], eres[10][3];
 matRes_LinearSolid_TetP2(emat,eres,
 vol,lambda,myu,
 0,1,0, 1,
 dldx,disp);
 
 for(unsigned int ino=0;ino<10;ino++){
 for(unsigned int jno=0;jno<10;jno++){
 unsigned int ipo = aPo[ino];
 unsigned int jpo = aPo[jno];
 for(unsigned int idim=0;idim<3;idim++){
 for(unsigned int jdim=0;jdim<3;jdim++){
 Mat(ipo*3+idim,jpo*3+jdim) += emat[ino][jno][idim][jdim];
 //        Mat.coeffRef(ipo*3+idim, jpo*3+jdim) += emat[ino][jno][idim][jdim];
 }
 }
 }
 }
 for(unsigned int ino=0;ino<10;ino++){
 unsigned int ipo = aPo[ino];
 for(unsigned int idim=0;idim<3;idim++){
 Lhs(ipo*3+idim) += eres[ino][idim];
 }
 }
 //    std::cout << itet << " " << aTet.size()/4 << std::endl;
 }
 
 std::vector<unsigned int> aFix;
 {
 for(unsigned int ipo=0;ipo<aXYZ.size()/3;ipo++){
 double p[3] = {aXYZ[ipo*3+0], aXYZ[ipo*3+1], aXYZ[ipo*3+2]};
 if( p[0] < 0 ){
 aFix.push_back(ipo);
 }
 }
 for(unsigned int iedge=0;iedge<aEdge.size()/2;iedge++){
 unsigned int i0 = aEdge[iedge*2+0];
 unsigned int i1 = aEdge[iedge*2+1];
 double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
 double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
 if( p0[0] < 0 && p1[0] < 0 ){
 aFix.push_back(iedge+ncorner);
 }
 }
 }
 for(unsigned int iifix=0;iifix<aFix.size();iifix++){
 unsigned int ifix0 = aFix[iifix];
 for(int i=0;i<3;i++){
 Lhs(ifix0*3+i) = 0;
 }
 for(unsigned int jpo=0;jpo<npoint;jpo++){
 for(int i=0;i<3;i++){
 for(int j=0;j<3;j++){
 Mat(ifix0*3+i,jpo*3+j) = 0;
 //          Mat.coeffRef(ifix0*3+i,jpo*3+j) = 0;
 }
 }
 }
 for(unsigned int jpo=0;jpo<npoint;jpo++){
 for(int i=0;i<3;i++){
 for(int j=0;j<3;j++){
 Mat(jpo*3+i,ifix0*3+j) = 0;
 //          Mat.coeffRef(jpo*3+i,ifix0*3+j) = 0;
 }
 }
 }
 {
 for(int i=0;i<3;i++){
 Mat(ifix0*3+i,ifix0*3+i) = 1;
 //        Mat.coeffRef(ifix0*3+i,ifix0*3+i) = 1;
 }
 }
 }
 //  Mat.makeCompressed();
 
 VectorXd Rhs;
 Rhs = Mat.partialPivLu().solve(Lhs);
 
 //  std::cout << ndof << std::endl;
 //  ConjugateGradient<MatrixXd> cg;
 //  cg.setMaxIterations(3000);
 //  cg.compute(Mat);
 //  Rhs = cg.solve(Lhs);
 for(unsigned int ip=0;ip<aXYZ.size()/3;ip++){
 for(unsigned int i=0;i<3;i++){
 aDisp[ip*3+i] += Rhs(ip*3+i);
 }
 std::cout << aDisp[ip*3+0] << " " << aDisp[ip*3+1] << " " << aDisp[ip*3+2] << std::endl;
 }
 }
 */

/*
 // this was used
 void Solve_LinearSolid_TetP2_MyMat()
 {
 unsigned int ncorner = aXYZ.size()/3;
 unsigned int nedge = aEdge.size()/2;
 unsigned int npoint = ncorner + nedge;
 
 ////
 MatVec::CMatDia_BlkCrs matFEM;
 MatVec::CVector_Blk vecRes;
 MatVec::CVector_Blk vecUpd;
 MatVec::CBCFlag bc_flag(npoint,3);
 {
 std::vector<unsigned int> aEl;
 aEl.resize(aTet.size()/4*10);
 for(unsigned int itet=0;itet<aTet.size()/4;itet++){
 const unsigned int aPo[10] = {
 aTet[itet*4+0],
 aTet[itet*4+1],
 aTet[itet*4+2],
 aTet[itet*4+3],
 aTetEdge[itet*6+0] + ncorner,
 aTetEdge[itet*6+1] + ncorner,
 aTetEdge[itet*6+2] + ncorner,
 aTetEdge[itet*6+3] + ncorner,
 aTetEdge[itet*6+4] + ncorner,
 aTetEdge[itet*6+5] + ncorner };
 for(unsigned int ino=0;ino<10;ino++){
 aEl[itet*10+ino] = aPo[ino];
 }
 }
 matFEM.Initialize(npoint, 3);
 CIndexedArray crs_edge;
 crs_edge.SetEdge(npoint, aTet.size()/4, 10, aEl);
 matFEM.AddPattern(crs_edge);
 }
 vecRes.Initialize(npoint,3);
 vecUpd.Initialize(npoint,3);
 
 matFEM.SetZero();
 vecRes.SetVectorZero();
 
 unsigned int ndof = npoint*3;
 std::cout << "ndof: " << ndof << std::endl;
 
 for(unsigned int itet=0;itet<aTet.size()/4;itet++){
 const unsigned int aPo[10] = {
 aTet[itet*4+0],
 aTet[itet*4+1],
 aTet[itet*4+2],
 aTet[itet*4+3],
 aTetEdge[itet*6+0] + ncorner,
 aTetEdge[itet*6+1] + ncorner,
 aTetEdge[itet*6+2] + ncorner,
 aTetEdge[itet*6+3] + ncorner,
 aTetEdge[itet*6+4] + ncorner,
 aTetEdge[itet*6+5] + ncorner };
 unsigned int i0=aPo[0], i1=aPo[1], i2=aPo[2], i3=aPo[3];
 double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
 double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
 double p2[3] = {aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2]};
 double p3[3] = {aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2]};
 ////
 const double vol = TetVolume3D(p0,p1,p2,p3);
 double dldx[4][3];
 double zero_order_term[4];
 TetDlDx(dldx, zero_order_term,   p0,p1,p2,p3);
 ////
 double disp[10][3];
 for(unsigned int i=0;i<3;i++){
 disp[0][i] = aDisp[ aPo[0]*3+i];
 disp[1][i] = aDisp[ aPo[1]*3+i];
 disp[2][i] = aDisp[ aPo[2]*3+i];
 disp[3][i] = aDisp[ aPo[3]*3+i];
 disp[4][i] = aDispEdge[ aTetEdge[itet*6+0]*3+i ];
 disp[5][i] = aDispEdge[ aTetEdge[itet*6+1]*3+i ];
 disp[6][i] = aDispEdge[ aTetEdge[itet*6+2]*3+i ];
 disp[7][i] = aDispEdge[ aTetEdge[itet*6+3]*3+i ];
 disp[8][i] = aDispEdge[ aTetEdge[itet*6+4]*3+i ];
 disp[9][i] = aDispEdge[ aTetEdge[itet*6+5]*3+i ];
 }
 ////
 double emat[10][10][3][3], eres[10][3];
 matRes_LinearSolid_TetP2(emat,eres,
 vol,lambda,myu,
 g_x,g_y,g_z, 1,
 dldx,disp);
 
 matFEM.Mearge(10, aPo, 10, aPo, 9, &emat[0][0][0][0]);
 for(unsigned int ino=0;ino<10;ino++){
 unsigned int ipo = aPo[ino];
 for(unsigned int idim=0;idim<3;idim++){
 vecRes.AddValue(ipo,idim,eres[ino][idim]);
 }
 }
 }
 {
 vecRes.AddValue(ino_exforce,0,exforce_vector.x);
 vecRes.AddValue(ino_exforce,1,exforce_vector.y);
 vecRes.AddValue(ino_exforce,2,exforce_vector.z);
 }
 
 {
 for(unsigned int ipo=0;ipo<aXYZ.size()/3;ipo++){
 double p[3] = {aXYZ[ipo*3+0], aXYZ[ipo*3+1], aXYZ[ipo*3+2]};
 double h = p[0]*norm.x + p[1]*norm.y + p[2]*norm.z;
 if( h < -clearance_height ){
 bc_flag.SetBC(ipo, 0);
 bc_flag.SetBC(ipo, 1);
 bc_flag.SetBC(ipo, 2);
 }
 }
 for(unsigned int iedge=0;iedge<aEdge.size()/2;iedge++){
 unsigned int i0 = aEdge[iedge*2+0];
 unsigned int i1 = aEdge[iedge*2+1];
 double p0[3] = {aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2]};
 double p1[3] = {aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2]};
 double h0 = p0[0]*norm.x + p0[1]*norm.y + p0[2]*norm.z;
 double h1 = p1[0]*norm.x + p1[1]*norm.y + p1[2]*norm.z;
 if( p0[0] < -clearance_height && p1[0] < -clearance_height ){
 bc_flag.SetBC(iedge+ncorner,0);
 bc_flag.SetBC(iedge+ncorner,1);
 bc_flag.SetBC(iedge+ncorner,2);
 }
 }
 matFEM.SetBoundaryCondition(bc_flag);
 bc_flag.SetZeroToBCDof(vecRes);
 }
 
 double conv_ratio = 1.0e-5;
 unsigned int iteration = 1000;
 MatVec::Sol::Solve_CG(conv_ratio, iteration, matFEM, vecRes, vecUpd);
 //  MatVec::Sol::Solve_PCG(conv_ratio, iteration, matFEM, vecRes, vecUpd, matPrec,matFEM);
 for(unsigned int ipo=0;ipo<aXYZ.size()/3;ipo++){
 aDisp[ipo*3+0] += vecUpd.GetValue(ipo, 0);
 aDisp[ipo*3+1] += vecUpd.GetValue(ipo, 1);
 aDisp[ipo*3+2] += vecUpd.GetValue(ipo, 2);
 }
 for(unsigned int ipo=0;ipo<aEdge.size()/2;ipo++){
 aDispEdge[ipo*3+0] = vecUpd.GetValue(ipo+ncorner,0);
 aDispEdge[ipo*3+1] = vecUpd.GetValue(ipo+ncorner,1);
 aDispEdge[ipo*3+2] = vecUpd.GetValue(ipo+ncorner,2);
 }
 
 }
 
 */


void dfm2::MergeLinSys_ShellStaticPlateBendingMITC3_MeshTri2D
(CMatrixSparse<double>& mat_A,
 double* vec_b,
 const double thick,
 const double lambda,
 const double myu,
 const double rho,
 const double gravity_z,
 const double* aXY1, int nXY,
 const unsigned int* aTri1, unsigned int nTri,
 const double* aVal)
{
  const int np = nXY;
  std::vector<int> tmp_buffer(np, -1);
  for(unsigned int iel=0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double P[3][2]; FetchData(&P[0][0],  3,2,aIP, aXY1);
    double u[3][3]; FetchData(&u[0][0],   3,3,aIP, aVal);
    ////
    double W=0.0, dW[3][3], ddW[3][3][3][3];
    for(int i=0;i<9;++i){ (&dW[0][0])[i] = 0.0; }
    for(int i=0;i<81;++i){ (&ddW[0][0][0][0])[i] = 0.0; }
    WdWddW_PlateBendingMITC3(W,dW,ddW,
                             P, u,
                             thick, lambda, myu);
    {
      const double A = TriArea2D(P[0],P[1],P[2]);
      dW[0][0] = rho*A*thick/3.0*gravity_z;
      dW[1][0] = rho*A*thick/3.0*gravity_z;
      dW[2][0] = rho*A*thick/3.0*gravity_z;
    }
    for (unsigned int ino = 0; ino<3; ino++){
      const unsigned int ip = aIP[ino];
      vec_b[ip*3+0] += dW[ino][0];
      vec_b[ip*3+1] += dW[ino][1];
      vec_b[ip*3+2] += dW[ino][2];
    }
    // marge dde
    mat_A.Mearge(3, aIP, 3, aIP, 9, &ddW[0][0][0][0], tmp_buffer);
  }
}

void dfm2::MassLumped_ShellPlateBendingMITC3
(double* aM,
 double rho, double thick,
 const double* aXY, unsigned int nXY,
 const unsigned int* aTri, unsigned int nTri)
{
  const unsigned int nDoF = nXY*3;
  for(unsigned int i=0;i<nDoF;++i){ aM[i] = 0.0; }
  for(unsigned int it=0;it<nTri;++it){
    const unsigned int i0 = aTri[it*3+0]; assert(i0<nXY);
    const unsigned int i1 = aTri[it*3+1]; assert(i1<nXY);
    const unsigned int i2 = aTri[it*3+2]; assert(i2<nXY);
    const double* p0 = aXY+i0*2;
    const double* p1 = aXY+i1*2;
    const double* p2 = aXY+i2*2;
    const double a012 = TriArea2D(p0, p1, p2);
    double m0 = a012/3.0*rho*thick;
    double m1 = a012/3.0*rho*thick*thick*thick/12.0;
    double m2 = m1;
    aM[i0*3+0] += m0;  aM[i1*3+0] += m0;  aM[i2*3+0] += m0;
    aM[i0*3+1] += m1;  aM[i1*3+1] += m1;  aM[i2*3+1] += m1;
    aM[i0*3+2] += m2;  aM[i1*3+2] += m2;  aM[i2*3+2] += m2;
  }
}
