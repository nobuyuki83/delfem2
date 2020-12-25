/*
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

#include <complex>
#include "delfem2/femem2.h"
#include "delfem2/femem3.h"
#include "delfem2/lsmats.h"
//
#include "delfem2/fem_emats.h"

namespace delfem2 {
namespace fem_emats {

DFM2_INLINE void FetchData(
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


// area of a triangle
DFM2_INLINE double TriArea2D(const double p0[], const double p1[], const double p2[]){
  return 0.5*((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]));
}

DFM2_INLINE double TetVolume3D
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
DFM2_INLINE void TetDlDx(
    double dldx[][3], double a[],
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


DFM2_INLINE void MatVec3
(double y[3],
 const double m[9], const double x[3]){
  y[0] = m[0]*x[0] + m[1]*x[1] + m[2]*x[2];
  y[1] = m[3]*x[0] + m[4]*x[1] + m[5]*x[2];
  y[2] = m[6]*x[0] + m[7]*x[1] + m[8]*x[2];
}

DFM2_INLINE void MatMat3
(double* C,
 const double* A, const double* B)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      C[i*3+j] = A[i*3+0]*B[0*3+j] + A[i*3+1]*B[1*3+j] + A[i*3+2]*B[2*3+j];
    }
  }
}

DFM2_INLINE void MatMatTrans3
(double* C,
 const double* A, const double* B)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      C[i*3+j] = A[i*3+0]*B[j*3+0] + A[i*3+1]*B[j*3+1] + A[i*3+2]*B[j*3+2];
    }
  }
}

}
}

// ==================================================

void delfem2::MergeLinSys_Poission_MeshTri2D(
    CMatrixSparse<double>& mat_A,
    double* vec_b,
    const double alpha,
    const double source,
    const double* aXY1,
    unsigned int np,
    const unsigned int* aTri1,
    unsigned int nTri,
    const double* aVal)
{
  namespace lcl = ::delfem2::fem_emats;
  const unsigned int nDoF = np;
  //
  std::vector<int> tmp_buffer(nDoF, -1);
  for (unsigned int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; lcl::FetchData(&coords[0][0],3,2,aIP, aXY1);
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

void delfem2::MergeLinSys_Poission_MeshTet3D(
    CMatrixSparse<double>& mat_A,
    double* vec_b,
    const double alpha,
    const double source,
    const double* aXYZ,
    unsigned int nXYZ,
    const unsigned int* aTet,
    unsigned int nTet,
    const double* aVal)
{
  namespace lcl = ::delfem2::fem_emats;
  const unsigned int np = nXYZ;
  std::vector<int> tmp_buffer(np, -1);
  for (unsigned int itet = 0; itet<nTet; ++itet){
    const unsigned int i0 = aTet[itet*4+0];
    const unsigned int i1 = aTet[itet*4+1];
    const unsigned int i2 = aTet[itet*4+2];
    const unsigned int i3 = aTet[itet*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; lcl::FetchData(&coords[0][0],4,3,aIP, aXYZ);
    const double value[4] = { aVal[i0], aVal[i1], aVal[i2], aVal[i3] };
    //
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

void delfem2::MergeLinSys_Diffusion_MeshTri2D(
    CMatrixSparse<double>& mat_A,
    double* vec_b,
    const double alpha,
    const double rho,
    const double source,
    const double dt_timestep,
    const double gamma_newmark,
    const double* aXY1,
    unsigned int nXY,
    const unsigned int* aTri1,
    unsigned int nTri,
    const double* aVal,
    const double* aVelo)
{
//  const int nDoF = nXY;
  ////
//  mat_A.SetZero();
//  for(int idof=0;idof<nDoF;++idof){ vec_b[idof] = 0.0; }
  namespace lcl = ::delfem2::fem_emats;
  std::vector<int> tmp_buffer(nXY, -1);
  for (unsigned int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; lcl::FetchData(&coords[0][0],3,2,aIP, aXY1);
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

void delfem2::MergeLinSys_Diffusion_MeshTet3D(
    CMatrixSparse<double>& mat_A,
    double* vec_b,
    const double alpha,
    const double rho,
    const double source,
    const double dt_timestep,
    const double gamma_newmark,
    const double* aXYZ,
    unsigned int nXYZ,
    const unsigned int* aTet,
    unsigned int nTet,
    const double* aVal,
    const double* aVelo)
{
  namespace lcl = ::delfem2::fem_emats;
  const int np = nXYZ;
  std::vector<int> tmp_buffer(np, -1);
  for (unsigned int iel = 0; iel<nTet; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; lcl::FetchData(&coords[0][0],4,3,aIP, aXYZ);
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

// above 2D
// -------------------------------------------------
// below 3D

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

