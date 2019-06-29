/**
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <stdio.h>

#include "delfem2/fem.h"
#include "delfem2/fem_ematrix.h"
#include "delfem2/matrix_sparse.h"
#include "delfem2/ilu_sparse.h"


static void FetchData
(double* val_to,
 int nno, int ndim,
 const unsigned int* aIP,
 const double* val_from,
 int nstride, int noffset)
{
  assert( nstride >= ndim );
  for(int ino=0;ino<nno;++ino){
    int ip = aIP[ino];
    for(int idim=0;idim<ndim;++idim){
      val_to[ino*ndim+idim] = val_from[ip*nstride+noffset+idim];
    }
  }
}

void XPlusAY
(std::vector<double>& X,
 const int nDoF,
 const std::vector<int>& aBCFlag,
 double alpha,
 const std::vector<double>& Y)
{
  for(int i=0;i<nDoF;++i ){
    if( aBCFlag[i] !=0 ) continue;
    X[i] += alpha*Y[i];
  }
}

void XPlusAYBZ
(std::vector<double>& X,
 const int nDoF,
 const std::vector<int>& aBCFlag,
 double alpha,
 const std::vector<double>& Y,
 double beta,
 const std::vector<double>& Z)
{
  for(int i=0;i<nDoF;++i ){
    if( aBCFlag[i] !=0 ) continue;
    X[i] += alpha*Y[i] + beta*Z[i];
  }
}

void XPlusAYBZCW
(std::vector<double>& X,
 const int nDoF,
 const std::vector<int>& aBCFlag,
 double alpha,
 const std::vector<double>& Y,
 double beta,
 const std::vector<double>& Z,
 double gamma,
 const std::vector<double>& W)
{
  for(int i=0;i<nDoF;++i ){
    if( aBCFlag[i] !=0 ) continue;
    X[i] += alpha*Y[i] + beta*Z[i] + gamma*W[i];
  }
}

// set boundary condition
void setRHS_Zero
(std::vector<double>& vec_b,
 const std::vector<int>& aBCFlag,
 int iflag_nonzero)
{
  const int ndof = (int)vec_b.size();
  for (int i=0;i<ndof;++i){
    if (aBCFlag[i]==iflag_nonzero) continue;
    vec_b[i] = 0;
  }
}

void setRHS_MasterSlave
(double* vec_b,
 int nDoF,
 const int* aMSFlag)
{
  for(int idof=0;idof<nDoF;++idof){
    int jdof = aMSFlag[idof];
    if( jdof == -1 ) continue;
    vec_b[jdof] += vec_b[idof];
    vec_b[idof] = 0;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




void SolveLinSys_PCG
(const CMatrixSparse& mat_A,
 std::vector<double>& vec_b,
 std::vector<double>& vec_x,
 CPreconditionerILU& ilu_A,
 double& conv_ratio,
 int& iteration)
{
  // set ILU preconditioner
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  // solve linear system
  //Solve_CG(conv_ratio, iteration, mat_A, vec_b, vec_x);
  //  Solve_BiCGSTAB(conv_ratio, iteration, mat_A, vec_b, vec_x);
  //  Solve_PBiCGSTAB(conv_ratio, iteration, mat_A, ilu_A, vec_b, vec_x);
  vec_x.resize(vec_b.size());
  Solve_PCG(vec_b.data(), vec_x.data(), conv_ratio, iteration, mat_A, ilu_A);
  std::cout<<"  conv_ratio:"<<conv_ratio<<"  iteration:"<<iteration<<std::endl;
}

bool SolveLinSys_BiCGStab
(CMatrixSparse& mat_A,
 std::vector<double>& vec_b,
 std::vector<double>& vec_x,
 CPreconditionerILU& ilu_A,
 double& conv_ratio,
 int& iteration)
{
  // set ILU preconditioner
  ilu_A.SetValueILU(mat_A);
  bool res_ilu = ilu_A.DoILUDecomp();
  if( !res_ilu ){ return false; }
  // solve linear system
  //  double conv_ratio = 1.0e-4;
  //  int iteration = 1000;
  //  Solve_CG(conv_ratio, iteration, mat_A, vec_b, vec_x);
  //  Solve_BiCGSTAB(conv_ratio, iteration, mat_A, vec_b, vec_x);
  vec_x.resize(vec_b.size());
  Solve_PBiCGStab(vec_b.data(), vec_x.data(), conv_ratio, iteration, mat_A, ilu_A);
  /// Solve_PCG(conv_ratio, iteration, mat_A, ilu_A, vec_b, vec_x);
  //  std::cout<<"  interative solver --- conv_ratio:"<<conv_ratio<<"  iteration:"<<iteration<<std::endl;
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void MergeLinSys_Poission_MeshTri2D
(CMatrixSparse& mat_A,
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
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1,2,0);
    const double value[3] = { aVal[i0], aVal[i1], aVal[i2] };
    ////
    double eres[3];
    double emat[3][3];
    MakeMat_Poisson2D_P1
    (alpha, source,
     coords, value,
     eres,emat);
    for (int ino = 0; ino<3; ino++){
      const int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    mat_A.Mearge(3, aIP, 3, aIP, 1, &emat[0][0], tmp_buffer);
  }
}

void MergeLinSys_Poission_MeshTet3D
(CMatrixSparse& mat_A,
 double* vec_b,
 const double alpha,
 const double source,
 const double* aXYZ, int nXYZ,
 const unsigned int* aTet, int nTet,
 const double* aVal)
{
  const int np = nXYZ;
  /////
  std::vector<int> tmp_buffer(np, -1);
  for (int itet = 0; itet<nTet; ++itet){
    const unsigned int i0 = aTet[itet*4+0];
    const unsigned int i1 = aTet[itet*4+1];
    const unsigned int i2 = aTet[itet*4+2];
    const unsigned int i3 = aTet[itet*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; FetchData(&coords[0][0],4,3,aIP, aXYZ,3,0);
    const double value[4] = { aVal[i0], aVal[i1], aVal[i2], aVal[i3] };
    ////
    double eres[4], emat[4][4];
    MakeMat_Poisson3D_P1
    (alpha, source,
     coords, value,
     eres,emat);
    for (int ino = 0; ino<4; ino++){
      const int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    // marge dde
    mat_A.Mearge(4, aIP, 4, aIP, 1, &emat[0][0], tmp_buffer);
  }
}

void MergeLinSys_Diffusion_MeshTri2D
(CMatrixSparse& mat_A,
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
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1, 2,0);
    const double value[3] = { aVal[ i0], aVal[ i1], aVal[ i2] };
    const double velo[ 3] = { aVelo[i0], aVelo[i1], aVelo[i2] };
    ////
    double eres[3];
    double emat[3][3];
    MakeMat_Diffusion2D_P1
    (alpha, source,
     dt_timestep, gamma_newmark, rho,
     coords, value, velo,
     eres,emat);
    for (int ino = 0; ino<3; ino++){
      const int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    mat_A.Mearge(3, aIP, 3, aIP, 1, &emat[0][0], tmp_buffer);
  }
}

void MergeLinSys_Diffusion_MeshTet3D
(CMatrixSparse& mat_A,
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
//  const int nDoF = np;
  ////
//  mat_A.SetZero();
//  vec_b.assign(nDoF, 0.0);
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTet; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; FetchData(&coords[0][0],4,3,aIP, aXYZ,3,0);
    const double value[4] = { aVal[ i0], aVal[ i1], aVal[ i2], aVal[ i3] };
    const double velo[ 4] = { aVelo[i0], aVelo[i1], aVelo[i2], aVelo[i3] };
    ////
    double eres[4];
    double emat[4][4];
    MakeMat_Diffusion3D_P1
    (alpha, source,
     dt_timestep, gamma_newmark, rho,
     coords, value, velo,
     eres,emat);
    for (int ino = 0; ino<4; ino++){
      const int ip = aIP[ino];
      vec_b[ip] += eres[ino];
    }
    mat_A.Mearge(4, aIP, 4, aIP, 1, &emat[0][0], tmp_buffer);
  }
}

void MergeLinSys_SolidStaticLinear_MeshTri2D
(CMatrixSparse& mat_A,
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
//  const int nDoF = np*2;
  ////
//  mat_A.SetZero();
//  for(int i=0;i<nDoF;++i){ vec_b[i] = 0.0; }
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTri; ++iel){
    const unsigned int i0 = aTri1[iel*3+0];
    const unsigned int i1 = aTri1[iel*3+1];
    const unsigned int i2 = aTri1[iel*3+2];
    const unsigned int aIP[3] = {i0,i1,i2};
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1,2,0);
    double disps[3][2]; FetchData(&disps[0][0],3,2,aIP, aVal,2,0);
    ////
    double eres[3][2];
    double emat[3][3][2][2];
    MakeMat_LinearSolid2D_Static_P1
    (myu, lambda,
     rho, g_x, g_y,
     disps, coords,
     eres,emat);
    for (int ino = 0; ino<3; ino++){
      const int ip = aIP[ino];
      vec_b[ip*2+0] += eres[ino][0];
      vec_b[ip*2+1] += eres[ino][1];
    }
    // marge dde
    mat_A.Mearge(3, aIP, 3, aIP, 4, &emat[0][0][0][0], tmp_buffer);
  }
}

void MergeLinSys_SolidDynamicLinear_MeshTri2D
(CMatrixSparse& mat_A,
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
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1,2,0);
    double disps[3][2]; FetchData(&disps[0][0],3,2,aIP, aVal,2,0);
    double velos[3][2]; FetchData(&velos[0][0],3,2,aIP, aVelo,2,0);
    double accs[3][2]; FetchData(&accs[0][0],3,2,aIP, aAcc,2,0);
    ////
    double eres[3][2];
    double emat[3][3][2][2];
    MakeMat_LinearSolid2D_Dynamic_P1
    (myu, lambda,
     rho, g_x, g_y,
     dt_timestep, gamma_newmark, beta_newmark,
     disps, velos, accs, coords,
     eres,emat,
     true);
    for (int ino = 0; ino<3; ino++){
      const int ip = aIP[ino];
      vec_b[ip*2+0] += eres[ino][0];
      vec_b[ip*2+1] += eres[ino][1];
    }
    // marge dde
    mat_A.Mearge(3, aIP, 3, aIP, 4, &emat[0][0][0][0], tmp_buffer);
  }
}

void MergeLinSys_StokesStatic2D
(CMatrixSparse& mat_A,
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
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1,2,0);
    double velo_press[3][3]; FetchData(&velo_press[0][0],3,3,aIP, aVal,3,0);
    ////
    double eres[3][3];
    double emat[3][3][3][3];
    ////
    MakeMat_Stokes2D_Static_P1(myu, g_x, g_y, coords, velo_press, emat, eres);
    for (int ino = 0; ino<3; ino++){
      const int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0];
      vec_b[ip*3+1] += eres[ino][1];
      vec_b[ip*3+2] += eres[ino][2];
    }
    // marge dde
    mat_A.Mearge(3, aIP, 3, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
}

void MergeLinSys_StokesDynamic2D
(CMatrixSparse& mat_A,
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
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1,2,0);
    double velo_press[3][3]; FetchData(&velo_press[0][0],3,3,aIP, aVal, 3,0);
    double acc_apress[3][3]; FetchData(&acc_apress[0][0],3,3,aIP, aVelo, 3,0);
    ////
    double eres[3][3];
    double emat[3][3][3][3];
    ////
    MakeMat_Stokes2D_Dynamic_P1(myu, rho,  g_x, g_y,
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

void MergeLinSys_NavierStokes2D
(CMatrixSparse& mat_A,
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
    double coords[3][2]; FetchData(&coords[0][0],3,2,aIP, aXY1,2,0);
    double velo_press[3][3]; FetchData(&velo_press[0][0],3,3,aIP, aVal,3,0);
    double acc_apress[3][3]; FetchData(&acc_apress[0][0],3,3,aIP, aVelo, 3,0);
    ////
    double eres[3][3], emat[3][3][3][3];
    MakeMat_NavierStokes2D_Dynamic_P1(myu, rho,  g_x, g_y,
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


// compute total energy and its first and second derivatives
double MergeLinSys_Cloth
(CMatrixSparse& ddW, // (out) second derivative of energy
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
      for(int i=0;i<3;i++){ c[ino][i] = aXYZ [ip*3+i]; }
    }
    double e, de[3][3], dde[3][3][3][3];
    WdWddW_CST( e,de,dde, C,c, lambda,myu );
    W += e;  // marge energy
    // marge de
    for(int ino=0;ino<3;ino++){
      const int ip = aIP[ino];
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
      const int ip = aIP[ino];
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




double MergeLinSys_Contact
(CMatrixSparse& ddW,
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


void MergeLinSys_SolidStaticLinear_MeshTet3D
(CMatrixSparse& mat_A,
 double* vec_b,
 const double myu,
 const double lambda,
 const double rho,
 const double g_x,
 const double g_y,
 const double g_z,
 const double* aXYZ, int nXYZ,
 const unsigned int* aTet, int nTet,
 const double* aVal)
{
  const int np = nXYZ;
//  const int nDoF = np*3;
  //////
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTet; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = { i0, i1, i2, i3 };
    double coords[4][3]; FetchData(&coords[0][0], 4, 3, aIP, aXYZ, 3, 0);
    double disps[4][3]; FetchData(&disps[0][0], 4, 3, aIP, aVal, 3, 0);
    ////
    double eres[4][3];
    double emat[4][4][3][3];
    MakeMat_LinearSolid3D_Static_P1(myu, lambda,
                                    rho, g_x, g_y, g_z,
                                    coords, disps,
                                    emat,eres);
    for (int ino = 0; ino<4; ino++){
      const int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0];
      vec_b[ip*3+1] += eres[ino][1];
      vec_b[ip*3+2] += eres[ino][2];
    }
    // marge dde
    mat_A.Mearge(4, aIP, 4, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
}

void MergeLinSys_LinearSolid3D_Static_Q1
(CMatrixSparse& mat_A,
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
    double coords[8][3]; FetchData(&coords[0][0], 8, 3, aIP, aXYZ.data(), 3, 0);
    double disps[8][3]; FetchData(&disps[0][0], 8, 3, aIP, aVal.data(), 3, 0);
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

void MergeLinSys_SolidDynamicLinear_MeshTet3D
(CMatrixSparse& mat_A,
 double* vec_b,
 const double myu,
 const double lambda,
 const double rho,
 const double g_x,
 const double g_y,
 const double g_z,
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
//  const int nDoF = np*3;
  ////
  std::vector<int> tmp_buffer(np, -1);
  for (int iel = 0; iel<nTet; ++iel){
    const unsigned int i0 = aTet[iel*4+0];
    const unsigned int i1 = aTet[iel*4+1];
    const unsigned int i2 = aTet[iel*4+2];
    const unsigned int i3 = aTet[iel*4+3];
    const unsigned int aIP[4] = {i0,i1,i2,i3};
    double coords[4][3]; FetchData(&coords[0][0],4,3,aIP, aXYZ, 3,0);
    double disps[4][3];  FetchData(&disps[0][0], 4,3,aIP, aVal, 3,0);
    double velos[4][3];  FetchData(&velos[0][0], 4,3,aIP, aVelo, 3,0);
    double accs[4][3];   FetchData(&accs[0][0],  4,3,aIP, aAcc, 3,0);
    ////
    double eres[4][3], emat[4][4][3][3];
    MakeMat_LinearSolid3D_Dynamic_P1
    (myu, lambda,
     rho, g_x, g_y, g_z,
     dt_timestep, gamma_newmark, beta_newmark,
     disps, velos, accs, coords,
     eres,emat,
     true);
    for (int ino = 0; ino<4; ino++){
      const int ip = aIP[ino];
      vec_b[ip*3+0] += eres[ino][0];
      vec_b[ip*3+1] += eres[ino][1];
      vec_b[ip*3+2] += eres[ino][2];
    }
    // marge dde
    mat_A.Mearge(4, aIP, 4, aIP, 9, &emat[0][0][0][0], tmp_buffer);
  }
}



void MergeLinSys_Stokes3D_Static
(CMatrixSparse& mat_A,
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
    double coords[4][3]; FetchData(&coords[0][0],4,3,aIP, aXYZ.data(),3,0);
    double velo_press[4][4]; FetchData(&velo_press[0][0],4,4,aIP, aVal.data(),4,0);
    ////
    double eres[4][4];
    double emat[4][4][4][4];
    ////
    MakeMat_Stokes3D_Static_P1(myu, g_x, g_y, g_z,
                               coords, velo_press,
                               emat, eres);
    for (int ino = 0; ino<4; ino++){
      const int ip = aIP[ino];
      vec_b[ip*4+0] += eres[ino][0];
      vec_b[ip*4+1] += eres[ino][1];
      vec_b[ip*4+2] += eres[ino][2];
      vec_b[ip*4+3] += eres[ino][3];
    }
    // marge dde
    mat_A.Mearge(4, aIP, 4, aIP,16, &emat[0][0][0][0], tmp_buffer);
  }
}

void MergeLinSys_Stokes3D_Dynamic
(CMatrixSparse& mat_A,
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
    double coords[4][3]; FetchData(&coords[0][0],4,3,aIP, aXYZ.data(),3,0);
    double velo_press[4][4]; FetchData(&velo_press[0][0],4,4,aIP, aVal.data(), 4,0);
    double acc_apress[4][4]; FetchData(&acc_apress[0][0],4,4,aIP, aVelo.data(),4,0);
    ////
    double eres[4][4];
    double emat[4][4][4][4];
    MakeMat_Stokes3D_Dynamic_P1(myu, rho,  g_x, g_y,g_z,
                                dt_timestep, gamma_newmark,
                                coords, velo_press, acc_apress,
                                emat, eres);
    for (int ino = 0; ino<4; ino++){
      const int ip = aIP[ino];
      vec_b[ip*4+0] += eres[ino][0];
      vec_b[ip*4+1] += eres[ino][1];
      vec_b[ip*4+2] += eres[ino][2];
      vec_b[ip*4+3] += eres[ino][3];
    }
    mat_A.Mearge(4, aIP, 4, aIP,16, &emat[0][0][0][0], tmp_buffer);
  }
}

void MergeLinSys_NavierStokes3D_Dynamic
(CMatrixSparse& mat_A,
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
    double coords[4][3]; FetchData(&coords[0][0],4,3,aIP, aXYZ.data(),3,0);
    double velo_press[4][4]; FetchData(&velo_press[0][0],4,4,aIP, aVal.data(), 4,0);
    double acc_apress[4][4]; FetchData(&acc_apress[0][0],4,4,aIP, aVelo.data(),4,0);
    ////
    double eres[4][4], emat[4][4][4][4];
    MakeMat_NavierStokes3D_Dynamic_P1(myu, rho,  g_x, g_y,g_z,
                                      dt_timestep, gamma_newmark,
                                      coords, velo_press, acc_apress,
                                      emat, eres);
    for (int ino = 0; ino<4; ino++){
      const int ip = aIP[ino];
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
