/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/v3q.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/lsilu_mats.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/lsmats.h"
#include "delfem2/lsvecx.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/femnavierstokes.h"
#include "delfem2/femstokes.h"
#include "delfem2/fempoisson.h"
#include "delfem2/femsolidlinear.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/iss.h"
#include "delfem2/jagarray.h"
#include <GLFW/glfw3.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>


namespace dfm2 = delfem2;

// ------------------------

std::vector<double> aVal;
std::vector<double> aVelo;
std::vector<double> aAcc;

double cur_time = 0.0;
double dt_timestep = 0.01;
//double gamma_newmark = 0.6;
//double beta_newmark = 0.36;
double gamma_newmark = 1;
double beta_newmark = 0.5;

// ----------------------------

void InitializeProblem_Poisson(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  const unsigned int np = aXYZ.size()/3;
  aVal.assign(np, 0.0);
  aBCFlag.assign(np, 0);
  for(int ip=0;ip<np;++ip){
//    const double px = aXYZ[ip*3+0];
    const double py = aXYZ[ip*3+1];
//    const double pz = aXYZ[ip*3+2];
    if( py > +0.4 ){
      aBCFlag[ip] = 1;
      aVal[ip] = 1.0;
    }
    if( py < -0.4 ){
      aBCFlag[ip] = 1;
      aVal[ip] = 0.0;
    }
  }
  //
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                             aTet.data(), aTet.size()/4, 4,
                             (int)aXYZ.size()/3);
  dfm2::JArray_Sort(psup_ind, psup);
  ////
  mat_A.Initialize(np, 1, true); // diagonal
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  ilu_A.Initialize_ILUk(mat_A, 0);
}

void SolveProblem_Poisson(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    std::vector<double>& vec_b,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np;
  // -----------------------
  const double alpha = 1.0;
  const double source = 0.0;
  mat_A.setZero();
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_Poission_MeshTet3D(
      mat_A,vec_b.data(),
      alpha,source,
      aXYZ.data(), aXYZ.size()/3,
      aTet.data(), aTet.size()/4,
      aVal.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag,0);
  // ------------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  vec_x.resize(vec_b.size());
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    dfm2::Solve_PCG(
        dfm2::CVecXd(vec_b),
        dfm2::CVecXd(vec_x),
        dfm2::CVecXd(tmp0),
        dfm2::CVecXd(tmp1),
        conv_ratio, iteration, mat_A, ilu_A);
  }
  // ------------------------------
  dfm2::XPlusAY(aVal,nDoF,aBCFlag,
                1.0,vec_x);
}


void InitializeProblem_Diffusion(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  //double len = 1.1;
  const int np = (int)aXYZ.size()/3;
  aVal.assign(np, 0.0);
  aVelo.assign(np, 0.0);
  aBCFlag.assign(np, 0);
  for(int ip=0;ip<np;++ip){
//    const double px = aXYZ[ip*3+0];
    const double py = aXYZ[ip*3+1];
//    const double pz = aXYZ[ip*3+2];
    if( py > +0.4 ){
      aBCFlag[ip] = 1;
      aVal[ip] = 1.0;
    }
    if( py < -0.4 ){
      aBCFlag[ip] = 1;
      aVal[ip] = 0.0;
    }
  }
  //
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                             aTet.data(), aTet.size()/4, 4,
                             (int)aXYZ.size()/3);
  dfm2::JArray_Sort(psup_ind, psup);
  ////
  mat_A.Initialize(np, 1, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  ilu_A.Initialize_ILU0(mat_A);
}

void SolveProblem_Diffusion(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    std::vector<double>& vec_b,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np;
  // -----------------
  const double alpha = 1.0;
  const double rho = 1.0;
  const double source = 1.0;
  mat_A.setZero();
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_Diffusion_MeshTet3D(
      mat_A,vec_b.data(),
      alpha, rho, source,
      dt_timestep, gamma_newmark,
      aXYZ.data(), aXYZ.size()/3,
      aTet.data(), aTet.size()/4,
      aVal.data(),aVelo.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag,0);
  // ------------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  vec_x.resize(vec_b.size());
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    Solve_PCG(dfm2::CVecXd(vec_b),dfm2::CVecXd(vec_x),dfm2::CVecXd(tmp0),dfm2::CVecXd(tmp1),
        conv_ratio, iteration,
        mat_A, ilu_A);
  }
  // ----------------------
  dfm2::XPlusAYBZ(aVal,nDoF,aBCFlag,
                  dt_timestep*gamma_newmark,vec_x,
                  dt_timestep,aVelo);
  dfm2::XPlusAY(aVelo,nDoF,aBCFlag,
                1.0,vec_x);
  // -----------------------
  dt_timestep = 0.03;
}


void InitializeProblem_ShellEigenPB(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  // set boundary condition
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np*3;
  aVal.assign(nDoF, 0.0);
  aBCFlag.assign(nDoF, 0);
  for (int ip = 0; ip<np; ++ip){
    const double px = aXYZ[ip*3+0];
    const double py = aXYZ[ip*3+1];
    const double pz = aXYZ[ip*3+2];
    if (py > 0.45){
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aBCFlag[ip*3+2] = 1;
    }
  }
  //
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                             aTet.data(), aTet.size()/4, 4,
                             (int)aXYZ.size()/3);
  dfm2::JArray_Sort(psup_ind, psup);
  //
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(
      psup_ind.data(), psup_ind.size(),
      psup.data(),     psup.size());
  ilu_A.Initialize_ILU0(mat_A);
}


void SolveProblem_LinearSolid_Static(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    std::vector<double>& vec_b,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np*3;
  //
  double myu = 1.0;
  double lambda = 1.0;
  double rho = 1.0;
  double g[3] = {0.0, -0.5, 0.0};
  mat_A.setZero();
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_SolidLinear_Static_MeshTet3D(
      mat_A, vec_b.data(),
      myu, lambda, rho, g,
      aXYZ.data(), aXYZ.size()/3,
      aTet.data(), aTet.size()/4,
      aVal.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag,0);
  // --------------------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  vec_x.resize(vec_b.size());
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    Solve_PCG(dfm2::CVecXd(vec_b),dfm2::CVecXd(vec_x),dfm2::CVecXd(tmp0),dfm2::CVecXd(tmp1),
        conv_ratio, iteration,
        mat_A, ilu_A);
  }
  // ------------------------------
  dfm2::XPlusAY(aVal, nDoF, aBCFlag,
    1.0, vec_x);
}


void InitializeProblem_LinearSolid_Dynamic(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  const double len = 1.1;
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np*3;
  ////
  aVal.assign(nDoF, 0.0);
  aVelo.assign(nDoF, 0.0);
  aAcc.assign(nDoF, 0.0);
  aBCFlag.assign(nDoF, 0);
  for(int ip=0;ip<np;++ip){
    const double px = aXYZ[ip*3+0];
    const double py = aXYZ[ip*3+1];
    const double pz = aXYZ[ip*3+2];
    if (py > 0.45){
      aBCFlag[ip*3+0] = 1;
      aBCFlag[ip*3+1] = 1;
      aBCFlag[ip*3+2] = 1;
    }
  }
  //
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(
      psup_ind, psup,
      aTet.data(), aTet.size()/4, 4,
      (int)aXYZ.size()/3);
  dfm2::JArray_Sort(psup_ind, psup);
  mat_A.Initialize(np, 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  ilu_A.Initialize_ILU0(mat_A);
  dt_timestep = 0.03;
}

void SolveProblem_LinearSolid_Dynamic(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    std::vector<double>& vec_b,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np*3;
  // --------
  double myu = 10.0;
  double lambda = 1.0;
  double rho = 1.0;
  const double g[3] = {0.0, -3, 0.0};
  mat_A.setZero();
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_SolidLinear_NewmarkBeta_MeshTet3D(
      mat_A,vec_b.data(),
      myu,lambda,rho,g,
      dt_timestep,gamma_newmark,beta_newmark,
      aXYZ.data(), aXYZ.size()/3,
      aTet.data(), aTet.size()/4,
      aVal.data(),aVelo.data(),aAcc.data());
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag,0);
  // ----------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  vec_x.resize(vec_b.size());
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    dfm2::Solve_PCG(
        dfm2::CVecXd(vec_b),dfm2::CVecXd(vec_x),dfm2::CVecXd(tmp0),dfm2::CVecXd(tmp1),
        conv_ratio, iteration,
        mat_A, ilu_A);
  }
  // -----------------------
  dfm2::XPlusAYBZCW(aVal, nDoF, aBCFlag,
                    dt_timestep,aVelo,
                    0.5*dt_timestep*dt_timestep,aAcc,
                    dt_timestep*dt_timestep*beta_newmark,vec_x);
  dfm2::XPlusAYBZ(aVelo,nDoF, aBCFlag,
                  dt_timestep*gamma_newmark,vec_x,
                  dt_timestep,aAcc);
  dfm2::XPlusAY(aAcc, nDoF, aBCFlag,
                1.0, vec_x);
}


void InitializeProblem_Stokes_Static(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  // set boundary condition
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np*4;
  //
  aVal.assign(np*4, 0.0);
  aBCFlag.assign(nDoF, 0);
  assert(aIsSurf.size() == np);
  for(unsigned int ip=0;ip<np;++ip){
//    const double px = aXYZ[ip*3+0];
    const double py = aXYZ[ip*3+1];
//    const double pz = aXYZ[ip*3+2];
    if( aIsSurf[ip] == 1 ){
      aBCFlag[ip*4+0] = 1;
      aBCFlag[ip*4+1] = 1;
      aBCFlag[ip*4+2] = 1;
    }
    if( py > +0.45 ){
      aVal[ip*4+0] = 10;
    }
    /*
    if( py > +0.45 ){
      aBCFlag[ip*4+0] = 1;
      aBCFlag[ip*4+1] = 1;
      aBCFlag[ip*4+2] = 1;
      aVal[ip*4+1] = 10;
    }
     */
  }
  /*
  CJaggedArray crs;
  crs.SetEdgeOfElem(aTet, (int)aTet.size()/4, 4, (int)aXYZ.size()/3, false);
  crs.Sort();
   */
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                             aTet.data(), aTet.size()/4, 4,
                             (int)aXYZ.size()/3);
  dfm2::JArray_Sort(psup_ind, psup);
  //
  mat_A.Initialize(np, 4, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  ilu_A.Initialize_ILU0(mat_A);
  //
}


void SolveProblem_Stokes_Static(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    std::vector<double>& vec_b,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np*4;
  // -------------
  double myu = 1.0;
  double rho = 1.0;
  double g_x = 0.0;
  double g_y = -0.0;
  double g_z = -0.0;
  mat_A.setZero();
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_Stokes3D_Static(
      mat_A,vec_b,
      myu,rho,g_x,g_y,g_z,
      aXYZ,aTet,
      aVal,aVelo);
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag,0);
  // ----------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  vec_x.resize(vec_b.size());
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    Solve_PCG(
        dfm2::CVecXd(vec_b),
        dfm2::CVecXd(vec_x),
        dfm2::CVecXd(tmp0),
        dfm2::CVecXd(tmp1),
        conv_ratio, iteration,
        mat_A, ilu_A);
  }
  // ----------------------
  dfm2::XPlusAY(aVal, nDoF, aBCFlag, 1.0, vec_x);
}

void InitializeProblem_Stokes_Dynamic(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  // set boundary condition
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np*4;
  aVal.assign(nDoF, 0.0);
  aVelo.assign(nDoF, 0.0);
  aBCFlag.assign(nDoF, 0);
  assert(aIsSurf.size() == np);
  for(unsigned int ip=0;ip<np;++ip){
//    const double px = aXYZ[ip*3+0];
    const double py = aXYZ[ip*3+1];
//    const double pz = aXYZ[ip*3+2];
    if( aIsSurf[ip] == 1 ){
      aBCFlag[ip*4+0] = 1;
      aBCFlag[ip*4+1] = 1;
      aBCFlag[ip*4+2] = 1;
    }
    if( py > +0.45 ){
      aVal[ip*4+0] = 10;
    }
    /*
     if( py > +0.45 ){
     aBCFlag[ip*4+0] = 1;
     aBCFlag[ip*4+1] = 1;
     aBCFlag[ip*4+2] = 1;
     aVal[ip*4+1] = 10;
     }
     */
  }
  //////
  /*
  CJaggedArray crs;
  crs.SetEdgeOfElem(aTet, (int)aTet.size()/4, 4, (int)aXYZ.size()/3, false);
  crs.Sort();
   */
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                             aTet.data(), aTet.size()/4, 4,
                             (int)aXYZ.size()/3);
  dfm2::JArray_Sort(psup_ind, psup);
  //
  mat_A.Initialize(np, 4, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  ilu_A.Initialize_ILU0(mat_A);
  //
}


void SolveProblem_Stokes_Dynamic(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    std::vector<double>& vec_b,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np*4;
  double myu = 1;
  double rho = 100.0;
  double g_x = 0.0;
  double g_y = -0.0;
  double g_z = -0.0;
  mat_A.setZero();
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_Stokes3D_Dynamic(mat_A,vec_b,
                                     myu,rho,g_x,g_y,g_z,
                                     dt_timestep,gamma_newmark,
                                     aXYZ,aTet,
                                     aVal,aVelo);
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag,0);
  // ------------------------
  std::vector<double> vec_x;
  double conv_ratio = 1.0e-4;
  int iteration = 1000;
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  vec_x.resize(vec_b.size());
  {
    const std::size_t n = vec_b.size();
    std::vector<double> tmp0(n), tmp1(n);
    Solve_PCG(
        dfm2::CVecXd(vec_b),
        dfm2::CVecXd(vec_x),
        dfm2::CVecXd(tmp0),
        dfm2::CVecXd(tmp1),
        conv_ratio, iteration, mat_A, ilu_A);
  }
  // -----------------------
  dfm2::XPlusAYBZ(
      aVal,nDoF,aBCFlag,
      dt_timestep*gamma_newmark,vec_x,
      dt_timestep,aVelo);
  dfm2::XPlusAY(
      aVelo,nDoF,aBCFlag,
      1.0,vec_x);
}


void InitializeProblem_NavierStokes_Dynamic(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  // set boundary condition
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np*4;
  // ----------
  aVal.assign(np*4, 0.0);
  aVelo.assign(np*4, 0.0);
  aBCFlag.assign(nDoF, 0);
  assert(aIsSurf.size() == np);
  for(unsigned int ip=0;ip<np;++ip){
    const double px = aXYZ[ip*3+0];
    const double py = aXYZ[ip*3+1];
    const double pz = aXYZ[ip*3+2];
    if( aIsSurf[ip] != 1 ){ continue; }
    aBCFlag[ip*4+0] = 1;
    aBCFlag[ip*4+1] = 1;
    aBCFlag[ip*4+2] = 1;
    int ishape = 0;
    if( ishape == 0 || ishape == 1 ){
      if( py > +0.48 ){
        aVal[ip*4+0] = 10;
      }
    }
    if( ishape == 2 ){
      if( fabs(pz) < 0.3-0.03 || fabs(px) < 0.3-0.03 ){
        if( py > +0.48 ){ aVal[ip*4+1] = 1; }
        if( py < -0.48 ){ 
          aBCFlag[ip*4+1] = 0;
//          aVal[ip*4+1] = 1; 
        }
      }
    }
  }
  //
  /*
  CJaggedArray crs;
  crs.SetEdgeOfElem(aTet, (int)aTet.size()/4, 4, (int)aXYZ.size()/3, false);
  crs.Sort();
   */
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                             aTet.data(), aTet.size()/4, 4,
                             (int)aXYZ.size()/3);
  dfm2::JArray_Sort(psup_ind, psup);
  //
  mat_A.Initialize(np, 4, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
//  ilu_A.Initialize_ILU0(mat_A);
  ilu_A.Initialize_ILUk(mat_A, 0);
  //
}


void SolveProblem_NavierStokes_Dynamic(
    std::vector<int>& aBCFlag,
    dfm2::CMatrixSparse<double>& mat_A,
    delfem2::CPreconditionerILU<double>& ilu_A,
    std::vector<double>& vec_b,
    //
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  const unsigned int np = aXYZ.size()/3;
  const unsigned int nDoF = np*4;
  //
  double myu = 1;
  double rho = 1000.0;
  double g_x = 0.0;
  double g_y = -0.0;
  double g_z = -0.0;
  mat_A.setZero();
  vec_b.assign(nDoF, 0.0);
  dfm2::MergeLinSys_NavierStokes3D_Dynamic(mat_A,vec_b,
                                           myu,rho,g_x,g_y,g_z,
                                           dt_timestep,gamma_newmark,
                                           aXYZ,aTet,
                                           aVal,aVelo);
  mat_A.SetFixedBC(aBCFlag.data());
  dfm2::setRHS_Zero(vec_b, aBCFlag,0);
  // --------------------------------------
  std::vector<double> vec_x;
//  SolveLinSys_PCG(mat_A,vec_b,vec_x,ilu_A,aBCFlag);
  double conv_ratio = 1.0e-5;
  int iteration = 1000;
  ilu_A.SetValueILU(mat_A);
  ilu_A.DoILUDecomp();
  vec_x.resize(vec_b.size());
  Solve_PBiCGStab(vec_b.data(),vec_x.data(),
                  conv_ratio,iteration,mat_A,ilu_A);
//  Solve_BiCGStab(mat_A,vec_b,vec_x,ilu_A, conv_ratio,iteration);
  // ----------------------------------------
  dfm2::XPlusAYBZ(aVal,nDoF,aBCFlag,
                  dt_timestep*gamma_newmark,vec_x,
                  dt_timestep,aVelo);
  dfm2::XPlusAY(aVelo,nDoF,aBCFlag,
                1.0,vec_x);
  /////
  dt_timestep = 0.0025;
}


class CInSphere : public dfm2::CInput_IsosurfaceStuffing
{
public:
  double SignedDistance(double x, double y, double z) const override {
    double n[3];
    return sdf.Projection(n,
                          x, y, z);
  }
  virtual void Level(int& ilevel_vol, int& ilevel_srf, int& nlayer, double& sdf,
                     double px, double py, double pz) const
  {
    sdf = this->SignedDistance(px, py, pz);
    ilevel_vol = -1;
    ilevel_srf = 2;
    nlayer = 3;
  }

public:
  delfem2::CSphere<double> sdf;
};

class CInBox : public dfm2::CInput_IsosurfaceStuffing
{
public:
  double SignedDistance(double x, double y, double z) const override
  {
    double n[3];
    return sdf.Projection(n,
                          x, y, z);
  }
public:
  delfem2::CBox<double> sdf;
};

void SetMesh(
    std::vector<unsigned int>& aTet,
    std::vector<double>& aXYZ,
    std::vector<int>& aIsSurf,
    int ishape)
{
  ::glMatrixMode(GL_MODELVIEW);
  
  if(ishape==0){
    const double rad = 0.5;
    CInSphere sphere;
    sphere.sdf.is_out_ = true;
    sphere.sdf.radius_ = rad;
    double cent[3] = {0,0,0};
    IsoSurfaceStuffing(aXYZ, aTet, aIsSurf, sphere, 0.30, rad*2.1, cent);
  }
  else if( ishape == 1 ){
    const double hwx = 0.5;
    const double hwy = 0.5;
    const double hwz = 0.5;
    CInBox box;
    box.sdf.hwx = hwx;
    box.sdf.hwy = hwy;
    box.sdf.hwz = hwz;
    double cent[3] = {0,0,0};
    dfm2::IsoSurfaceStuffing(aXYZ, aTet, aIsSurf, box, 0.2, 1.1, cent);
  }
  else if( ishape == 2 ){
    class CCavSphere : public dfm2::CInput_IsosurfaceStuffing
    {
    public:
      CCavSphere(){
        const double hwx = 0.3;
        const double hwy = 0.5;
        const double hwz = 0.3;
        box.sdf.hwx = hwx;
        box.sdf.hwy = hwy;
        box.sdf.hwz = hwz;
        ////
        const double rad = 0.1;
        sphere.sdf.radius_ = rad;
        sphere.sdf.is_out_ = true;
      }
      virtual double SignedDistance(double x, double y, double z ) const override {
        double dist0 = -sphere.SignedDistance(x, y, z);
        double cx = 0.0;
        double cy = 0.0;
        double cz = 0.0;
        double dist1 = box.SignedDistance(x-cx, y-cy, z-cz);
        return (dist0<dist1) ? dist0 : dist1;
      }
    public:
      CInBox box;
      CInSphere sphere;
    } cav_sphere;
    double cent[3] = {0,0,0};
    dfm2::IsoSurfaceStuffing(aXYZ, aTet, aIsSurf, cav_sphere, 0.05, 1.1, cent);
  }
}

// ---------------------------------------------

static void myGlVertex3d
(unsigned int ixyz,
 const std::vector<double>& aXYZ1 )
{
  ::glVertex3d(aXYZ1[ixyz*3+0], aXYZ1[ixyz*3+1], aXYZ1[ixyz*3+2] );
}

//static void myGlVertex3d(const CVec3& v){
//  ::glVertex3d(v.x, v.y, v.z);
//}

void myGlutDisplay(
    int iphysics,
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
//	::glClearColor(0.2f, 0.7f, 0.7f ,1.0f);
	::glClearColor(1.0f, 1.0f, 1.0f ,1.0f);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
  ::glPolygonOffset( 3.1f, 2.0f );

//  glEnable(GL_BLEND);
//  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  GLboolean is_lighting = glIsEnabled(GL_LIGHTING);
  {
    float color[4] = {200.0/256.0, 200.0/256.0, 200.0/256.0,1.0f};
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
    ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
//    glShadeModel(GL_SMOOTH);
    glShadeModel(GL_FLAT);
 }
  ::glDisable(GL_LIGHTING);
  
  if (iphysics==0 || iphysics==1){
    glShadeModel(GL_SMOOTH);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0,0,0);
    delfem2::opengl::DrawMeshTet3D_Edge(
        aXYZ.data(),aXYZ.size()/3,
        aTet.data(),aTet.size()/4);
    {
      std::vector< std::pair<double,delfem2::CColor> > colorMap;
      delfem2::ColorMap_BlueGrayRed(colorMap, 0, 1.0);
      delfem2::opengl::DrawMeshTet3D_ScalarP1(aXYZ.data(), aXYZ.size()/3,
                                              aTet.data(), aTet.size()/4,
                                              aVal.data(),
                                              colorMap);
    }
  }
  if (iphysics==2||iphysics==3){
    ::glColor3d(0,0,0);
    delfem2::opengl::DrawMeshTet3D_EdgeDisp(aXYZ.data(), aTet.data(), aTet.size()/4, aVal.data(), 1.0);
    ::glEnable(GL_LIGHTING);
    {
      float color[4] = {180.0/256.0, 180.0/256.0, 130.0/256.0,1.0f};
      ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
      ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
      glShadeModel(GL_FLAT);
    }
    delfem2::opengl::DrawMeshTet3D_FaceNormDisp(aXYZ.data(), aXYZ.size()/3,
                                       aTet.data(), aTet.size()/4,
                                       aVal.data());
  }
  if( iphysics == 4 || iphysics == 5 || iphysics == 6 ){
    ::glEnable(GL_LIGHTING);
    {
      float color[4] = {256.0/256.0, 0.0/256.0, 0.0/256.0,1.0f};
      ::glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color);
      ::glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,color);
      //    glShadeModel(GL_SMOOTH);
      glShadeModel(GL_FLAT);
    }
    for(unsigned int ip=0;ip<aXYZ.size()/3;++ip){
      const dfm2::CVec3d p(aXYZ[ip*3+0],aXYZ[ip*3+1],aXYZ[ip*3+2]);
      const dfm2::CVec3d v(aVal[ip*4+0],aVal[ip*4+1],aVal[ip*4+2]);
      double a = 0.1;
      delfem2::opengl::DrawArrow(p, a*v);
    }
  }

  
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }

  ::glColor3d(0,0,0);
}

void ProblemScalar(
    dfm2::glfw::CViewer3& viewer,
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  std::vector<int> aBCFlag;
  dfm2::CMatrixSparse<double> mat_A;
  std::vector<double> vec_b;
  delfem2::CPreconditionerILU<double>  ilu_A;
  //
  glfwSetWindowTitle(viewer.window, "Poisson");
  InitializeProblem_Poisson(
      aBCFlag,mat_A,ilu_A,
      aTet,aXYZ,aIsSurf);
  SolveProblem_Poisson(
      aBCFlag,mat_A,ilu_A,vec_b,
      aTet,aXYZ,aIsSurf);
  for(unsigned int iframe=0;iframe<50;++iframe){ // poisson
    viewer.DrawBegin_oldGL();
    myGlutDisplay(0,aTet,aXYZ,aIsSurf);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  // ---------------------------
  glfwSetWindowTitle(viewer.window, "Diffusion");
  InitializeProblem_Diffusion(
      aBCFlag,mat_A,ilu_A,
      aTet,aXYZ,aIsSurf);
  SolveProblem_Diffusion(
      aBCFlag,mat_A,ilu_A,vec_b,
      aTet,aXYZ,aIsSurf);
  for(unsigned int iframe=0;iframe<50;++iframe){
    SolveProblem_Diffusion(
        aBCFlag,mat_A,ilu_A,vec_b,
        aTet,aXYZ,aIsSurf);
    viewer.DrawBegin_oldGL();
    myGlutDisplay(1,aTet,aXYZ,aIsSurf);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
}

void ProblemSolidLinear(
    dfm2::glfw::CViewer3& viewer,
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  std::vector<int> aBCFlag;
  dfm2::CMatrixSparse<double> mat_A;
  std::vector<double> vec_b;
  delfem2::CPreconditionerILU<double>  ilu_A;
  //
  glfwSetWindowTitle(viewer.window, "SolidLinearStatic");
  InitializeProblem_ShellEigenPB(
      aBCFlag,mat_A,ilu_A,
      aTet,aXYZ,aIsSurf);
  SolveProblem_LinearSolid_Static(
      aBCFlag,mat_A,ilu_A,vec_b,
      aTet,aXYZ,aIsSurf);
  for(unsigned int iframe=0;iframe<50;++iframe){
    viewer.DrawBegin_oldGL();
    myGlutDisplay(2,aTet,aXYZ,aIsSurf);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  // ---------------------------
  glfwSetWindowTitle(viewer.window, "SolidLinearDynamic");
  InitializeProblem_LinearSolid_Dynamic(
      aBCFlag,mat_A,ilu_A,
      aTet,aXYZ,aIsSurf);
  for(unsigned int iframe=0;iframe<50;++iframe){
    SolveProblem_LinearSolid_Dynamic(
        aBCFlag,mat_A,ilu_A,vec_b,
        aTet,aXYZ,aIsSurf);
    viewer.DrawBegin_oldGL();
    myGlutDisplay(3,aTet,aXYZ,aIsSurf);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
}

void ProblemFluid(
    dfm2::glfw::CViewer3& viewer,
    const std::vector<unsigned int>& aTet,
    const std::vector<double>& aXYZ,
    const std::vector<int>& aIsSurf)
{
  std::vector<int> aBCFlag;
  dfm2::CMatrixSparse<double> mat_A;
  std::vector<double> vec_b;
  delfem2::CPreconditionerILU<double>  ilu_A;
  //
  glfwSetWindowTitle(viewer.window, "StokesStatic");
  InitializeProblem_Stokes_Static(
      aBCFlag,mat_A,ilu_A,
      aTet,aXYZ,aIsSurf);
  SolveProblem_Stokes_Static(
      aBCFlag,mat_A,ilu_A,vec_b,
      aTet,aXYZ,aIsSurf);
  for(unsigned int iframe=0;iframe<50;++iframe){
    viewer.DrawBegin_oldGL();
    myGlutDisplay(4,aTet,aXYZ,aIsSurf);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  // ---------------------------
  glfwSetWindowTitle(viewer.window, "StokesDynamic");
  InitializeProblem_Stokes_Dynamic(
      aBCFlag,mat_A,ilu_A,
      aTet,aXYZ,aIsSurf);
  for(unsigned int iframe=0;iframe<50;++iframe){
    SolveProblem_Stokes_Dynamic(
        aBCFlag,mat_A,ilu_A,vec_b,
        aTet,aXYZ,aIsSurf);
    viewer.DrawBegin_oldGL();
    myGlutDisplay(5,aTet,aXYZ,aIsSurf);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
  // ---------------------------
  glfwSetWindowTitle(viewer.window, "NavierStokesDynamic");
  InitializeProblem_NavierStokes_Dynamic(
      aBCFlag,mat_A,ilu_A,
      aTet,aXYZ,aIsSurf);
  for(unsigned int iframe=0;iframe<50;++iframe){
    SolveProblem_NavierStokes_Dynamic(
        aBCFlag,mat_A,ilu_A,vec_b,
        aTet,aXYZ,aIsSurf);
    viewer.DrawBegin_oldGL();
    myGlutDisplay(6,aTet,aXYZ,aIsSurf);
    viewer.SwapBuffers();
    glfwPollEvents();
    viewer.ExitIfClosed();
  }
}

int main(int argc,char* argv[])
{
  dfm2::glfw::CViewer3 viewer;
  dfm2::glfw::InitGLOld();
  viewer.InitGL();
  viewer.camera.view_height = 1.5;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();
  while (true){
    std::vector<unsigned int> aTet;
    std::vector<double> aXYZ;
    std::vector<int> aIsSurf;
    SetMesh(
        aTet,aXYZ,aIsSurf,
        0);
    ProblemScalar(viewer,aTet,aXYZ,aIsSurf);
    ProblemSolidLinear(viewer, aTet,aXYZ,aIsSurf);
    ProblemFluid(viewer, aTet,aXYZ,aIsSurf);
  }
}
