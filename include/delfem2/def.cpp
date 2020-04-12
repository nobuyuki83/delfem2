/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/def.h"
#include "delfem2/mat3.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vecxitrsol.h"
#include <cstdio>

namespace delfem2 {
namespace def {

void SetLinSys_LaplaceGraph_MeshTri3
(CMatrixSparse<double>& mat_A)
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

class CSysMat_DefLaplaceLinear{
public:
  CSysMat_DefLaplaceLinear(const CMatrixSparse<double>& A0,
                           const std::vector<int>& aFlg0,
                           double weight_bc0)
  : A(A0), aBCFlag(aFlg0), weight_bc(weight_bc0)
  {
    const unsigned int np = aFlg0.size()/3;
    vec_tmp.resize(np*3);
    
    // make jacobi preconditioner
    aDiaInv.assign(np*9,0.0);
    for(unsigned int ip=0;ip<np;++ip){
      for(unsigned int icrs=A.colInd[ip];icrs<A.colInd[ip+1];++icrs){
        unsigned int jp0 = A.rowPtr[icrs];
        MatTMat3_ScaleAdd(aDiaInv.data()+jp0*9,
                          A.valCrs.data()+icrs*9,
                          A.valCrs.data()+icrs*9,
                          1.0, 0.0); // del. prev. value and set new vaue
      }
      {
        MatTMat3_ScaleAdd(aDiaInv.data()+ip*9,
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
      Inverse_Mat3(aDiaInv.data()+ip*9);
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
  void Solve(double* v) const { // for preconditioner
    const unsigned int np = aBCFlag.size()/3;
    for(int ip=0;ip<np;++ip){
      double tmp[3];
      MatVec3(tmp, aDiaInv.data()+ip*9, v+ip*3);
      v[ip*3+0] = tmp[0];
      v[ip*3+1] = tmp[1];
      v[ip*3+2] = tmp[2];
    }
  }
public:
  const CMatrixSparse<double>& A;
  const std::vector<int>& aBCFlag;
  const double weight_bc;
  std::vector<double> aDiaInv;
  mutable std::vector<double> vec_tmp;
};

}
}


// ============================================

void delfem2::CDef_SingleLaplacian::Init
 (const std::vector<double>& aXYZ0,
  const std::vector<unsigned int>& aTri)
{
  std::vector<unsigned int> psup_ind, psup;
  JArray_PSuP_MeshElem(psup_ind, psup,
                       aTri.data(), aTri.size()/3, 3,
                       (int)aXYZ0.size()/3);
  JArray_Sort(psup_ind, psup);
  mat_A.Initialize(aXYZ0.size()/3, 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  // ---
  def::SetLinSys_LaplaceGraph_MeshTri3(mat_A);
  aRhs0.resize(aXYZ0.size());
  mat_A.MatVec(aRhs0.data(),
               1.0, aXYZ0.data(), 0.0);
}

void delfem2::CDef_SingleLaplacian::Solve
 (std::vector<double>& aXYZ1,
  const std::vector<double>& aXYZ0,
  const std::vector<int>& aBCFlag)
{    // ----------
  aRhs1 = aRhs0;
  for(int i=0;i<aBCFlag.size();++i){
    if( aBCFlag[i] == 0 ){ continue; }
    aRhs1[i] = aXYZ1[i];
  }
  mat_A.SetFixedBC_Dia(aBCFlag.data(), 1.0);
  mat_A.SetFixedBC_Row(aBCFlag.data());
  aXYZ1 = aXYZ0;
  aHistConv = Solve_BiCGStab(aRhs1, aXYZ1,
                             1.0e-5, 100, mat_A);
}


// ----------------------------

void delfem2::CDef_LaplacianLinear::Init
 (const std::vector<double>& aXYZ0,
  const std::vector<unsigned int>& aTri,
  bool is_preconditioner)
{
  this->is_preconditioner = is_preconditioner;
  std::vector<unsigned int> psup_ind, psup;
  JArray_PSuP_MeshElem(psup_ind, psup,
                       aTri.data(), aTri.size()/3, 3,
                       (int)aXYZ0.size()/3);
  JArray_Sort(psup_ind, psup);
  mat_A.Initialize(aXYZ0.size()/3, 3, true);
  mat_A.SetPattern(psup_ind.data(), psup_ind.size(),
                   psup.data(),     psup.size());
  def::SetLinSys_LaplaceGraph_MeshTri3(mat_A);
}

void delfem2::CDef_LaplacianLinear::Solve
 (std::vector<double>& aXYZ1,
  const std::vector<double>& aXYZ0,
  const std::vector<int>& aBCFlag)
{
  def::SetLinSys_LaplaceGraph_MeshTri3(mat_A);
  const double weight_bc = 100.0;
  def::CSysMat_DefLaplaceLinear mat_AtA(mat_A, aBCFlag, weight_bc);
  // ----------
  std::vector<double> aRhs(aXYZ0.size(),0.0);
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
  /*
   { // making RHS vector for fixed boundary condition
   const unsigned int np = aXYZ0.size()/3;
   std::vector<double> aGoal(np*3);
   SetPositionAtFixedBoundary(aGoal,
   iframe,aXYZ0,aBCFlag);
   for(int i=0;i<np*3;++i){
   if( aBCFlag[i] == 0 ){ continue; }
   aRhs[i] += (aGoal[i]-aXYZ1[i])*weight_bc;
   }
   }
   */
  std::vector<double> aUpd(aXYZ0.size(),0.0);
  if( is_preconditioner ){
    std::vector<double> aRes = Solve_PCG(aRhs.data(), aUpd.data(),
                                         aRhs.size(), 1.0e-7, 300, mat_AtA, mat_AtA);
    std::cout << aRes.size() << std::endl;
  }
  else{
    std::vector<double> aRes = Solve_CG(aRhs.data(), aUpd.data(),
                                        aRhs.size(), 1.0e-7, 300, mat_AtA);
    std::cout << aRes.size() << std::endl;
  }
  for(int i=0;i<aBCFlag.size();++i){ aXYZ1[i] += aUpd[i]; }
}
