/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/def.h"
#include "delfem2/mat3.h"
#include "delfem2/quat.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vecxitrsol.h"
#include <cstdio>

namespace delfem2 {
namespace def {

template <typename REAL>
DFM2_INLINE void Add3
 (REAL vo[3],
  const REAL vi[3])
{
  vo[0] += vi[0];
  vo[1] += vi[1];
  vo[2] += vi[2];
}
#ifdef DFM2_STATIC_LIBRARY
template void Add3(float vo[3], const float vi[3]);
template void Add3(double vo[3], const double vi[3]);
#endif

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

}
}
 


// ============================================

void delfem2::CDef_SingleLaplacianDisponly::Init
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

void delfem2::CDef_SingleLaplacianDisponly::Deform
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

void delfem2::CDef_LaplacianDisponly::Init
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

void delfem2::CDef_LaplacianDisponly::MakeLinearSystem()
{
  const unsigned int np = mat_A.nblk_col;
  vec_tmp.resize(np*3);
  
  // make jacobi preconditioner
  aDiaInv.assign(np*9,0.0);
  for(unsigned int ip=0;ip<np;++ip){
    for(unsigned int icrs=mat_A.colInd[ip];icrs<mat_A.colInd[ip+1];++icrs){
      unsigned int jp0 = mat_A.rowPtr[icrs];
      MatTMat3_ScaleAdd(aDiaInv.data()+jp0*9,
                        mat_A.valCrs.data()+icrs*9,
                        mat_A.valCrs.data()+icrs*9,
                        1.0, 0.0); // del. prev. value and set new vaue
    }
    {
      MatTMat3_ScaleAdd(aDiaInv.data()+ip*9,
                        mat_A.valDia.data()+ip*9,
                        mat_A.valDia.data()+ip*9,
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

void delfem2::CDef_LaplacianDisponly::Deform
 (std::vector<double>& aXYZ1,
  const std::vector<double>& aXYZ0,
  const std::vector<int>& aBCFlag)
{
  weight_bc = 100.0;
  this->aBCFlag = aBCFlag;
  this->MakeLinearSystem();
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
                                         aRhs.size(), 1.0e-7, 300, *this, *this);
    std::cout << "pcg: " << aRes.size() << std::endl;
  }
  else{
    std::vector<double> aRes = Solve_CG(aRhs.data(), aUpd.data(),
                                        aRhs.size(), 1.0e-7, 300, *this);
    std::cout << "cg: " << aRes.size() << std::endl;
  }
  for(int i=0;i<aBCFlag.size();++i){ aXYZ1[i] += aUpd[i]; }
}

void delfem2::CDef_LaplacianDisponly::MatVec
 (double* y,
  double alpha, const double* vec,  double beta) const {
  mat_A.MatVec(vec_tmp.data(),
               1, vec, 0.0);
  mat_A.MatTVec(y,
                alpha, vec_tmp.data(), beta);
  // add diagonal for fixed boundary condition
  for(int i=0;i<aBCFlag.size();++i){
    if( aBCFlag[i] == 0 ){ continue; }
    y[i] += weight_bc*vec[i];
  }
}

// for preconditioner
void delfem2::CDef_LaplacianDisponly::Solve(double* v) const
{
  const unsigned int np = aBCFlag.size()/3;
  for(int ip=0;ip<np;++ip){
    double tmp[3];
    MatVec3(tmp, aDiaInv.data()+ip*9, v+ip*3);
    v[ip*3+0] = tmp[0];
    v[ip*3+1] = tmp[1];
    v[ip*3+2] = tmp[2];
  }
}


// ======================================

delfem2::CDef_ArapEdgeLinearDisponly::CDef_ArapEdgeLinearDisponly
 (const std::vector<double>& aXYZ0,
  const std::vector<unsigned int>& aTri,
  double weight_bc0,
  const std::vector<int>& aBCFlag0) :
weight_bc(weight_bc0),
aBCFlag(aBCFlag0)
{
  const unsigned int np = aXYZ0.size()/3;
  JArray_PSuP_MeshElem(psup_ind, psup,
                       aTri.data(), aTri.size()/3, 3,
                       (int)aXYZ0.size()/3);
  JArray_Sort(psup_ind, psup);
  // ------
  const unsigned int ne = psup.size();
  // -----
  aMatEdge.resize(ne*9*2);
  assert(psup_ind.size()==np+1);
  for(unsigned int ip=0;ip<np;++ip){
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      //          unsigned int jp = psup[ipsup];
      Mat3_Identity(aMatEdge.data()+ipsup*18,   +1.0);
      Mat3_Identity(aMatEdge.data()+ipsup*18+9, -1.0);
    }
  }
  // ---
  vec_tmp.resize(ne*3);
}

void delfem2::CDef_ArapEdgeLinearDisponly::JacobiTVecTmp
 (double*y ,
  double alpha, double beta) const
{
  const unsigned int np = aBCFlag.size()/3;
  for(int i=0;i<np*3;++i){ y[i] *= beta; }
  for(unsigned int ip=0;ip<np;++ip){
    for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      unsigned int jp0 = psup[ipsup];
      MatTVec3_ScaleAdd(y+ip*3,
                        aMatEdge.data()+ipsup*18,
                        vec_tmp.data()+ipsup*3,
                        alpha, 1.0);
      MatTVec3_ScaleAdd(y+jp0*3,
                        aMatEdge.data()+ipsup*18+9,
                        vec_tmp.data()+ipsup*3,
                        alpha, 1.0);
    }
  }
}

void delfem2::CDef_ArapEdgeLinearDisponly::MakeLinearSystem
 (double* aRhs,
  const double* aXYZ0,
  const double* aXYZ1) const
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
  /*
  // making RHS vector for fixed boundary condition
  for(int i=0;i<np*3;++i){
    if( aBCFlag[i] == 0 ){ continue; }
    aRhs[i] += (aGoal[i]-aXYZ1[i])*weight_bc;
  }
   */
}

void delfem2::CDef_ArapEdgeLinearDisponly::MatVec
 (double* y,
  double alpha,
  const double* vec,
  double beta) const
{
  const unsigned int np = aBCFlag.size()/3;
  std::fill(vec_tmp.begin(),vec_tmp.end(), 0.0);
  for(unsigned int ip=0;ip<np;++ip){
    for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      unsigned int jp0 = psup[ipsup];
      MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
                       aMatEdge.data()+ipsup*18,
                       vec+ip*3,
                       1.0, 1.0);
      MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
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

void delfem2::CDef_ArapEdgeLinearDisponly::Deform
(std::vector<double>& aXYZ1,
 const std::vector<double>& aXYZ0)
{
  const unsigned int np = aBCFlag.size()/3;
  std::vector<double> aRhs(np*3,0.0);
  this->MakeLinearSystem(aRhs.data(),
                        aXYZ0.data(), aXYZ1.data());
  std::vector<double> aUpd(np*3,0.0);
  std::vector<double> aRes = Solve_CG(aRhs.data(), aUpd.data(),
                                      np*3, 1.0e-4, 300, *this);
//  std::cout << "iframe: " << iframe << "   nitr:" << aRes.size() << std::endl;
  for(int i=0;i<np*3;++i){ aXYZ1[i] += aUpd[i]; }
}

// ======================================================

void delfem2::CDef_ArapEdge::Init
 (const std::vector<double>& aXYZ0,
  const std::vector<unsigned int>& aTri,
  double weight_bc0,
  const std::vector<int>& aBCFlag,
  bool is_preconditioner0)
{
  this->weight_bc = weight_bc0;
  this->is_preconditioner = is_preconditioner0;
  this->aBCFlag = aBCFlag;
  const unsigned int np = aXYZ0.size()/3;
  JArray_PSuP_MeshElem(psup_ind, psup,
                       aTri.data(), aTri.size()/3, 3,
                       (int)aXYZ0.size()/3);
  JArray_Sort(psup_ind, psup);
  // ---------
  const unsigned int ne = psup.size();
  // -----
  aMatEdge.resize(ne*27);
  assert(psup_ind.size()==np+1);
  for(unsigned int ip=0;ip<np;++ip){
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      //          unsigned int jp = psup[ipsup];
      Mat3_Identity(aMatEdge.data()+ipsup*27+0, +1.0);
      Mat3_Identity(aMatEdge.data()+ipsup*27+9, -1.0);
    }
  }
  // ---
  vec_tmp.resize(ne*3);
}

void delfem2::CDef_ArapEdge::JacobiTVecTmp
 (double*y ,
  double alpha,
  double beta) const
{
  const unsigned int np = psup_ind.size()-1;
  for(int i=0;i<np*6;++i){ y[i] *= beta; }
  for(unsigned int ip=0;ip<np;++ip){
    for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      unsigned int jp0 = psup[ipsup];
      MatTVec3_ScaleAdd(y+ip*3,
                        aMatEdge.data()+ipsup*27+0,
                        vec_tmp.data()+ipsup*3,
                        alpha, 1.0);
      MatTVec3_ScaleAdd(y+jp0*3,
                        aMatEdge.data()+ipsup*27+9,
                        vec_tmp.data()+ipsup*3,
                        alpha, 1.0);
      MatTVec3_ScaleAdd(y+np*3+ip*3,
                        aMatEdge.data()+ipsup*27+18,
                        vec_tmp.data()+ipsup*3,
                        alpha, 1.0);
    }
  }
}

void delfem2::CDef_ArapEdge::MatVec
 (double* y,
  double alpha,
  const double* vec,
  double beta) const
{
  const unsigned int np = psup_ind.size()-1;
  std::fill(vec_tmp.begin(),vec_tmp.end(), 0.0);
  for(unsigned int ip=0;ip<np;++ip){
    for(int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      unsigned int jp0 = psup[ipsup];
      MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
                       aMatEdge.data()+ipsup*27+0,
                       vec+ip*3,
                       1.0, 1.0);
      MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
                       aMatEdge.data()+ipsup*27+9,
                       vec+jp0*3,
                       1.0, 1.0);
      MatVec3_ScaleAdd(vec_tmp.data()+ipsup*3,
                       aMatEdge.data()+ipsup*27+18,
                       vec+(np+ip)*3,
                       1.0, 1.0);
    }
  }
  this->JacobiTVecTmp(y,
                      alpha, beta);
  // add diagonal for fixed boundary condition
  for(int i=0;i<aBCFlag.size();++i){
    if( aBCFlag[i] == 0 ){ continue; }
    y[i] += weight_bc*vec[i];
    //      y[np*3+i] += weight_bc*vec[np*3+i];
  }
}

void delfem2::CDef_ArapEdge::MakeLinearSystem
 (double* aRhs,
  const double* aXYZ0,
  const double* aXYZ1,
  const double* aQuat)
{
  const unsigned int np = psup_ind.size()-1;
  const unsigned int ne = psup.size();
  vec_tmp.assign(ne*3,0);
  for(unsigned int ip=0;ip<np;++ip){
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const unsigned int jp0 = psup[ipsup];
      const double* q0 = aQuat+ip*4;
      const double d0[3] = { aXYZ0[jp0*3+0]-aXYZ0[ip*3+0], aXYZ0[jp0*3+1 ]-aXYZ0[ip*3+1], aXYZ0[jp0*3+2]-aXYZ0[ip*3+2] };
      const double d1[3] = { aXYZ1[jp0*3+0]-aXYZ1[ip*3+0], aXYZ1[jp0*3+1]-aXYZ1[ip*3+1], aXYZ1[jp0*3+2]-aXYZ1[ip*3+2] };
      double Rd0[3]; QuatVec(Rd0, q0,d0);
      vec_tmp[ipsup*3+0] += +(Rd0[0] - d1[0]);
      vec_tmp[ipsup*3+1] += +(Rd0[1] - d1[1]);
      vec_tmp[ipsup*3+2] += +(Rd0[2] - d1[2]);
      Mat3_Spin_ScaleAdd(aMatEdge.data()+ipsup*27+18,
                         Rd0,
                         -1.0, 0.0);
    }
  }
  this->JacobiTVecTmp(aRhs,
                      -1.0, 0.0);
  /*
  // making RHS vector for fixed boundary condition
  for(int i=0;i<np*3;++i){
    if( aBCFlag[i] == 0 ){ continue; }
    aRhs[i] += (aGoal[i]-aXYZ1[i])*weight_bc;
    //aRhs[i+np*3] = 0.0;
  }
   */
}

void delfem2::CDef_ArapEdge::MakePreconditionerJacobi()
{
  const unsigned int np = psup_ind.size()-1;
  aDiaInv.assign(np*2*9, 0.0);
  for(unsigned int ip=0;ip<np;++ip){
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup) {
      const unsigned int jp0 = psup[ipsup];
      MatTMat3_ScaleAdd(aDiaInv.data()+ip*9,
                        aMatEdge.data()+ipsup*27+0,
                        aMatEdge.data()+ipsup*27+0,
                        1.0,1.0);
      MatTMat3_ScaleAdd(aDiaInv.data()+jp0*9,
                        aMatEdge.data()+ipsup*27+9,
                        aMatEdge.data()+ipsup*27+9,
                        1.0,1.0);
      MatTMat3_ScaleAdd(aDiaInv.data()+(np+ip)*9,
                        aMatEdge.data()+ipsup*27+18,
                        aMatEdge.data()+ipsup*27+18,
                        1.0,1.0);
    }
  }
  for(unsigned int ip=0;ip<np;++ip){
    for(int idim=0;idim<3;++idim) {
      if (aBCFlag[ip*3+idim] == 0) { continue; }
      aDiaInv[ip*9+idim*3+idim] += weight_bc;
    }
  }
  for(unsigned int ip=0;ip<np*2;++ip){
    Inverse_Mat3(aDiaInv.data()+ip*9);
  }
}

void delfem2::CDef_ArapEdge::Solve(double* v) const
{
  const unsigned int np = psup_ind.size()-1;
  for(int ip=0;ip<np*2;++ip){
    double tmp[3];
    MatVec3(tmp, aDiaInv.data()+ip*9, v+ip*3);
    v[ip*3+0] = tmp[0];
    v[ip*3+1] = tmp[1];
    v[ip*3+2] = tmp[2];
  }
}


void delfem2::CDef_ArapEdge::Deform(std::vector<double>& aXYZ1,
                                std::vector<double>& aQuat,
                                const std::vector<double>& aXYZ0)
{
  const unsigned int np = psup_ind.size()-1;
  std::vector<double> aRhs(np*6,0.0);
  this->MakeLinearSystem(aRhs.data(),
                        aXYZ0.data(), aXYZ1.data(), aQuat.data());
  std::vector<double> aUpd(np*6,0.0);
  std::vector<double> aConvHist;
  if( is_preconditioner ){
    this->MakePreconditionerJacobi();
    aConvHist = Solve_PCG(aRhs.data(), aUpd.data(),
                          np*6, 1.0e-4, 400, *this, *this);
  }
  else{
    aConvHist = Solve_CG(aRhs.data(), aUpd.data(),
                         np*6, 1.0e-4, 400, *this);
  }
  for(int ip=0;ip<np;++ip){
    def::Add3(aXYZ1.data()+ip*3, aUpd.data()+ip*3);
    double q0[4]; Quat_CartesianAngle(q0, aUpd.data()+np*3+ip*3);
    double q1[4]; QuatQuat(q1, q0, aQuat.data()+ip*4);
    Copy_Quat(aQuat.data()+ip*4, q1);
  }
}
