/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/def.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/geo3_v23m34q.h" // update rotation by matching cluster

namespace delfem2 {
namespace def {

void SetLinSys_LaplaceGraph_MeshTri3
(CMatrixSparse<double>& mat_A)
{
  mat_A.SetZero();
  for(unsigned int ip=0;ip<mat_A.nblk_col;++ip){
    const double dn = (double)(mat_A.colInd[ip+1] - mat_A.colInd[ip]);
    for(unsigned int icrs=mat_A.colInd[ip];icrs<mat_A.colInd[ip+1];++icrs){
      mat_A.valCrs[icrs*9+0*3+0] = -1.0;
      mat_A.valCrs[icrs*9+1*3+1] = -1.0;
      mat_A.valCrs[icrs*9+2*3+2] = -1.0;
    }
    mat_A.valDia[ip*9+0*3+0] = dn;
    mat_A.valDia[ip*9+1*3+1] = dn;
    mat_A.valDia[ip*9+2*3+2] = dn;
  }
}

DFM2_INLINE void dWddW_ArapEnergy
(std::vector<double>& eM,
 std::vector<double>& eR,
 const double* Minv,
 const std::vector<unsigned int>& aIP,
 const std::vector<double>& aXYZ0,
 const std::vector<double>& aXYZ1,
 const std::vector<double>& aQuat1)
{
  const unsigned int nIP = aIP.size();
  const unsigned int nNg = nIP-1; // number of neighbor
  unsigned int ip = aIP[nNg];
  const CVec3d Pi(aXYZ0.data()+ip*3);
  const CMat3d LMi(Minv);
  const CMat3d Mrot = CMat3d::Quat(aQuat1.data()+ip*4);
  eM.assign(nIP*nIP*9, 0.0);
  for(unsigned int jjp=0;jjp<nNg;++jjp){
    for(unsigned int kkp=0;kkp<nNg;++kkp){
      const CVec3d vj = (CVec3d(aXYZ0.data()+aIP[jjp]*3)-Pi);
      const CVec3d vk = (CVec3d(aXYZ0.data()+aIP[kkp]*3)-Pi);
      CMat3d L1 = Mrot*CMat3d::Spin(vk.p)*LMi*CMat3d::Spin(vj.p)*Mrot.Trans();
      L1.AddToScale(eM.data()+(kkp*nIP+jjp)*9, -1.0);
      L1.AddToScale(eM.data()+(nNg*nIP+nNg)*9, -1.0);
      L1.AddToScale(eM.data()+(nNg*nIP+jjp)*9, +1.0);
      L1.AddToScale(eM.data()+(kkp*nIP+nNg)*9, +1.0);
    }
    {
      CMat3d L1 = CMat3d::Identity();
      L1.AddToScale(eM.data()+(jjp*nIP+jjp)*9, +1.0);
      L1.AddToScale(eM.data()+(nNg*nIP+nNg)*9, +1.0);
      L1.AddToScale(eM.data()+(nNg*nIP+jjp)*9, -1.0);
      L1.AddToScale(eM.data()+(jjp*nIP+nNg)*9, -1.0);
    }
  }
  
  //
  eR.assign(nIP*3, 0.0);
  const CVec3d pi(aXYZ1.data()+ip*3);
  CMat3d LM; LM.SetZero();
  for(unsigned int jjp=0;jjp<nNg;++jjp){
    const unsigned int jp = aIP[jjp];
    const CVec3d v0 = Mrot*(CVec3d(aXYZ0.data()+jp*3)-Pi);
    CVec3d pj(aXYZ1.data()+jp*3);
    const CVec3d v1 = pj-pi;
    const CVec3d r = -(v1-v0);
    r.AddToScale(eR.data()+nNg*3, +1);
    r.AddToScale(eR.data()+jjp*3, -1);
  }
}

}
}
 
// ==================================================

void delfem2::CDef_LaplacianLinearAsym::Init
 (const std::vector<double>& aXYZ0,
  const std::vector<unsigned int>& aTri)
{
  std::vector<unsigned int> psup_ind, psup;
  JArray_PSuP_MeshElem(psup_ind, psup,
                       aTri.data(),
                       aTri.size()/3, 3,
                       aXYZ0.size()/3);
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

void delfem2::CDef_LaplacianLinearAsym::Deform
 (std::vector<double>& aXYZ1,
  const std::vector<double>& aXYZ0,
  const std::vector<int>& aBCFlag)
{    // ----------
  aRhs1 = aRhs0;
  for(unsigned int i=0;i<aBCFlag.size();++i){
    if( aBCFlag[i] == 0 ){ continue; }
    aRhs1[i] = aXYZ1[i];
  }
  mat_A.SetFixedBC_Dia(aBCFlag.data(), 1.0);
  mat_A.SetFixedBC_Row(aBCFlag.data());
  aXYZ1 = aXYZ0;
  aHistConv = Solve_BiCGStab(aRhs1, aXYZ1,
                             1.0e-5, 100, mat_A);
}


// ===================================================
// below: implementation of CDef_LaplacianLinearGram

void delfem2::CDef_LaplacianLinearGram::Init
 (const std::vector<double>& aXYZ0,
  const std::vector<unsigned int>& aTri,
  bool is_preconditioner0)
{
  this->is_preconditioner = is_preconditioner0;
  std::vector<unsigned int> psup_ind, psup;
  JArray_PSuP_MeshElem(psup_ind, psup,
                       aTri.data(), aTri.size()/3, 3,
                       aXYZ0.size()/3);
  JArray_Sort(psup_ind, psup);
  Mat.Initialize(aXYZ0.size()/3, 3, true);
  Mat.SetPattern(psup_ind.data(), psup_ind.size(),
                 psup.data(),     psup.size());
  def::SetLinSys_LaplaceGraph_MeshTri3(Mat);
  const unsigned int np = Mat.nblk_col;
  aRes0.assign(np*3,0.0);
  Mat.MatVec(aRes0.data(),
             -1.0, aXYZ0.data(), 0.0);
}

void delfem2::CDef_LaplacianLinearGram::SetBoundaryCondition
(const std::vector<int>& aBCFlag_)
{
  this->aBCFlag = aBCFlag_;
  if( !is_preconditioner ){ return; }
  // ---------
  // make jacobi preconditioner
  const unsigned int np = Mat.nblk_col;
  aDiaInv.assign(np*9,0.0);
  for(unsigned int ip=0;ip<np;++ip){
    for(unsigned int icrs=Mat.colInd[ip];icrs<Mat.colInd[ip+1];++icrs){
      unsigned int jp0 = Mat.rowPtr[icrs];
      MatTMat3_ScaleAdd(aDiaInv.data()+jp0*9,
                        Mat.valCrs.data()+icrs*9,
                        Mat.valCrs.data()+icrs*9,
                        1.0, 1.0); // del. prev. value and set new vaue
    }
    {
      MatTMat3_ScaleAdd(aDiaInv.data()+ip*9,
                        Mat.valDia.data()+ip*9,
                        Mat.valDia.data()+ip*9,
                        1.0, 1.0); // del. prev. value and set new vaue
    }
  }
  for(unsigned int ip=0;ip<np;++ip){
    for(int i=0;i<3;++i){
      if( aBCFlag[ip*3+i] == 0 ){ continue; }
      aDiaInv[ip*9+i*3+i] += weight_bc;
    }
  }
  for(unsigned int ip=0;ip<np;++ip){
    Inverse_Mat3(aDiaInv.data()+ip*9);
  }
}

void delfem2::CDef_LaplacianLinearGram::Deform
 (std::vector<double>& aXYZ1,
  const std::vector<double>& aXYZ0) const
{
  vec_tmp0.resize(aXYZ0.size());
  vec_tmp1.resize(aXYZ0.size());
  vec_tmp2.resize(aXYZ0.size());
  //
  std::vector<double>& aRhs = vec_tmp0;
  { // making RHS vector for elastic deformation
    vec_tmp2 = aRes0;
    Mat.MatVec(vec_tmp2.data(),
                 +1.0, aXYZ1.data(), 1.0);
    Mat.MatTVec(aRhs.data(),
                -1.0, vec_tmp2.data(), 0.0);
  }
  std::vector<double>& aUpd = vec_tmp1;
  aUpd.assign(aXYZ0.size(),0.0);
  if( is_preconditioner ){
    aConvHist = Solve_PCG(aRhs.data(), aUpd.data(),
                          aRhs.size(), 1.0e-7, 300, *this, *this);
  }
  else{
    aConvHist = Solve_CG(aRhs.data(), aUpd.data(),
                         aRhs.size(), 1.0e-7, 300, *this);
  }
  for(unsigned int i=0;i<aBCFlag.size();++i){ aXYZ1[i] += aUpd[i]; }
}

void delfem2::CDef_LaplacianLinearGram::MatVec
 (double* y,
  double alpha, const double* vec,  double beta) const
{
  Mat.MatVec(vec_tmp2.data(),
               1, vec, 0.0);
  Mat.MatTVec(y,
              alpha, vec_tmp2.data(), beta);
  // add diagonal for fixed boundary condition
  for(unsigned int i=0;i<aBCFlag.size();++i){
    if( aBCFlag[i] == 0 ){ continue; }
    y[i] += weight_bc*vec[i];
  }
}

// for preconditioner
void delfem2::CDef_LaplacianLinearGram::SolvePrecond(double* v) const
{
  const unsigned int np = aBCFlag.size()/3;
  for(unsigned int ip=0;ip<np;++ip){
    double tmp[3];
    MatVec3(tmp, aDiaInv.data()+ip*9, v+ip*3);
    v[ip*3+0] = tmp[0];
    v[ip*3+1] = tmp[1];
    v[ip*3+2] = tmp[2];
  }
}

// above: delfem2::CDef_LaplacianLinearGram
// =========================================================================
// below: delfem2::CDef_LaplacianLinear


namespace delfem2 {

void Hoge
(std::vector<double>& eM,
 const std::vector<unsigned int>& aIP,
 const std::vector<double>& aXYZ0)
{
  const unsigned int nIP = aIP.size();
  const unsigned int nNg = nIP-1; // number of neighbor
  double dn = (double)nNg;
  eM.assign(nIP*nIP*9, 0.0);
  const CMat3d L1 = CMat3d::Identity();
  L1.AddToScale(eM.data()+(nNg*nIP+nNg)*9, +dn*dn);
  for(unsigned int jjp=0;jjp<nNg;++jjp){
    L1.AddToScale(eM.data()+(nNg*nIP+jjp)*9, -dn);
    L1.AddToScale(eM.data()+(jjp*nIP+nNg)*9, -dn);
    for(unsigned int kkp=0;kkp<nNg;++kkp){
      L1.AddToScale(eM.data()+(jjp*nIP+kkp)*9, +1.0);
    }
  }
}

}

void delfem2::CDef_LaplacianLinear::Init
(const std::vector<double>& aXYZ0,
 const std::vector<unsigned int>& aTri,
 bool is_preconditioner_)
{
  const unsigned int np = aXYZ0.size()/3;
  this->is_preconditioner = is_preconditioner_;
  std::vector<unsigned int> psup_ind, psup;
  JArray_PSuP_MeshElem(psup_ind, psup,
                       aTri.data(), aTri.size()/3, 3,
                       np);
  JArray_Sort(psup_ind, psup);
  {
    std::vector<unsigned int> psup_ind1, psup1;
    JArray_Extend(psup_ind1, psup1,
                  psup_ind, psup);
    JArray_Sort(psup_ind1, psup1);
    Mat.Initialize(np, 3, true);
    assert( psup_ind1.size() == np+1 );
    Mat.SetPattern(psup_ind1.data(), psup_ind1.size(),
                   psup1.data(), psup1.size());
  }
    
  std::vector<int> tmp_buffer;
  Mat.SetZero();
  for(unsigned int ip=0;ip<np;++ip){
    std::vector<unsigned int> aIP;
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      aIP.push_back(psup[ipsup]);
    }
    aIP.push_back(ip);
    std::vector<double> eM;
    delfem2::Hoge(eM, aIP, aXYZ0);
    Mat.Mearge(aIP.size(), aIP.data(),
               aIP.size(), aIP.data(),
               9, eM.data(),
               tmp_buffer);
  }
  
  aRes0.resize(aXYZ0.size());
  Mat.MatVec(aRes0.data(),
             -1.0, aXYZ0.data(), 0.0);
  
  this->Prec.Clear();
  if( is_preconditioner ){
    this->Prec.Initialize_ILU0(Mat);
  }
}


void delfem2::CDef_LaplacianLinear::SetBoundaryCondition(const std::vector<int> &aBCFlag_)
{
  this->aBCFlag = aBCFlag_;
  if( !is_preconditioner ){ return; }
  //
  const unsigned int np = Mat.nblk_col;
  assert( aBCFlag.size() == np*3 );
  for(int ip=0;ip<np;++ip){
    for(int idim=0;idim<3;++idim){
      if( aBCFlag[ip*3+idim] == 0 ){ continue; }
      Mat.valDia[ip*9+idim*3+idim] += weight_bc;
    }
  }
  this->Prec.SetValueILU(Mat);
  this->Prec.DoILUDecomp();
  for(int ip=0;ip<np;++ip){
    for(int idim=0;idim<3;++idim){
      if( aBCFlag[ip*3+idim] == 0 ){ continue; }
      Mat.valDia[ip*9+idim*3+idim] -= weight_bc;
    }
  }
}



void delfem2::CDef_LaplacianLinear::Deform
(std::vector<double>& aXYZ1,
 const std::vector<double>& aXYZ0) const
{
  // ----------
  vec_tmp0.resize(aXYZ0.size());
  vec_tmp1.resize(aXYZ0.size());
  //
  std::vector<double>& aRhs = vec_tmp0;
  memcpy(aRhs.data(), aRes0.data(), aRes0.size()*sizeof(double) );
  Mat.MatVec(aRhs.data(),
             -1.0, aXYZ1.data(), -1.0);
  std::vector<double>& aUpd = vec_tmp1;
  aUpd.assign(aXYZ0.size(),0.0);
  if( is_preconditioner ){
    aConvHist = Solve_PCG(aRhs.data(), aUpd.data(),
                          aRhs.size(), 1.0e-7, 300, *this, Prec);
  }
  else{
    aConvHist = Solve_CG(aRhs.data(), aUpd.data(),
                         aRhs.size(), 1.0e-7, 300, *this);
  }
  for(unsigned int i=0;i<aBCFlag.size();++i){ aXYZ1[i] += aUpd[i]; }
}

void delfem2::CDef_LaplacianLinear::MatVec(
    double* y,
    double alpha, const double* vec,  double beta) const
{
  Mat.MatTVec(y,
              alpha, vec, beta);
  // add diagonal for fixed boundary condition
  for(unsigned int i=0;i<aBCFlag.size();++i){
    if( aBCFlag[i] == 0 ){ continue; }
    y[i] += weight_bc*vec[i];
  }
}

// for preconditioner
void delfem2::CDef_LaplacianLinear::SolvePrecond(double* v) const
{
  /*
  const unsigned int np = aBCFlag.size()/3;
  for(unsigned int ip=0;ip<np;++ip){
    double tmp[3];
    MatVec3(tmp, aDiaInv.data()+ip*9, v+ip*3);
    v[ip*3+0] = tmp[0];
    v[ip*3+1] = tmp[1];
    v[ip*3+2] = tmp[2];
  }
   */
}





// ============================================

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
  for(unsigned int i=0;i<np*3;++i){ y[i] *= beta; }
  for(unsigned int ip=0;ip<np;++ip){
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
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
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
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
  for(unsigned int i=0;i<aBCFlag.size();++i){
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
  for(unsigned int i=0;i<np*3;++i){ aXYZ1[i] += aUpd[i]; }
}

// ======================================================

void delfem2::CDef_ArapEdge::Init
 (const std::vector<double>& aXYZ0,
  const std::vector<unsigned int>& aTri,
  double weight_bc0,
  const std::vector<int>& aBCFlag0,
  bool is_preconditioner0)
{
  this->weight_bc = weight_bc0;
  this->is_preconditioner = is_preconditioner0;
  this->aBCFlag = aBCFlag0;
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
  for(unsigned int i=0;i<np*6;++i){ y[i] *= beta; }
  for(unsigned int ip=0;ip<np;++ip){
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
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
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
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
  for(unsigned int i=0;i<aBCFlag.size();++i){
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

void delfem2::CDef_ArapEdge::SolvePrecond(double* v) const
{
  const unsigned int np = psup_ind.size()-1;
  for(unsigned int ip=0;ip<np*2;++ip){
    double tmp[3];
    MatVec3(tmp, aDiaInv.data()+ip*9, v+ip*3);
    v[ip*3+0] = tmp[0];
    v[ip*3+1] = tmp[1];
    v[ip*3+2] = tmp[2];
  }
}


void delfem2::CDef_ArapEdge::Deform
 (std::vector<double>& aXYZ1,
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
  for(unsigned int ip=0;ip<np;++ip){
    aXYZ1[ip*3+0] += aUpd[ip*3+0];
    aXYZ1[ip*3+1] += aUpd[ip*3+1];
    aXYZ1[ip*3+2] += aUpd[ip*3+2];
    double q0[4]; Quat_CartesianAngle(q0, aUpd.data()+np*3+ip*3);
    double q1[4]; QuatQuat(q1, q0, aQuat.data()+ip*4);
    Copy_Quat(aQuat.data()+ip*4, q1);
  }
}


// ===========================================================
// below: implementation of CDef_Arap class

void delfem2::CDef_Arap::Init
 (const std::vector<double>& aXYZ0,
  const std::vector<unsigned int>& aTri,
  bool is_preconditioner_)
{
  this->is_preconditioner = is_preconditioner_;
  const unsigned int np = aXYZ0.size()/3;
  JArray_PSuP_MeshElem(psup_ind, psup,
                       aTri.data(), aTri.size()/3, 3,
                       aXYZ0.size()/3);
  JArray_Sort(psup_ind, psup);
  {
    std::vector<unsigned int> psup_ind1, psup1;
    JArray_Extend(psup_ind1, psup1,
                  psup_ind, psup);
    JArray_Sort(psup_ind1, psup1);
    Mat.Initialize(np, 3, true);
    assert( psup_ind1.size() == np+1 );
    Mat.SetPattern(psup_ind1.data(), psup_ind1.size(),
                   psup1.data(), psup1.size());
  }
  
  Precomp.resize(np*9);
  for(unsigned int ip=0;ip<np;++ip){
    const CVec3d Pi(aXYZ0.data()+ip*3);
    CMat3d LM; LM.SetZero();
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const unsigned int jp = psup[ipsup];
      const CVec3d v0 = (CVec3d(aXYZ0.data()+jp*3)-Pi);
      LM += Mat3_CrossCross(v0);
    }
    CMat3d LMi = LM.Inverse();
    LMi.CopyTo(Precomp.data()+ip*9);
  }
  
  this->Prec.Clear();
  if( is_preconditioner ){
    this->Prec.Initialize_ILU0(Mat);
  }
  
}


void delfem2::CDef_Arap::Deform
(std::vector<double>& aXYZ1,
 std::vector<double>& aQuat1,
 const std::vector<double>& aXYZ0,
 const std::vector<int>& aBCFlag)
{
  const unsigned int np = aXYZ0.size()/3;
  Mat.SetZero();
  this->aRes1.assign(np*3, 0.0);
  std::vector<int> tmp_buffer;
  for(unsigned int ip=0;ip<np;++ip){
    std::vector<unsigned int> aIP;
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      aIP.push_back(psup[ipsup]);
    }
    aIP.push_back(ip);
    std::vector<double> eM, eR;
    def::dWddW_ArapEnergy(eM,eR,
                          Precomp.data()+ip*9,
                          aIP,aXYZ0,aXYZ1,aQuat1);
    Mat.Mearge(aIP.size(), aIP.data(),
               aIP.size(), aIP.data(),
               9, eM.data(),
               tmp_buffer);
    for(unsigned int iip=0;iip<aIP.size();++iip){
      const int jp0 = aIP[iip];
      aRes1[jp0*3+0] += eR[iip*3+0];
      aRes1[jp0*3+1] += eR[iip*3+1];
      aRes1[jp0*3+2] += eR[iip*3+2];
    }
  }
  Mat.AddDia(1.0e-8);
  
  Mat.SetFixedBC(aBCFlag.data());
  setRHS_Zero(aRes1, aBCFlag, 0);
  
  aUpd1.resize(aRes1.size());
  std::vector<double> aConvHist;
  
  if( is_preconditioner ){
    this->Prec.SetValueILU(Mat);
    this->Prec.DoILUDecomp();
    aConvHist = Solve_PCG(aRes1.data(), aUpd1.data(),
                         aRes1.size(), 1.0e-7, 300, Mat, Prec);
  }
  else{
    aConvHist = Solve_CG(aRes1.data(), aUpd1.data(),
                         aRes1.size(), 1.0e-7, 300, Mat);
  }
  
  std::cout << aConvHist.size() << std::endl;
  for(unsigned int i=0;i<np*3;++i){ aXYZ1[i] -= aUpd1[i]; }
  for(int itr=0;itr<1;++itr){
    UpdateRotationsByMatchingCluster(aQuat1,
                                     aXYZ0,aXYZ1,psup_ind,psup);
  }
}
