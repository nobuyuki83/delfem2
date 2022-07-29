/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/defarap.h"

#include <cstring>  // memcpy
#include <utility>

#include "delfem2/vec3.h"
#include "delfem2/quat.h"
#include "delfem2/mat3.h"
#include "delfem2/mat3_funcs.h"
#include "delfem2/svd3.h"
#include "delfem2/msh_topology_uniform.h"
#include "delfem2/view_vectorx.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/jagarray.h"

namespace delfem2::deflap {

DFM2_INLINE void dWddW_ArapEnergy(
    std::vector<double> &eM,
    std::vector<double> &eR,
    const double *Minv,
    const std::vector<unsigned int> &aIP,
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aXYZ1,
    const std::vector<double> &aQuat1) {
  const size_t nIP = aIP.size();
  const size_t nNg = nIP - 1; // number of neighbor
  unsigned int ip = aIP[nNg];
  const CVec3d Pi(aXYZ0.data() + ip * 3);
  const CMat3d LMi(Minv);
  const CMat3d Mrot = CMat3d::Quat(aQuat1.data() + ip * 4);
  eM.assign(nIP * nIP * 9, 0.0);
  for (unsigned int jjp = 0; jjp < nNg; ++jjp) {
    for (unsigned int kkp = 0; kkp < nNg; ++kkp) { // nonlinear component
      const CVec3d vj = (CVec3d(aXYZ0.data() + aIP[jjp] * 3) - Pi);
      const CVec3d vk = (CVec3d(aXYZ0.data() + aIP[kkp] * 3) - Pi);
      CMat3d L1 = Mrot * CMat3d::Skew(vk) * LMi * CMat3d::Skew(vj) * Mrot.transpose();
      L1.AddToScale(eM.data() + (kkp * nIP + jjp) * 9, -1.0);
      L1.AddToScale(eM.data() + (nNg * nIP + nNg) * 9, -1.0);
      L1.AddToScale(eM.data() + (nNg * nIP + jjp) * 9, +1.0);
      L1.AddToScale(eM.data() + (kkp * nIP + nNg) * 9, +1.0);
    }
    { // linear component
      CMat3d L1 = CMat3d::Identity();
      L1.AddToScale(eM.data() + (jjp * nIP + jjp) * 9, +1.0);
      L1.AddToScale(eM.data() + (nNg * nIP + nNg) * 9, +1.0);
      L1.AddToScale(eM.data() + (nNg * nIP + jjp) * 9, -1.0);
      L1.AddToScale(eM.data() + (jjp * nIP + nNg) * 9, -1.0);
    }
  }

  //
  eR.assign(nIP * 3, 0.0);
  const CVec3d pi(aXYZ1.data() + ip * 3);
  CMat3d LM;
  LM.setZero();
  for (unsigned int jjp = 0; jjp < nNg; ++jjp) {
    const unsigned int jp = aIP[jjp];
    const CVec3d v0 = Mrot * (CVec3d(aXYZ0.data() + jp * 3) - Pi);
    CVec3d pj(aXYZ1.data() + jp * 3);
    const CVec3d v1 = pj - pi;
    const CVec3d r = -(v1 - v0);
    r.AddToScale(eR.data() + nNg * 3, +1.);
    r.AddToScale(eR.data() + jjp * 3, -1.);
  }
}

}

// ============================================

delfem2::CDef_ArapEdgeLinearDisponly::CDef_ArapEdgeLinearDisponly(
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aTri,
    double weight_bc0,
    std::vector<int> aBCFlag0) :
    weight_bc(weight_bc0),
    aBCFlag(std::move(aBCFlag0)) {
  const size_t np = aXYZ0.size() / 3;
  JArray_PSuP_MeshElem(
      psup_ind, psup,
      aTri.data(), aTri.size() / 3, 3,
      aXYZ0.size() / 3);
  JArray_Sort(psup_ind, psup);
  // ------
  const size_t ne = psup.size();
  // -----
  aMatEdge.resize(ne * 9 * 2);
  assert(psup_ind.size() == np + 1);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      //          unsigned int jp = psup[ipsup];
      Mat3_Identity(aMatEdge.data() + ipsup * 18, +1.0);
      Mat3_Identity(aMatEdge.data() + ipsup * 18 + 9, -1.0);
    }
  }
  // ---
  vec_tmp.resize(ne * 3);
}

void delfem2::CDef_ArapEdgeLinearDisponly::JacobiTVecTmp(
    double *y,
    double alpha, double beta) const {
  const size_t np = aBCFlag.size() / 3;
  for (unsigned int i = 0; i < np * 3; ++i) { y[i] *= beta; }
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      unsigned int jp0 = psup[ipsup];
      MatTVec3_ScaleAdd(y + ip * 3,
                        aMatEdge.data() + ipsup * 18,
                        vec_tmp.data() + ipsup * 3,
                        alpha, 1.0);
      MatTVec3_ScaleAdd(y + jp0 * 3,
                        aMatEdge.data() + ipsup * 18 + 9,
                        vec_tmp.data() + ipsup * 3,
                        alpha, 1.0);
    }
  }
}

void delfem2::CDef_ArapEdgeLinearDisponly::MakeLinearSystem(
    double *aRhs,
    const double *aXYZ0,
    const double *aXYZ1) const {
  const size_t np = aBCFlag.size() / 3;
  const size_t ne = psup.size();
  vec_tmp.assign(ne * 3, 0);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      const unsigned int jp0 = psup[ipsup];
      const double d0[3] = {aXYZ0[jp0 * 3 + 0] - aXYZ0[ip * 3 + 0],
                            aXYZ0[jp0 * 3 + 1] - aXYZ0[ip * 3 + 1],
                            aXYZ0[jp0 * 3 + 2] - aXYZ0[ip * 3 + 2]};
      const double d1[3] = {aXYZ1[jp0 * 3 + 0] - aXYZ1[ip * 3 + 0],
                            aXYZ1[jp0 * 3 + 1] - aXYZ1[ip * 3 + 1],
                            aXYZ1[jp0 * 3 + 2] - aXYZ1[ip * 3 + 2]};
      vec_tmp[ipsup * 3 + 0] += +(d0[0] - d1[0]);
      vec_tmp[ipsup * 3 + 1] += +(d0[1] - d1[1]);
      vec_tmp[ipsup * 3 + 2] += +(d0[2] - d1[2]);
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

void delfem2::CDef_ArapEdgeLinearDisponly::MatVec(
    double *y,
    double alpha,
    const double *vec,
    double beta) const {
  const size_t np = aBCFlag.size() / 3;
  std::fill(vec_tmp.begin(), vec_tmp.end(), 0.0);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      unsigned int jp0 = psup[ipsup];
      MatVec3_ScaleAdd(vec_tmp.data() + ipsup * 3,
                       aMatEdge.data() + ipsup * 18,
                       vec + ip * 3,
                       1.0, 1.0);
      MatVec3_ScaleAdd(vec_tmp.data() + ipsup * 3,
                       aMatEdge.data() + ipsup * 18 + 9,
                       vec + jp0 * 3,
                       1.0, 1.0);
    }
  }
  this->JacobiTVecTmp(y,
                      alpha, beta);
  // add diagonal for fixed boundary condition
  for (unsigned int i = 0; i < aBCFlag.size(); ++i) {
    if (aBCFlag[i] == 0) { continue; }
    y[i] += weight_bc * vec[i];
  }
}

void delfem2::CDef_ArapEdgeLinearDisponly::Deform(
    std::vector<double> &aXYZ1,
    const std::vector<double> &aXYZ0) {
  const size_t np = aBCFlag.size() / 3;
  std::vector<double> aRhs(np * 3, 0.0);
  this->MakeLinearSystem(aRhs.data(),
                         aXYZ0.data(), aXYZ1.data());
  std::vector<double> aUpd(np * 3, 0.0);
  {
    const std::size_t n = np * 3;
    std::vector<double> tmp0(n), tmp1(n);
    std::vector<double> aRes = Solve_CG(
        ViewAsVectorXd(aRhs),
        ViewAsVectorXd(aUpd),
        ViewAsVectorXd(tmp0),
        ViewAsVectorXd(tmp1),
        1.0e-4, 300, *this);
  }
//  std::cout << "iframe: " << iframe << "   nitr:" << aRes.size() << std::endl;
  for (unsigned int i = 0; i < np * 3; ++i) { aXYZ1[i] += aUpd[i]; }
}

// ======================================================

void delfem2::CDef_ArapEdge::Init(
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aTri,
    double weight_bc0,
    const std::vector<int> &aBCFlag0,
    bool is_preconditioner0) {
  this->weight_bc = weight_bc0;
  this->is_preconditioner = is_preconditioner0;
  this->aBCFlag = aBCFlag0;
  const size_t np = aXYZ0.size() / 3;
  JArray_PSuP_MeshElem(
      psup_ind, psup,
      aTri.data(), aTri.size() / 3, 3,
      aXYZ0.size() / 3);
  JArray_Sort(psup_ind, psup);
  // ---------
  const size_t ne = psup.size();
  // -----
  aMatEdge.resize(ne * 27);
  assert(psup_ind.size() == np + 1);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      //          unsigned int jp = psup[ipsup];
      Mat3_Identity(aMatEdge.data() + ipsup * 27 + 0, +1.0);
      Mat3_Identity(aMatEdge.data() + ipsup * 27 + 9, -1.0);
    }
  }
  // ---
  vec_tmp.resize(ne * 3);
}

void delfem2::CDef_ArapEdge::JacobiTVecTmp(
    double *y,
    double alpha,
    double beta) const {
  assert(!psup_ind.empty());
  const size_t np = psup_ind.size() - 1;
  for (unsigned int i = 0; i < np * 6; ++i) { y[i] *= beta; }
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      unsigned int jp0 = psup[ipsup];
      MatTVec3_ScaleAdd(y + ip * 3,
                        aMatEdge.data() + ipsup * 27 + 0,
                        vec_tmp.data() + ipsup * 3,
                        alpha, 1.0);
      MatTVec3_ScaleAdd(y + jp0 * 3,
                        aMatEdge.data() + ipsup * 27 + 9,
                        vec_tmp.data() + ipsup * 3,
                        alpha, 1.0);
      MatTVec3_ScaleAdd(y + np * 3 + ip * 3,
                        aMatEdge.data() + ipsup * 27 + 18,
                        vec_tmp.data() + ipsup * 3,
                        alpha, 1.0);
    }
  }
}

void delfem2::CDef_ArapEdge::MatVec(
    double *y,
    double alpha,
    const double *vec,
    double beta) const {
  const size_t np = psup_ind.size() - 1;
  std::fill(vec_tmp.begin(), vec_tmp.end(), 0.0);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      unsigned int jp0 = psup[ipsup];
      MatVec3_ScaleAdd(vec_tmp.data() + ipsup * 3,
                       aMatEdge.data() + ipsup * 27 + 0,
                       vec + ip * 3,
                       1.0, 1.0);
      MatVec3_ScaleAdd(vec_tmp.data() + ipsup * 3,
                       aMatEdge.data() + ipsup * 27 + 9,
                       vec + jp0 * 3,
                       1.0, 1.0);
      MatVec3_ScaleAdd(vec_tmp.data() + ipsup * 3,
                       aMatEdge.data() + ipsup * 27 + 18,
                       vec + (np + ip) * 3,
                       1.0, 1.0);
    }
  }
  this->JacobiTVecTmp(y,
                      alpha, beta);
  // add diagonal for fixed boundary condition
  for (unsigned int i = 0; i < aBCFlag.size(); ++i) {
    if (aBCFlag[i] == 0) { continue; }
    y[i] += weight_bc * vec[i];
    //      y[np*3+i] += weight_bc*vec[np*3+i];
  }
}

void delfem2::CDef_ArapEdge::MakeLinearSystem(
    double *aRhs,
    const double *aXYZ0,
    const double *aXYZ1,
    const double *aQuat) {
  const size_t np = psup_ind.size() - 1;
  const size_t ne = psup.size();
  vec_tmp.assign(ne * 3, 0);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      const unsigned int jp0 = psup[ipsup];
      const double *q0 = aQuat + ip * 4;
      const double d0[3] = {aXYZ0[jp0 * 3 + 0] - aXYZ0[ip * 3 + 0],
                            aXYZ0[jp0 * 3 + 1] - aXYZ0[ip * 3 + 1],
                            aXYZ0[jp0 * 3 + 2] - aXYZ0[ip * 3 + 2]};
      const double d1[3] = {aXYZ1[jp0 * 3 + 0] - aXYZ1[ip * 3 + 0],
                            aXYZ1[jp0 * 3 + 1] - aXYZ1[ip * 3 + 1],
                            aXYZ1[jp0 * 3 + 2] - aXYZ1[ip * 3 + 2]};
      double Rd0[3];
      QuatVec(Rd0, q0, d0);
      vec_tmp[ipsup * 3 + 0] += +(Rd0[0] - d1[0]);
      vec_tmp[ipsup * 3 + 1] += +(Rd0[1] - d1[1]);
      vec_tmp[ipsup * 3 + 2] += +(Rd0[2] - d1[2]);
      Mat3_Spin_ScaleAdd(aMatEdge.data() + ipsup * 27 + 18,
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

void delfem2::CDef_ArapEdge::MakePreconditionerJacobi() {
  assert(!psup_ind.empty());
  const size_t np = psup_ind.size() - 1;
  aDiaInv.assign(np * 2 * 9, 0.0);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      const unsigned int jp0 = psup[ipsup];
      MatTMat3_ScaleAdd(aDiaInv.data() + ip * 9,
                        aMatEdge.data() + ipsup * 27 + 0,
                        aMatEdge.data() + ipsup * 27 + 0,
                        1.0, 1.0);
      MatTMat3_ScaleAdd(aDiaInv.data() + jp0 * 9,
                        aMatEdge.data() + ipsup * 27 + 9,
                        aMatEdge.data() + ipsup * 27 + 9,
                        1.0, 1.0);
      MatTMat3_ScaleAdd(aDiaInv.data() + (np + ip) * 9,
                        aMatEdge.data() + ipsup * 27 + 18,
                        aMatEdge.data() + ipsup * 27 + 18,
                        1.0, 1.0);
    }
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (int idim = 0; idim < 3; ++idim) {
      if (aBCFlag[ip * 3 + idim] == 0) { continue; }
      aDiaInv[ip * 9 + idim * 3 + idim] += weight_bc;
    }
  }
  for (unsigned int ip = 0; ip < np * 2; ++ip) {
    Inverse_Mat3(aDiaInv.data() + ip * 9);
  }
}

void delfem2::CDef_ArapEdge::SolvePrecond(double *v) const {
  assert(!psup_ind.empty());
  const size_t np = psup_ind.size() - 1;
  for (unsigned int ip = 0; ip < np * 2; ++ip) {
    double tmp[3];
    MatVec3(tmp, aDiaInv.data() + ip * 9, v + ip * 3);
    v[ip * 3 + 0] = tmp[0];
    v[ip * 3 + 1] = tmp[1];
    v[ip * 3 + 2] = tmp[2];
  }
}

void delfem2::CDef_ArapEdge::Deform(
    std::vector<double> &aXYZ1,
    std::vector<double> &aQuat,
    const std::vector<double> &aXYZ0) {
  const size_t np = psup_ind.size() - 1;
  std::vector<double> aRhs(np * 6, 0.0);
  this->MakeLinearSystem(
      aRhs.data(),
      aXYZ0.data(), aXYZ1.data(), aQuat.data());
  std::vector<double> aUpd(np * 6, 0.0);
  std::vector<double> aConvHist;
  if (is_preconditioner) {
    const std::size_t n = np * 6;
    std::vector<double> tmp0(n), tmp1(n);
    this->MakePreconditionerJacobi();
    aConvHist = Solve_PCG(
        ViewAsVectorXd(aRhs),
        ViewAsVectorXd(aUpd),
        ViewAsVectorXd(tmp0),
        ViewAsVectorXd(tmp1),
        1.0e-4, 400, *this, *this);
  } else {
    const std::size_t n = np * 6;
    std::vector<double> tmp0(n), tmp1(n);
    aConvHist = Solve_CG(
        ViewAsVectorXd(aRhs),
        ViewAsVectorXd(aUpd),
        ViewAsVectorXd(tmp0),
        ViewAsVectorXd(tmp1),
        1.0e-4, 400, *this);
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    aXYZ1[ip * 3 + 0] += aUpd[ip * 3 + 0];
    aXYZ1[ip * 3 + 1] += aUpd[ip * 3 + 1];
    aXYZ1[ip * 3 + 2] += aUpd[ip * 3 + 2];
    double q0[4];
    Quat_CartesianAngle(q0, aUpd.data() + np * 3 + ip * 3);
    double q1[4];
    QuatQuat(q1, q0, aQuat.data() + ip * 4);
    Copy_Quat(aQuat.data() + ip * 4, q1);
  }
}


// ===========================================================
// below: implementation of CDef_Arap class

void delfem2::CDef_Arap::Init(
    const std::vector<double> &aXYZ0,
    const std::vector<unsigned int> &aTri,
    bool is_preconditioner_) {
  this->is_preconditioner = is_preconditioner_;
  const size_t np = aXYZ0.size() / 3;
  JArray_PSuP_MeshElem(
      psup_ind, psup,
      aTri.data(), aTri.size() / 3, 3,
      aXYZ0.size() / 3);
  JArray_Sort(psup_ind, psup);
  {
    std::vector<unsigned int> psup_ind1, psup1;
    JArray_Extend(
        psup_ind1, psup1,
        psup_ind.data(), psup_ind.size(), psup.data());
    JArray_Sort(psup_ind1, psup1);
    Mat.Initialize(
        static_cast<unsigned int>(np), 3, true);
    assert(psup_ind1.size() == np + 1);
    Mat.SetPattern(psup_ind1.data(), psup_ind1.size(),
                   psup1.data(), psup1.size());
  }

  Precomp.resize(np * 9);
  for (unsigned int ip = 0; ip < np; ++ip) {
    const CVec3d Pi(aXYZ0.data() + ip * 3);
    CMat3d LM;
    LM.setZero();
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      const unsigned int jp = psup[ipsup];
      const CVec3d v0 = (CVec3d(aXYZ0.data() + jp * 3) - Pi);
      LM += Mat3_CrossCross(v0);
    }
    CMat3d LMi = LM.Inverse();
    LMi.CopyTo(Precomp.data() + ip * 9);
  }

  this->Prec.Clear();
  if (is_preconditioner) {
    this->Prec.SetPattern0(Mat);
  }

}

void delfem2::CDef_Arap::Deform(
    std::vector<double> &aXYZ1,
    std::vector<double> &aQuat1,
    const std::vector<double> &aXYZ0,
    const std::vector<int> &aBCFlag) {
  const size_t np = aXYZ0.size() / 3;
  Mat.setZero();
  this->aRes1.assign(np * 3, 0.0);
  std::vector<unsigned int> tmp_buffer;
  for (unsigned int ip = 0; ip < np; ++ip) {
    std::vector<unsigned int> aIP;
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      aIP.push_back(psup[ipsup]);
    }
    aIP.push_back(ip);
    std::vector<double> eM, eR;
    deflap::dWddW_ArapEnergy(
        eM, eR,
        Precomp.data() + ip * 9,
        aIP, aXYZ0, aXYZ1, aQuat1);
    Mearge(Mat,
           aIP.size(), aIP.data(),
           aIP.size(), aIP.data(),
           9, eM.data(),
           tmp_buffer);
    for (unsigned int iip = 0; iip < aIP.size(); ++iip) {
      const unsigned int jp0 = aIP[iip];
      aRes1[jp0 * 3 + 0] += eR[iip * 3 + 0];
      aRes1[jp0 * 3 + 1] += eR[iip * 3 + 1];
      aRes1[jp0 * 3 + 2] += eR[iip * 3 + 2];
    }
  }
  Mat.AddDia(1.0e-8);

  // for(unsigned int icrs=0;icrs<Mat.val_crs_.size()/9;++icrs){
  {
    unsigned int icrs = 1352;
    std::cout << icrs << std::endl;
    std::cout << "   " << Mat.val_crs_[icrs*9+0] << " " << Mat.val_crs_[icrs*9+1] << " " << Mat.val_crs_[icrs*9+2] << std::endl;
    std::cout << "   " << Mat.val_crs_[icrs*9+3] << " " << Mat.val_crs_[icrs*9+4] << " " << Mat.val_crs_[icrs*9+5] << std::endl;
    std::cout << "   " << Mat.val_crs_[icrs*9+6] << " " << Mat.val_crs_[icrs*9+7] << " " << Mat.val_crs_[icrs*9+8] << std::endl;
  }

  Mat.SetFixedBC(aBCFlag.data());
  setRHS_Zero(aRes1, aBCFlag, 0);

  aUpd1.resize(aRes1.size());
  if (is_preconditioner) {
    this->Prec.CopyValue(Mat);
    this->Prec.Decompose();
    const std::size_t n = np * 3;
    std::vector<double> tmp0(n), tmp1(n);
    aConvHist = Solve_PCG(
        ViewAsVectorXd(aRes1),
        ViewAsVectorXd(aUpd1),
        ViewAsVectorXd(tmp0),
        ViewAsVectorXd(tmp1),
        1.0e-7, 300, Mat, Prec);
  } else {
    const std::size_t n = np * 3;
    assert(aRes1.size() == n && aUpd1.size() == n);
    std::vector<double> tmp0(n), tmp1(n);
    aConvHist = Solve_CG(
        ViewAsVectorXd(aRes1),
        ViewAsVectorXd(aUpd1),
        ViewAsVectorXd(tmp0),
        ViewAsVectorXd(tmp1),
        1.0e-7, 300, Mat);
  }

  for (unsigned int i = 0; i < np * 3; ++i) { aXYZ1[i] -= aUpd1[i]; }
  // ----
  /*
  for(int itr=0;itr<1;++itr){
    UpdateRotationsByMatchingCluster_Linear(aQuat1,
                                            aXYZ0,aXYZ1,psup_ind,psup);
  }
   */

}

void delfem2::CDef_Arap::UpdateQuaternions_Svd(
    std::vector<double> &aQuat1,
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aXYZ1) const {
  std::size_t np = aXYZ1.size() / 3;
  for (unsigned int ip = 0; ip < np; ++ip) {
    UpdateRotationsByMatchingCluster_SVD(
        aQuat1,
        ip, aXYZ0, aXYZ1, psup_ind, psup);
  }
}

DFM2_INLINE void delfem2::UpdateRotationsByMatchingCluster_Linear(
    std::vector<double> &aQuat1,
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aXYZ1,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup) {
  const size_t np = aXYZ0.size() / 3;
  for (unsigned int ip = 0; ip < np; ++ip) {
    const CVec3d Pi(aXYZ0.data() + ip * 3);
    const CVec3d pi(aXYZ1.data() + ip * 3);
    const CQuatd Qi(aQuat1.data() + ip * 4);
    CMat3d Mat;
    Mat.setZero();
    CVec3d rhs;
    rhs.setZero();
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      const unsigned int jp = psup[ipsup];
      const CVec3d v0 = Qi.RotateVector(CVec3d(aXYZ0.data() + jp * 3) - Pi);
      const CVec3d d01 = CVec3d(aXYZ1.data() + jp * 3) - pi - v0;
      Mat += Mat3_CrossCross(v0);
      rhs += d01.cross(v0);
    }
    CVec3d sol = Mat.Inverse() * rhs;
    CQuatd q0 = Quat_CartesianAngle(sol);
    CQuatd q1 = q0 * Qi;
    q1.CopyTo(aQuat1.data() + ip * 4);
  }
}

DFM2_INLINE void delfem2::UpdateRotationsByMatchingCluster_SVD(
    std::vector<double> &aQuat1,
    unsigned int ip,
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aXYZ1,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup) {
  const CVec3d Pi(aXYZ0.data() + ip * 3);
  const CVec3d pi(aXYZ1.data() + ip * 3);
  const CMat3d R0i = CMat3d::Quat(aQuat1.data() + ip * 4);
  double A[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
    const unsigned int jp = psup[ipsup];
    const CVec3d &v0 = R0i * (CVec3d(aXYZ0.data() + jp * 3) - Pi);
    const CVec3d v1 = CVec3d(aXYZ1.data() + jp * 3) - pi;
    const double *dp = v1.p;
    const double *dq = v0.p;
    A[0 * 3 + 0] += dp[0] * dq[0];
    A[0 * 3 + 1] += dp[0] * dq[1];
    A[0 * 3 + 2] += dp[0] * dq[2];
    A[1 * 3 + 0] += dp[1] * dq[0];
    A[1 * 3 + 1] += dp[1] * dq[1];
    A[1 * 3 + 2] += dp[1] * dq[2];
    A[2 * 3 + 0] += dp[2] * dq[0];
    A[2 * 3 + 1] += dp[2] * dq[1];
    A[2 * 3 + 2] += dp[2] * dq[2];
  }
  CMat3d dRi;
  GetRotPolarDecomp(dRi.p_, A, 40);
  const CMat3d R1 = dRi * R0i;
  const CQuatd q1 = R1.GetQuaternion();
  q1.CopyTo(aQuat1.data() + ip * 4);
}
