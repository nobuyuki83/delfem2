/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/rig_geo3.h"

#include <cassert>
#include <climits>

#ifndef M_PI
#  define M_PI 3.1415926535
#endif

// ------------------------------------------------------------

/*
namespace delfem2::rig_v3q {

DFM2_INLINE bool isActive_AABB(const double aabb[6]){
    return aabb[0] <= aabb[1];
}

DFM2_INLINE void copy_AABB(double aabb[6], const double aabb0[6]){
  if( aabb == aabb0 ) return;
  for(int i=0;i<6;++i){ aabb[i] = aabb0[i]; }
}

DFM2_INLINE void myAdd_AABB(double aabb[6], const double aabb0[6], const double aabb1[6])
{
  if( !isActive_AABB(aabb0) && !isActive_AABB(aabb1) ){ aabb[0]=1; aabb[1]=-1; return; }
  if( !isActive_AABB(aabb0) ){ copy_AABB(aabb, aabb1); return; }
  if( !isActive_AABB(aabb1) ){ copy_AABB(aabb, aabb0); return; }
  aabb[0] = ( aabb0[0] < aabb1[0] ) ? aabb0[0] : aabb1[0];
  aabb[1] = ( aabb0[1] > aabb1[1] ) ? aabb0[1] : aabb1[1];
  aabb[2] = ( aabb0[2] < aabb1[2] ) ? aabb0[2] : aabb1[2];
  aabb[3] = ( aabb0[3] > aabb1[3] ) ? aabb0[3] : aabb1[3];
  aabb[4] = ( aabb0[4] < aabb1[4] ) ? aabb0[4] : aabb1[4];
  aabb[5] = ( aabb0[5] > aabb1[5] ) ? aabb0[5] : aabb1[5];
}

DFM2_INLINE double MyDotX(
    const double* va,
    const double* vb,
    unsigned int n)
{
  double r = 0.0;
  for(unsigned int i=0;i<n;i++){ r += va[i]*vb[i]; }
  return r;
}


DFM2_INLINE void MyMatMatTX(
    double* M, // [ni, nj]
    unsigned int ni,
    unsigned int nj,
    const double* A, // [ni, nk]
    unsigned int nk,
    const double* B) // [nj, nk]
{
  for(unsigned int i=0;i<ni;++i){
    for(unsigned int j=0;j<nj;++j){
      M[i*nj+j] = 0.0;
      for(unsigned int k=0;k<nk;++k){
        M[i*nj+j] += A[i*nk+k] * B[j*nk+k];
      }
    }
  }
}

}
*/

void delfem2::SparsifyMatrixRow(
    std::vector<double> &aWBone_RigSparse,
    std::vector<unsigned int> &aIdBone_RigSparse,
    const double *aW,
    unsigned int nrow,
    unsigned int ncol,
    double thres) {
  unsigned int nbone_nonzeroweight = 0;
  for (unsigned int ip = 0; ip < nrow; ++ip) {
    unsigned icnt = 0;
    for (unsigned int ib = 0; ib < ncol; ++ib) {
      if (aW[ip * ncol + ib] < thres) { continue; }
      icnt++;
    }
    if (icnt > nbone_nonzeroweight) { nbone_nonzeroweight = icnt; }
  }
  aWBone_RigSparse.resize(nbone_nonzeroweight * nrow);
  aIdBone_RigSparse.resize(nbone_nonzeroweight * nrow);
  for (unsigned int ip = 0; ip < nrow; ++ip) {
    unsigned icnt = 0;
    double w_sum = 0.0;
    for (unsigned int ib = 0; ib < ncol; ++ib) {
      if (aW[ip * ncol + ib] < thres) { continue; }
      w_sum += aW[ip * ncol + ib];
      aWBone_RigSparse[ip * nbone_nonzeroweight + icnt] = aW[ip * ncol + ib];
      aIdBone_RigSparse[ip * nbone_nonzeroweight + icnt] = ib;
      icnt++;
    }
    if (icnt > nbone_nonzeroweight) { nbone_nonzeroweight = icnt; }
  }
}

DFM2_INLINE void delfem2::Transpose_Mat(
    std::vector<double> &At,
    const std::vector<double> &A,
    unsigned int nrow,
    unsigned int ncol) {
  At.resize(A.size());
  for (unsigned int i = 0; i < nrow; ++i) {
    for (unsigned int j = 0; j < ncol; ++j) {
      At[j * nrow + i] = A[i * ncol + j];
    }
  }
}

DFM2_INLINE void delfem2::Points3_WeighttranspPosition(
    std::vector<double> &aPos,
    const std::vector<double> &Weighttransp,
    const std::vector<double> &aXYZ0) {
  const size_t nXYZ = aXYZ0.size() / 3;
  const size_t nPos = Weighttransp.size() / nXYZ;
  aPos.assign(nPos * 3, 0.0);
  for (unsigned int ib = 0; ib < nPos; ++ib) {
    aPos[ib * 3 + 0] = 0;
    aPos[ib * 3 + 1] = 0;
    aPos[ib * 3 + 2] = 0;
    for (unsigned int ip = 0; ip < nXYZ; ++ip) {
      aPos[ib * 3 + 0] += Weighttransp[ip * nPos + ib] * aXYZ0[ip * 3 + 0];
      aPos[ib * 3 + 1] += Weighttransp[ip * nPos + ib] * aXYZ0[ip * 3 + 1];
      aPos[ib * 3 + 2] += Weighttransp[ip * nPos + ib] * aXYZ0[ip * 3 + 2];
    }
  }
}

DFM2_INLINE void delfem2::Points3_WeightsparsePosition(
    std::vector<double> &aPos0,
    unsigned int nPos,
    const std::vector<double> &aSparseW,
    const std::vector<unsigned int> &aSparseIdp,
    const std::vector<double> &aXYZ0) {
  assert(aSparseW.size() == aSparseIdp.size());
  const size_t np_sparse = aSparseIdp.size() / nPos;
  assert(aSparseIdp.size() % np_sparse == 0);
  aPos0.assign(nPos * 3, 0.0);
  for (unsigned int ib = 0; ib < nPos; ++ib) {
    aPos0[ib * 3 + 0] = 0;
    aPos0[ib * 3 + 1] = 0;
    aPos0[ib * 3 + 2] = 0;
    for (unsigned int ipb = 0; ipb < np_sparse; ++ipb) {
      const unsigned int ip = aSparseIdp[ib * np_sparse + ipb];
      const double w0 = aSparseW[ib * np_sparse + ipb];
      aPos0[ib * 3 + 0] += w0 * aXYZ0[ip * 3 + 0];
      aPos0[ib * 3 + 1] += w0 * aXYZ0[ip * 3 + 1];
      aPos0[ib * 3 + 2] += w0 * aXYZ0[ip * 3 + 2];
    }
  }
}

// above: matrix operation
// ------------------------------------------------------------
// below: skinning

DFM2_INLINE void delfem2::CRigBone::SetRotationBryant(
    double rx, double ry, double rz) {
  Quat_Bryant(quatRelativeRot, rx, ry, rz);
}

DFM2_INLINE void delfem2::CRigBone::DeformSkin(
    double pos2[3],
    const double pos0[3]) const {
  const double pos0a[4] = {pos0[0], pos0[1], pos0[2], 1.0};
  double pos1a[4];
  MatVec4(pos1a, invBindMat, pos0a);
  double pos2a[4];
  MatVec4(pos2a, affmat3Global, pos1a);
  pos2[0] = pos2a[0];
  pos2[1] = pos2a[1];
  pos2[2] = pos2a[2];
}

DFM2_INLINE void delfem2::CRigBone::SetTranslation(
    double tx, double ty, double tz) {
  this->transRelative[0] = tx;
  this->transRelative[1] = ty;
  this->transRelative[2] = tz;
}

DFM2_INLINE void delfem2::UpdateBoneRotTrans(
    std::vector<CRigBone> &aBone) {
  for (std::size_t ibone = 0; ibone < aBone.size(); ++ibone) {
    CMat4d m01 = CMat4d::Translation(aBone[ibone].transRelative);
    m01 = m01 * CMat4d::Quat(aBone[ibone].quatRelativeRot);
    m01 = m01 * CMat4d::AffineScale(aBone[ibone].scale);
    // m01 = T*Q*S
    const int ibone_p = aBone[ibone].ibone_parent;
    if (ibone_p < 0 || ibone_p >= (int) aBone.size()) { // root bone
      Copy_Mat4(aBone[ibone].affmat3Global, m01.mat);
      continue;
    }
    assert(ibone_p < (int) ibone);
    MatMat4(aBone[ibone].affmat3Global,
            aBone[ibone_p].affmat3Global, m01.mat);
  }
}

void delfem2::SetCurrentBoneRotationAsDefault(
    std::vector<CRigBone> &aBone) {
  const size_t nb = aBone.size();
  std::vector<double> aRot(nb * 16);
  for (unsigned int ib = 0; ib < nb; ++ib) {
    const int ibp = aBone[ib].ibone_parent;
    if (ibp == -1) {
       Mat4_Identity(aRot.data() + ib * 16);
      continue;
    } else {
      assert(ibp < (int) ib);
      double R0[16];
      Mat4_AffineQuaternion(R0, aBone[ibp].quatRelativeRot);
      MatMat4(aRot.data() + ib * 16, aRot.data() + ibp * 16, R0);
    }
  }
  // change from here
  for (unsigned int ib = 0; ib < nb; ++ib) {
    double Rv[3];
    Vec3_Mat4Vec3_Affine(Rv, aRot.data() + ib * 16, aBone[ib].transRelative);
    aBone[ib].transRelative[0] = Rv[0];
    aBone[ib].transRelative[1] = Rv[1];
    aBone[ib].transRelative[2] = Rv[2];
  }
  for (unsigned int ib = 0; ib < nb; ++ib) {
    double R1[16], R1B[16];
    Mat4_AffineQuaternion(R1, aBone[ib].quatRelativeRot);
    MatMat4(R1B, R1, aBone[ib].invBindMat);
    MatMat4(aBone[ib].invBindMat, aRot.data() + ib * 16, R1B);
  }
  for (unsigned int ib = 0; ib < nb; ++ib) {
    Quat_Identity(aBone[ib].quatRelativeRot);
  }
  UpdateBoneRotTrans(aBone);
}

DFM2_INLINE void delfem2::Skinning_LBS(
    std::vector<double> &aXYZ1,
    const std::vector<double> &aXYZ0,
    const std::vector<CRigBone> &aBone,
    const std::vector<double> &aW) {
  const size_t nBone = aBone.size();
  const size_t nP = aXYZ0.size() / 3;
  aXYZ1.resize(aXYZ0.size());
  assert(aW.size() == nBone * nP);
  for (unsigned int ip = 0; ip < nP; ++ip) {
    const double *p0 = aXYZ0.data() + ip * 3;
    double *p1 = aXYZ1.data() + ip * 3;
    p1[0] = 0.0;
    p1[1] = 0.0;
    p1[2] = 0.0;
    for (unsigned int ibone = 0; ibone < nBone; ++ibone) {
      double p2[3];
      aBone[ibone].DeformSkin(p2, p0);
      p1[0] += aW[ip * nBone + ibone] * p2[0];
      p1[1] += aW[ip * nBone + ibone] * p2[1];
      p1[2] += aW[ip * nBone + ibone] * p2[2];
    }
  }
}

// ---------

DFM2_INLINE void delfem2::SkinningSparse_LBS(
    std::vector<double> &aXYZ1,
    const std::vector<double> &aXYZ0,
    const std::vector<CRigBone> &aBone,
    const std::vector<double> &aWBoneSparse,
    const std::vector<unsigned> &aIdBoneSparse) {
//  const size_t nBone = aBone.size();
  const size_t nP = aXYZ0.size() / 3;
  const size_t nBW = aWBoneSparse.size() / nP;
  assert(aWBoneSparse.size() == nBW * nP);
  assert(aIdBoneSparse.size() == nBW * nP);
  aXYZ1.resize(aXYZ0.size());
  for (unsigned int ip = 0; ip < nP; ++ip) {
    const double *p0 = aXYZ0.data() + ip * 3;
    double *p1 = aXYZ1.data() + ip * 3;
    p1[0] = 0.0;
    p1[1] = 0.0;
    p1[2] = 0.0;
    for (unsigned int ibw = 0; ibw < nBW; ++ibw) {
      const unsigned int ib0 = aIdBoneSparse[ip * nBW + ibw];
      const double w0 = aWBoneSparse[ip * nBW + ibw];
      double p2[3];
      aBone[ib0].DeformSkin(p2, p0);
      p1[0] += w0 * p2[0];
      p1[1] += w0 * p2[1];
      p1[2] += w0 * p2[2];
    }
  }
}

DFM2_INLINE void delfem2::Skinning_LBS_LocalWeight(
    double *aXYZ,
    const double *aXYZ0,
    size_t nXYZ,
    const std::vector<CRigBone> &aBone,
    const double *aRigWeight,
    const unsigned int *aRigJoint) {
  for (unsigned int ip = 0; ip < nXYZ; ++ip) {
    double pos0[4] = {aXYZ0[ip * 3 + 0], aXYZ0[ip * 3 + 1], aXYZ0[ip * 3 + 2], 1.0};
    double pos1[3] = {0, 0, 0};
    double sum_w = 0.0;
    for (int iij = 0; iij < 4; ++iij) {
      double w = aRigWeight[ip * 4 + iij];
      if (w < 1.0e-30) { continue; }
      unsigned int ij = aRigJoint[ip * 4 + iij];
      // std::cout << ip << " " << iij << " " << ij << " " << w << std::endl;
      sum_w += w;
      assert (ij < aBone.size());
      double pos0a[4], pos0b[4];
      MatVec4(pos0a, aBone[ij].invBindMat, pos0);
      MatVec4(pos0b, aBone[ij].affmat3Global, pos0a);
      pos1[0] += w * pos0b[0];
      pos1[1] += w * pos0b[1];
      pos1[2] += w * pos0b[2];
    }
    assert(fabs(sum_w) > 1.0e-10);
    pos1[0] /= sum_w;
    pos1[1] /= sum_w;
    pos1[2] /= sum_w;
    aXYZ[ip * 3 + 0] = pos1[0];
    aXYZ[ip * 3 + 1] = pos1[1];
    aXYZ[ip * 3 + 2] = pos1[2];
  }
}

// ----------------------------------

DFM2_INLINE void
delfem2::InitBones_JointPosition(
    std::vector<CRigBone> &aBone,
    unsigned int nBone,
    const unsigned int *aIndBoneParent,
    const double *aJntPos0) {
  aBone.resize(nBone);
  for (unsigned int ib = 0; ib < nBone; ++ib) {
    unsigned int ibp = aIndBoneParent[ib];
    aBone[ib].ibone_parent = ibp;
    aBone[ib].invBindMat[3] = -aJntPos0[ib * 3 + 0];
    aBone[ib].invBindMat[7] = -aJntPos0[ib * 3 + 1];
    aBone[ib].invBindMat[11] = -aJntPos0[ib * 3 + 2];
    if (ibp != UINT_MAX) {
      aBone[ib].transRelative[0] = +aJntPos0[ib * 3 + 0] - aJntPos0[ibp * 3 + 0];
      aBone[ib].transRelative[1] = +aJntPos0[ib * 3 + 1] - aJntPos0[ibp * 3 + 1];
      aBone[ib].transRelative[2] = +aJntPos0[ib * 3 + 2] - aJntPos0[ibp * 3 + 2];
    } else {
      aBone[ib].transRelative[0] = +aJntPos0[ib * 3 + 0];
      aBone[ib].transRelative[1] = +aJntPos0[ib * 3 + 1];
      aBone[ib].transRelative[2] = +aJntPos0[ib * 3 + 2];
    }
  }
  UpdateBoneRotTrans(aBone);
}

// -------------------------------------------

DFM2_INLINE void
delfem2::SetMat4AffineBone_FromJointRelativeRotation(
    std::vector<double> &aMat4AffineBone,
    const double trans_root[3],
    const std::vector<double> &aQuatRelativeRot,
    const std::vector<int> &aIndBoneParent,
    const std::vector<double> &aJntPos0) {
  const size_t nBone = aIndBoneParent.size();
  assert(nBone >= 1);
  assert(aMat4AffineBone.size() == nBone * 16);
  Mat4_ScaleRotTrans(aMat4AffineBone.data(),
                     1.0, aQuatRelativeRot.data(), trans_root);
  for (unsigned int ibone = 1; ibone < nBone; ++ibone) {
    int ibp = aIndBoneParent[ibone];
    assert(ibp >= 0 && ibp < (int) nBone);
    // inv binding mat
    const double p1[3] = {
        aJntPos0[ibone * 3 + 0],
        aJntPos0[ibone * 3 + 1],
        aJntPos0[ibone * 3 + 2]};
    CMat4<double> M0, M1, M2;
    M0.SetAffineTranslate(-p1[0], -p1[1], -p1[2]);
    Mat4_AffineQuaternion(M1.mat,
                          aQuatRelativeRot.data() + ibone * 4);
    M2.SetAffineTranslate(+p1[0], +p1[1], +p1[2]);
    const CMat4<double> M4 = M2 * M1 * M0;
    MatMat4(aMat4AffineBone.data() + ibone * 16,
            aMat4AffineBone.data() + ibp * 16,
            M4.mat);
  }
}

DFM2_INLINE void delfem2::Rig_SkinReferncePositionsBoneWeighted(
    std::vector<double> &aRefPosAff,  // [ np, nBone*4 ]
    const std::vector<delfem2::CRigBone> &aBone1,
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aW) {
  const size_t np = aXYZ0.size() / 3;
  const size_t nb = aBone1.size();
  aRefPosAff.resize(np * nb * 4);
  for (unsigned int ip = 0; ip < np; ++ip) {
    double p0a[4] = {aXYZ0[ip * 3 + 0], aXYZ0[ip * 3 + 1], aXYZ0[ip * 3 + 2], 1.0};
    for (unsigned int ib = 0; ib < nb; ++ib) {
      double p0b[4];
      delfem2::MatVec4(p0b,
                       aBone1[ib].invBindMat, p0a);
      aRefPosAff[ip * (nb * 4) + ib * 4 + 0] = aW[ip * nb + ib] * p0b[0];
      aRefPosAff[ip * (nb * 4) + ib * 4 + 1] = aW[ip * nb + ib] * p0b[1];
      aRefPosAff[ip * (nb * 4) + ib * 4 + 2] = aW[ip * nb + ib] * p0b[2];
      aRefPosAff[ip * (nb * 4) + ib * 4 + 3] = aW[ip * nb + ib];
    }
  }
}



