/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_RIG_GEO3_H
#define DFM2_RIG_GEO3_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>

#include "delfem2/geo3_v23m34q.h"
#include "delfem2/dfm2_inline.h"

namespace delfem2 {

DFM2_INLINE void Transpose_Mat(
    std::vector<double> &At,
    const std::vector<double> &A,
    unsigned int ncol,
    unsigned int nrow);

/**
 *
 * @param[out] aWBone_RigSparse
 * @param[out] aIdBone_RigSparse
 * @param[in] aW
 * @param[in] np
 * @param[in] nb
 * @param[in] thres
 */
DFM2_INLINE void SparsifyMatrixRow(
    std::vector<double> &aSparseA,
    std::vector<unsigned int> &aSparseInd,
    const double *pA,
    unsigned int nrow,
    unsigned int ncol,
    double thres);

DFM2_INLINE void Points3_WeighttranspPosition(
    std::vector<double> &aPos0,
    const std::vector<double> &aWeighttransp,
    const std::vector<double> &aXYZ0);

DFM2_INLINE void Points3_WeightsparsePosition(
    std::vector<double> &aPos0,
    unsigned int nPos0,
    const std::vector<double> &aSparseW,
    const std::vector<unsigned int> &aSparseIp,
    const std::vector<double> &aXYZ0);

/**
 * @brief Set 3D affine matrix that transfrom from intial position from the deformed poisition for each bones.
 * @details this funcition is for rigging without the class "CRigBone"
 */
DFM2_INLINE void SetMat4AffineBone_FromJointRelativeRotation(
    std::vector<double> &aMat4AffineBone,
    const double trans_root[3],
    const std::vector<double> &aQuatRelativeRot,
    const std::vector<int> &aIndBoneParent,
    const std::vector<double> &aJntPos0);

// ----------------------------------------------
// CRigBone

/**
 *@brief articulated rigid body for character rigging
 *@details this class support rig in GlTF and BioVision BVH
 */
class CRigBone {
 public:
  CRigBone()
      : invBindMat{}, transRelative{}, quatRelativeRot{}, affmat3Global{} {
    Mat4_Identity(invBindMat);
    Quat_Identity(quatRelativeRot);
    scale = 1;
    transRelative[0] = 0;
    transRelative[1] = 0;
    transRelative[2] = 0;
    ibone_parent = -1;
  }

  // TODO: rename this function PosRoot
  delfem2::CVec3d Pos() const {
    // invBindMat will map the global position to the origin.
    // the translation component of affmat3Global represent the root position after skeletal deformation
    return delfem2::CVec3d(affmat3Global[3], affmat3Global[7], affmat3Global[11]);
  }

  void DeformSkin(double pos2[3],
                  const double pos0[3]) const;

  void SetRotationBryant(
      double rx, double ry, double rz);

  void SetTranslation(
      double tx, double ty, double tz);

 public:
  std::string name; // initialized and stay constant

  /** 
   * @brief Inverse of Affine matrix to send the skin to the bone reference config.
   * @details The joint position of this bone will be mapped to the origin.
   * The data format is column-major
   */
  double invBindMat[16];

  /**
   * @brief index of parent bone
   * @brief initialized and stay constant
   */
  int ibone_parent;

  /**
   * position of the joint position from parent joint
   */
  double transRelative[3];

  double scale; // scale

  /**
   * @brief rotation of the joint from parent joint (quaternion x,y,z,w).
   * @details this value will be changed  when pose is edited
   */
  double quatRelativeRot[4];

  /**
   * @brief affine matrix  to send bone from the origin to the deformed pose
   * @details this value will be set when the pose is edited using the function
   * "void UpdateBoneRotTrans(std::vector<CRigBone>& aBone)"
   */
  double affmat3Global[16];
};

/**
 * @brief set CRigBone.affmat3Global for all the bones
 * @details transformation goes in the order of "scale(S)" -> rotation(R) -> translation(T)
 * Affine matrix of i-th bone will be computed as A_i = T_i*R_i*S_i * T_i-1*R_i-1*S_i-1 * .... * T_0*R_0*S_0
 * The point is transfromed for i-th bone as V = A_i*B_i*v where B_i is the "inverseBindingMatrix"
 */
DFM2_INLINE void UpdateBoneRotTrans(
    std::vector<CRigBone> &aBone);

DFM2_INLINE void SetCurrentBoneRotationAsDefault(
    std::vector<CRigBone> &aBone);

DFM2_INLINE void PickBone(
    int &ibone_selected,
    int &ielem_selected,
    const std::vector<CRigBone> &aBone,
    const delfem2::CVec3d &src,
    const delfem2::CVec3d &dir,
    double rad_hndlr,
    double tol);


// ----------------------------------

DFM2_INLINE void Skinning_LBS_LocalWeight(
    double *aXYZ,
    const double *aXYZ0,
    size_t nXYZ,
    const std::vector<CRigBone> &aBone,
    const double *aRigWeight,
    const unsigned int *aRigJoint);

/**
 * @params aW rigging weight [np, nbone]
 */
DFM2_INLINE void Skinning_LBS(
    std::vector<double> &aXYZ1,
    const std::vector<double> &aXYZ0,
    const std::vector<CRigBone> &aBone,
    const std::vector<double> &aW);

/**
 *
 * @param[out] aXYZ1
 * @param[in] aXYZ0
 * @param[in] aBone
 * @param[in] aWBoneSparse
 * @param[in] aIdBoneSparse
 */
DFM2_INLINE void SkinningSparse_LBS(
    std::vector<double> &aXYZ1,
    const std::vector<double> &aXYZ0,
    const std::vector<CRigBone> &aBone,
    const std::vector<double> &aWBoneSparse,
    const std::vector<unsigned> &aIdBoneSparse);

// --------------------------------------

DFM2_INLINE void InitBones_JointPosition(
    std::vector<CRigBone> &aBone,
    unsigned int nb,
    const unsigned int *aIndBoneParent,
    const double *aJntPos0);

// --------------------------------------

DFM2_INLINE void Rig_SkinReferncePositionsBoneWeighted(
    std::vector<double> &aRefPosAff,
    const std::vector<CRigBone> aBone1,
    const std::vector<double> &aXYZ0,
    const std::vector<double> &aW);

// ------------------------------------

class CChannel_BioVisionHierarchy {
 public:
  CChannel_BioVisionHierarchy(unsigned int ib, int ia, bool br) {
    this->ibone = ib;
    this->iaxis = ia;
    this->isrot = br;
  }

 public:
  unsigned int ibone;
  int iaxis;
  bool isrot;
};

DFM2_INLINE void Read_BioVisionHierarchy(
    std::vector<CRigBone> &aBone,
    std::vector<CChannel_BioVisionHierarchy> &aChannelInfo,
    int &nframe,
    std::vector<double> &aChannelValue,
    const std::string &path_bvh);

/**
 * @brief set value to CRigBone.rot (bone rotation from parent bone)
 */
DFM2_INLINE void SetPose_BioVisionHierarchy(
    std::vector<CRigBone> &aBone,
    const std::vector<CChannel_BioVisionHierarchy> &aChannelInfo,
    const double *aVal);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/rig_geo3.cpp"
#endif

#endif

