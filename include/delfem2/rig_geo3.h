/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_RIG_V3Q_H
#define DFM2_RIG_V3Q_H

#include "delfem2/dfm2_inline.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include "delfem2/vec3.h"

namespace delfem2 {

/**
 * @brief Set 3D affine matrix that transfrom from intial position from the deformed poisition for each bones.
 * @details this funcition is for rigging without the class "CRigBone"
 */
DFM2_INLINE void SetMat4AffineBone_FromJointRelativeRotation(
    std::vector<double>& aMat4AffineBone,
    const double trans_root[3],
    const std::vector<double>& aQuatRelativeRot,
    const std::vector<int>& aIndBoneParent,
    const std::vector<double>& aJntPos0);

// ----------------------------------------------
// CRigBone

/**
 *@brief articulated rigid body for character rigging
 *@details this class support rig in GlTF and BioVision BVH
 */
class CRigBone
{
public:
    CRigBone(){
      for(int i=0;i<16;++i){ invBindMat[i]=0.0; }
      invBindMat[ 0] = 1.0;
      invBindMat[ 5] = 1.0;
      invBindMat[10] = 1.0;
      invBindMat[15] = 1.0;
      //
      scale = 1;
      quatRelativeRot[0] = 1;
      quatRelativeRot[1] = 0;
      quatRelativeRot[2] = 0;
      quatRelativeRot[3] = 0;
      transRelative[0] = 0;
      transRelative[1] = 0;
      transRelative[2] = 0;
      ibone_parent = -1;
    }
  delfem2::CVec3d Pos() const{
    return delfem2::CVec3d(affmat3Global[3],affmat3Global[7],affmat3Global[11]);
  }
  void DeformSkin(double pos2[3],
                  const double pos0[3]) const;
  void SetRotationBryant(double rx, double ry, double rz);
  void SetTranslation(double tx, double ty, double tz);
  /*
  int PickHandler(const delfem2::CVec3d& org,
                  const delfem2::CVec3d& dir,
                  double rad_handlr,
                  double tol) const;
   */
  void AffineJoint(const double a[16]) const;
public:
  std::string name; // initialized and stay constant
  
  /** 
   * @details Inverse of Affine matrix to send the skin to the bone reference config. The joint position of this bone will be mapped to the origin
   */
  double invBindMat[16];
  
  int ibone_parent; // initialized and stay constant
  
  double transRelative[3]; // position of the joint position from parent joint
  
  double scale; // scale
  
  /**
   * @brief rotation of the joint from parent joint (quaternion w,x,y,z).
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
 * @brief set "CRgidiBone.affmat3Global" based on "CRigidBone.quadRelativeRot"
 */
DFM2_INLINE void UpdateBoneRotTrans(std::vector<CRigBone>& aBone);

DFM2_INLINE void PickBone(
    int& ibone_selected,
    int& ielem_selected,
    const std::vector<CRigBone>& aBone,
    const delfem2::CVec3d& src,
    const delfem2::CVec3d& dir,
    double rad_hndlr,
    double tol);


// ----------------------------------

DFM2_INLINE void Skinning_LBS_LocalWeight(
    double* aXYZ,
    const double* aXYZ0,
    unsigned int nXYZ,
    const unsigned int* aTri,
    unsigned int nTri,
    const std::vector<CRigBone>& aBone,
    const double* aRigWeight,
    const unsigned int* aRigJoint);

/**
 * @params aW rigging weight [np, nbone]
 */
DFM2_INLINE void Skinning_LBS(
    std::vector<double>& aXYZ1,
    const std::vector<double>& aXYZ0,
    const std::vector<CRigBone>& aBone,
    const std::vector<double>& aW);


// --------------------------------------

DFM2_INLINE void Smpl2Rig(
    std::vector<CRigBone>& aBone,
    const std::vector<int>& aIndBoneParent,
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aJntRgrs);

// --------------------------------------

DFM2_INLINE void Rig_SkinReferncePositionsBoneWeighted(
    std::vector<double>& aRefPosAff,
    const std::vector<CRigBone> aBone1,
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aW);

// ------------------------------------



class CTarget
{
public:
  unsigned int ib;
  CVec3d pos;
public:
  void WdW
  (std::vector<double>& aW,
   std::vector<double>& adW,
   const std::vector<CRigBone>& aBone,
   std::vector<double>& aL) const; // [ [nb, 3],  [ndim(3), nBone, ndim(4)] ]
};

DFM2_INLINE void Rig_SensitivityBoneTransform(double* aL, // [ ndim(3), nBone, ndim(4) ]
                                  unsigned int ib_s,
                                  unsigned int idim_s,
                                  const std::vector<CRigBone> aBone1);

DFM2_INLINE void Rig_SensitivityBoneTransform_Eigen(
    std::vector<double>& Lx, // [ [nBone, ndim(3)], [nBone, ndim(4)] ]
    std::vector<double>& Ly, // [ [nBone, ndim(3)], [nBone, ndim(4)] ]
    std::vector<double>& Lz, // [ [nBone, ndim(3)], [nBone, ndim(4)] ]
    unsigned int ib_s,
    double idim_s,
    bool is_rot,
    const std::vector<CRigBone> aBone1);

DFM2_INLINE void Rig_WdW_Target_Eigen(
    std::vector<double>& aW,
    std::vector<double>& adW,
    const std::vector<CRigBone>& aBone,
    const CTarget& target,
    const std::vector<double>& Lx,  // [ nsns, nBone*4 ]
    const std::vector<double>& Ly,  // [ nsns, nBone*4 ]
    const std::vector<double>& Lz); // [ nsns, nBone*4 ]

DFM2_INLINE void Rig_SensitivitySkin_BoneRotation_Eigen(
    std::vector<double>& dSkinX,
    std::vector<double>& dSkinY,
    std::vector<double>& dSkinZ,
    const std::vector<CRigBone>& aBone1,
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aW,
    const std::vector<double>& Lx, // [ nBone*3, nBone*4 ]
    const std::vector<double>& Ly, // [ nBone*3, nBone*4 ]
    const std::vector<double>& Lz); // [ nBone*3, nBone*4 ]

DFM2_INLINE void Rig_SensitivitySkin_BoneRotation(
    std::vector<double>& aSns, // np*ndim(3) * nb*3
    const std::vector<CRigBone> aBone1,
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aW,
    const std::vector<double>& aL);


// ------------------------------------

class CChannel_BioVisionHierarchy
{
public:
  CChannel_BioVisionHierarchy(int ib, int ia, bool br){
    this->ibone = ib;
    this->iaxis = ia;
    this->isrot = br;
  }
public:
  int ibone;
  int iaxis;
  bool isrot;
};

DFM2_INLINE void Read_BioVisionHierarchy(
    std::vector<CRigBone>& aBone,
    std::vector<CChannel_BioVisionHierarchy>& aChannelInfo,
    int& nframe,
    std::vector<double>& aChannelValue,
    const std::string& path_bvh);

/**
 * @brief set value to CRigBone.rot (bone rotation from parent bone)
 */
DFM2_INLINE void SetPose_BioVisionHierarchy(
    std::vector<CRigBone>& aBone,
    const std::vector<CChannel_BioVisionHierarchy>& aChannelInfo,
    const double *aVal);

} // namespace delfem2

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/rig_geo3.cpp"
#endif

#endif // #define DFM2_RIG_V3Q_H

