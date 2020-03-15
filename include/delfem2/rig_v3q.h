/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_RIG_V3Q_H
#define DFM2_RIG_V3Q_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>

#include "delfem2/vec3.h"

namespace delfem2 {

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
      trans[0] = 0;
      trans[1] = 0;
      trans[2] = 0;
      ibone_parent = -1;
    }
  delfem2::CVec3d Pos() const{
    return delfem2::CVec3d(affmat3Global[3],affmat3Global[7],affmat3Global[11]);
  }
  void SetRotationBryant(double rx, double ry, double rz);
  void SetTranslation(double tx, double ty, double tz);
  int PickHandler(const delfem2::CVec3d& org,
                  const delfem2::CVec3d& dir,
                  double rad_handlr,
                  double tol) const;
  void AffineJoint(const double a[16]) const;
public:
  std::string name; // initialized and stay constant
  
  /** 
   * @details Inverse of Affine matrix to send this bone to the origin and reference config
   */
  double invBindMat[16];
  
  int ibone_parent; // initialized and stay constant
  
  double trans[3]; // position of the joint position from parent joint
  
  double scale; // scale
  
  /**
   * @brief rotation of the joint from parent joint (quaternion w,x,y,z).
   * @details this value will be changed  when pose is edited
   */
  double quatRelativeRot[4];

  /**
   * @brief affine matrix  to send bone in the origin to the deformed pose
   * @details this value will be set when the pose is edited using the function
   * "void UpdateBoneRotTrans(std::vector<CRigBone>& aBone)"
   */
  double affmat3Global[16];
};





/**
 * @brief set "CRgidiBone.affmat3Global" based on "CRigidBone.quadRelativeRot"
 */
void UpdateBoneRotTrans(std::vector<CRigBone>& aBone);

void PickBone(int& ibone_selected,
              int& ielem_selected,
              const std::vector<CRigBone>& aBone,
              const delfem2::CVec3d& src,
              const delfem2::CVec3d& dir,
              double rad_hndlr,
              double tol);


// ----------------------------------

void UpdateRigSkin(double* aXYZ,
                   const double* aXYZ0,
                   unsigned int nXYZ,
                   const unsigned int* aTri,
                   unsigned int nTri,
                   const std::vector<CRigBone>& aBone,
                   const double* aRigWeight,
                   const unsigned int* aRigJoint);

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

void Read_BioVisionHierarchy(std::vector<CRigBone>& aBone,
                             std::vector<CChannel_BioVisionHierarchy>& aChannelInfo,
                             int& nframe,
                             std::vector<double>& aChannelValue,
                             const std::string& path_bvh);

/**
 * @brief set value to CRigBone.rot (bone rotation from parent bone)
 */
void SetPose_BioVisionHierarchy(std::vector<CRigBone>& aBone,
                                const std::vector<CChannel_BioVisionHierarchy>& aChannelInfo,
                                const double *aVal);

/**
 * @brief Set 3D affine matrix that transfrom from intial position from the deformed poisition for each bones.
 */
void SetMat4AffineBone_FromJointRelativeRotation(
    std::vector<double>& aMat4AffineBone,
    const double trans_root[3],
    const std::vector<double>& aQuatRelativeRot,
    const std::vector<int>& aIndBoneParent,
    const std::vector<double>& aJntPos0);


} // namespace delfem2

#endif // #define DFM2_RIG_V3Q_H

