/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef RIG_V3Q_H
#define RIG_V3Q_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>

#include "delfem2/vec3.h"

class CVector3;

class CRigBone
{
public:
    CRigBone(){
      for(int i=0;i<16;++i){ invBindMat[i]=0.0; }
      invBindMat[ 0] = 1.0;
      invBindMat[ 5] = 1.0;
      invBindMat[10] = 1.0;
      invBindMat[15] = 1.0;
      //////////
      scale = 1;
      rot[0] = 1;
      rot[1] = 0;
      rot[2] = 0;
      rot[3] = 0;
      trans[0] = 0;
      trans[1] = 0;
      trans[2] = 0;
      ibone_parent = -1;
    }
  CVector3 Pos() const{
    return CVector3(Mat[3],Mat[7],Mat[11]);
  }
  void SetRotationBryant(double rx, double ry, double rz);
  void SetTranslation(double tx, double ty, double tz);
  int PickHandler(const CVector3& org,
                  const CVector3& dir,
                  double rad_handlr,
                  double tol) const;
  void AffineJoint(const double a[16]) const;
public:
  std::string name;
  int ibone_parent;
  /////
  double invBindMat[16];
  double Mat[16];
  //////
  double rot[4]; // rotation of the joint from parent joint (quaternion w,x,y,z)
  double trans[3]; // position of the joint position from parent joint
  double scale; // scale
};


void UpdateBoneRotTrans(std::vector<CRigBone>& aBone);

void PickBone(int& ibone_selected,
              int& ielem_selected,
              const std::vector<CRigBone>& aBone,
              const CVector3& src,
              const CVector3& dir,
              double rad_hndlr,
              double tol);


////////////////////////////////////////////////////////////////////////

void UpdateRigSkin(double* aXYZ,
                   const double* aXYZ0,
                   unsigned int nXYZ,
                   const unsigned int* aTri,
                   unsigned int nTri,
                   const std::vector<CRigBone>& aBone,
                   const double* aRigWeight,
                   const unsigned int* aRigJoint);

////////////////////////////////////////////////////////////////////////

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

void SetPose_BioVisionHierarchy(std::vector<CRigBone>& aBone,
                                const std::vector<CChannel_BioVisionHierarchy>& aChannelInfo,
                                const double *aVal);



/*
class CBoneGoal
{
public:
  CBoneGoal(){
    itype = 0; // 0:position, 1:on line
  }
  void GetGoalPos(double* pos_trg,
                  const double* org_rot,
                  const double* pos_cur) const;
public:
  int ibone;
  ////
  int itype;
  double pos[3];
  double dir[3];
};
 */

//void BoneOptimization(std::vector<CBone_RigMsh>& aBone,
//                      const std::vector<CBoneGoal>& aBoneGoal);

//void DrawBoneTarget(const std::vector<CBoneGoal>& aBoneGoal,
//                    const std::vector<CBone_RigMsh>& aBone);

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

#endif /* rigmesh_hpp */
