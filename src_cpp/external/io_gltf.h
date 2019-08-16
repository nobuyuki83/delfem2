#ifndef IO_GLTF_H
#define IO_GLTF_H

#include <vector>
#include <stdio.h>

#include "delfem2/rig_v3q.h" // CRigBone is defined here


void Print(const tinygltf::Model& model);

void GetMeshInfo(std::vector<double>& aXYZ,
                 std::vector<unsigned int>& aTri,
                 std::vector<double>& aRigWeight,
                 std::vector<unsigned int>& aRigJoint,
                 const tinygltf::Model& model,
                 int imsh, int iprimitive);

void GetBoneBinding(std::vector<CRigBone>& aBone,
                    const tinygltf::Model& model);

void SetBone(std::vector<CRigBone>& aBone,
             const tinygltf::Model& model,
             int inode_cur, int ibone_p,
             const std::vector<int>& mapNode2Bone);


#endif
