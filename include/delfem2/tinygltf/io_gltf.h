#ifndef IO_GLTF_H
#define IO_GLTF_H

#include <vector>
#include <stdio.h>
#include <memory>

#include "delfem2/rig_geo3.h" // CRigBone is defined here

namespace tinygltf{
  class Model;
};

namespace delfem2 {

void Print(const tinygltf::Model& model);

void GetMeshInfo(std::vector<double>& aXYZ,
                 std::vector<unsigned int>& aTri,
                 std::vector<double>& aRigWeight,
                 std::vector<unsigned int>& aRigJoint,
                 const tinygltf::Model& model,
                 int imsh, int iprimitive);

void GetBoneBinding(std::vector<delfem2::CRigBone>& aBone,
                    const tinygltf::Model& model);

/**
 * given a node idx of root bone, it recursively search the tree brunch to construct bone tree
 * @param[out] aBone
 * @param[in] model
 * @param[in] inode_cur
 * @param[in] ibone_p
 * @param[in] mapNode2Bone
 */
void SetBone(
    std::vector<delfem2::CRigBone>& aBone,
    const tinygltf::Model& model,
    unsigned int inode_cur,
    int ibone_p,
    const std::vector<unsigned int>& mapNode2Bone);

void GetBone(
    std::vector<delfem2::CRigBone>& aBone,
    const tinygltf::Model& model,
    unsigned int iskin);

class CGLTF
{
public:
  bool Read(const std::string& fpath);
  void Print() const;
  void GetMeshInfo(std::vector<double>& aXYZ0,
                   std::vector<unsigned int>& aTri,
                   std::vector<double>& aRigWeight,
                   std::vector<unsigned int>& aRigJoint,
                   int imesh, int iprimitive) const;
  void GetBone(
      std::vector<delfem2::CRigBone>& aBone,
      unsigned int iskin) const;
public:
  tinygltf::Model* pModel = nullptr;
};

}


#endif
