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

void SetBone(std::vector<delfem2::CRigBone>& aBone,
             const tinygltf::Model& model,
             int inode_cur, int ibone_p,
             const std::vector<int>& mapNode2Bone);

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
  void GetBone(std::vector<delfem2::CRigBone>& aBone,
               int iskin) const;
public:
  tinygltf::Model* model;
};

}


#endif
