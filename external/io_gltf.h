#ifndef IO_GLTF_H
#define IO_GLTF_H

#include <vector>
#include <stdio.h>
#include <memory>

#include "delfem2/rig_v3q.h" // CRigBone is defined here

namespace tinygltf{
  class Model;
};

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
  void GetBone(std::vector<CRigBone>& aBone,
               int iskin) const;
public:
  tinygltf::Model* model;
};


#endif
