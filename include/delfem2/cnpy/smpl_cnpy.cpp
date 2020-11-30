/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/cnpy/smpl_cnpy.h"
#include "cnpy/cnpy.h"
#include <vector>

namespace delfem2{
namespace smpl_cnpy{

DFM2_INLINE void LoadSmpl_Pose_FromCnpy(
    std::vector<double>& aXYZ0,
    std::vector<double>& aW,
    std::vector<unsigned int>& aTri,
    std::vector<unsigned int>& aIndBoneParent,
    std::vector<double>& aJntRgrs,
    ::cnpy::npz_t& my_npz)
{
  const unsigned int nP = my_npz["vertices_template"].shape[0]; // number of points
  assert( my_npz["vertices_template"].shape[1] == 3 );
  aXYZ0 = my_npz["vertices_template"].as_vec<double>();
  {
    ::cnpy::NpyArray& npT = my_npz["face_indices"];
    assert( npT.shape.size() == 2 );
    assert( npT.shape[1] == 3 );
    aTri = npT.as_vec<unsigned>();
    for(auto &indp: aTri){ indp -= 1; }
  }
  const unsigned int nBone = my_npz["weights"].shape[1];
  assert( my_npz["weights"].shape[0] == nP );
  aW = my_npz["weights"].as_vec<double>();
  { // bone tree information
    const ::cnpy::NpyArray& npT = my_npz["kinematic_tree"];
    assert( npT.shape.size() == 2 && npT.shape[0] == 2 && npT.shape[1] == nBone );
#ifndef NDEBUG
    for(unsigned int ib=0;ib<nBone;++ib){
      const int jb0 = npT.data<int>()[ib];
      const int ib0 = npT.data<int>()[ib+nBone];
      assert( ib0 == (int)ib );
      assert( jb0 == -1 || (jb0<(int)ib) );
    }
#endif
    const int* tree = npT.data<int>();
    aIndBoneParent.assign(tree, tree+nBone);
  }
  assert( my_npz["joint_regressor"].shape[0] == nBone );
  assert( my_npz["joint_regressor"].shape[1] == nP );
  assert( my_npz["joint_regressor"].fortran_order == 1 );
  aJntRgrs = my_npz["joint_regressor"].as_vec<double>();
  assert( aJntRgrs.size() == nBone*nP );
}

}
}

DFM2_INLINE void delfem2::cnpy::LoadSmpl_Bone(
    std::vector<double>& aXYZ0,
    std::vector<double>& aW,
    std::vector<unsigned int>& aTri,
    std::vector<unsigned int>& aIndBoneParent,
    std::vector<double>& aJntRgrs,
    const std::string& fpath)
{
  ::cnpy::npz_t my_npz = ::cnpy::npz_load(fpath);
  delfem2::smpl_cnpy::LoadSmpl_Pose_FromCnpy(
      aXYZ0,aW,aTri,aIndBoneParent,aJntRgrs,
      my_npz);
}


DFM2_INLINE void delfem2::cnpy::LoadSmpl_BoneBlendshape(
    std::vector<double>& aXYZ0,
    std::vector<double>& aW,
    std::vector<unsigned int>& aTri,
    std::vector<unsigned int>& aIndBoneParent,
    std::vector<double>& aJntRgrs,
    std::vector<double>& aBlendShape,
    std::vector<double>& aBlendPose,
    const std::string& fpath)
{
  ::cnpy::npz_t my_npz = ::cnpy::npz_load(fpath);
  delfem2::smpl_cnpy::LoadSmpl_Pose_FromCnpy(
      aXYZ0,aW,aTri,aIndBoneParent,aJntRgrs,
      my_npz);
  const unsigned int np0 = aXYZ0.size()/3;
  {
    ::cnpy::NpyArray& npS = my_npz["shape_blend_shapes"];
    assert( npS.shape.size() == 3 && npS.shape[0] == np0 && npS.shape[1] == 3 && npS.shape[2] == 10 );
    assert( npS.fortran_order == 0 );
    aBlendShape = npS.as_vec<double>();
  }
  {
    ::cnpy::NpyArray& npP = my_npz["pose_blend_shapes"];
    assert( npP.shape.size() == 3 && npP.shape[0] == np0 && npP.shape[1] == 3 && npP.shape[2] == 23*9 );
    aBlendPose = npP.as_vec<double>();
  }

}