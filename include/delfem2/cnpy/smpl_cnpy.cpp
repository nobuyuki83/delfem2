//
//  smplio.cpp
//  000_OpenWin
//
//  Created by Nobuyuki Umetani on 2020-03-11.
//
#include <iostream>
#include <vector>
#include "cnpy/cnpy.h"

#include "smpl_cnpy.h"

namespace dfm2 = delfem2;

void dfm2::cnpy::LoadSmpl(
    std::vector<double>& aXYZ0,
    std::vector<double>& aW,
    std::vector<unsigned int>& aTri,
    std::vector<int>& aIndBoneParent,
    std::vector<double>& aJntRgrs,
    const std::string& fpath)
{
  ::cnpy::npz_t my_npz = ::cnpy::npz_load(fpath);
  const unsigned int nP = my_npz["vertices_template"].shape[0];
  assert( my_npz["vertices_template"].shape[1] == 3 );
  aXYZ0 = my_npz["vertices_template"].as_vec<double>();
  {
    ::cnpy::NpyArray& npT = my_npz["face_indices"];
    assert( my_npz["face_indices"].shape[1] == 3 );
    aTri = npT.as_vec<unsigned>();
    for(auto &i: aTri){ i -= 1; }
  }
  const unsigned int nBone = my_npz["weights"].shape[1];
  assert( my_npz["weights"].shape[0] == nP );
  aW = my_npz["weights"].as_vec<double>();
  {
    const ::cnpy::NpyArray& npT = my_npz["kinematic_tree"];
    const int* tree = npT.data<int>();
    aIndBoneParent.assign(tree, tree+nBone);
  }
  assert( my_npz["joint_regressor"].shape[0] == nBone );
  assert( my_npz["joint_regressor"].shape[1] == nP );
  std::cout << my_npz["joint_regressor"].fortran_order << std::endl;
  aJntRgrs = my_npz["joint_regressor"].as_vec<double>();
  assert( aJntRgrs.size() == nBone*nP );
}


