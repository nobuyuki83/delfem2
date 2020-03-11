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
    std::vector<double>& aJntPos0)
{
  ::cnpy::npz_t my_npz = ::cnpy::npz_load(std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
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
  const std::vector<double> aJntRgrs = my_npz["joint_regressor"].as_vec<double>();
  assert( aJntRgrs.size() == nBone*nP );
  
  {
    const unsigned int nP = aXYZ0.size()/3;
    const unsigned int nBone = aIndBoneParent.size();
    aJntPos0.assign(nBone*3, 0.0);
    for(int ib=0;ib<nBone;++ib){
      aJntPos0[ib*3+0] = 0;
      aJntPos0[ib*3+1] = 0;
      aJntPos0[ib*3+2] = 0;
      for(int ip=0;ip<nP;++ip){
        aJntPos0[ib*3+0] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+0];
        aJntPos0[ib*3+1] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+1];
        aJntPos0[ib*3+2] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+2];
      }
    }
  }
}
