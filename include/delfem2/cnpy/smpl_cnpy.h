//
//  smplio.hpp
//  000_OpenWin
//
//  Created by Nobuyuki Umetani on 2020-03-11.
//

#ifndef DFM2_SMPL_CNPY_H
#define DFM2_SMPL_CNPY_H

#include "delfem2/dfm2_inline.h"
#include <cstdio>
#include <vector>
#include <string>

namespace delfem2{
namespace cnpy{

/**
 * @brief Load information for SMPL model
 * @param[out] aXYZ0 template position for vertices
 * @param[out] aW skinning weight
 * @param[out] aTri triangle index
 * @param[out] aIndBoneParent id of parent bone
 * @param[out] aJntRgrs wait of the body surface to obtain joint position
 * @param[in] fpath path
 */
DFM2_INLINE void LoadSmpl_Bone(
    std::vector<double>& aXYZ0,
    std::vector<double>& aW,
    std::vector<unsigned int>& aTri,
    std::vector<int>& aIndBoneParent,
    std::vector<double>& aJntRgrs,
    const std::string& fpath);

DFM2_INLINE void LoadSmpl_BoneBlendshape(
    std::vector<double>& aXYZ0,
    std::vector<double>& aW,
    std::vector<unsigned int>& aTri,
    std::vector<int>& aIndBoneParent,
    std::vector<double>& aJntRgrs,
    std::vector<double>& aBlendShape,
    std::vector<double>& aBlendPose,
    const std::string& fpath);
  
}
}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/cnpy/smpl_cnpy.cpp"
#endif

#endif /* DFM2_SMPL_CNPY_H */