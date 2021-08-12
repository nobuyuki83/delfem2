//
//  smplio.hpp
//  000_OpenWin
//
//  Created by Nobuyuki Umetani on 2020-03-11.
//

#ifndef DFM2_SMPL_CNPY_H
#define DFM2_SMPL_CNPY_H

#include <cstdio>
#include <vector>
#include <string>

#include "delfem2/dfm2_inline.h"

namespace delfem2{
namespace cnpy{

/**
 * @brief Load information for SMPL model
 * @param[out] aXYZ0 template position for vertices
 * @param[out] aW skinning weight
 * @param[out] aTri triangle index
 * @param[out] aIndBoneParent id of parent bone
 * @param[out] aJntRgrs weight of the body surface to obtain joint position (column-major)
 * @param[in] fpath path
 */
DFM2_INLINE void LoadSmpl_Bone(
    std::vector<double>& aXYZ0,
    std::vector<double>& aW,
    std::vector<unsigned int>& aTri,
    std::vector<unsigned int>& aIndBoneParent,
    std::vector<double>& aJntRgrs,
    const std::string& fpath);

/**
 *
 * @param aXYZ0
 * @param aW
 * @param aTri
 * @param aIndBoneParent
 * @param[out] aJntRgrs weight of the body surface to obtain joint position (column-major)
 * @param aBlendShape
 * @param aBlendPose
 * @param fpath
 */
DFM2_INLINE void LoadSmpl_BoneBlendshape(
    std::vector<double>& aXYZ0,
    std::vector<double>& aW,
    std::vector<unsigned int>& aTri,
    std::vector<unsigned int>& aIndBoneParent,
    std::vector<double>& aJntRgrs,
    std::vector<double>& aBlendShape,
    std::vector<double>& aBlendPose,
    const std::string& fpath);
  
}
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/cnpy/smpl_cnpy.cpp"
#endif

#endif /* DFM2_SMPL_CNPY_H */