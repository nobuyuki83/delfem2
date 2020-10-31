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

namespace delfem2{
namespace cnpy{

DFM2_INLINE void LoadSmpl(
    std::vector<double>& aXYZ0,
    std::vector<double>& aW,
    std::vector<unsigned int>& aTri,
    std::vector<int>& aIndBoneParent,
    std::vector<double>& aJntRgrs,
    const std::string& fpath);
  
}
}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/cnpy/smpl_cnpy.cpp"
#endif

#endif /* DFM2_SMPL_CNPY_H */