//
//  smplio.hpp
//  000_OpenWin
//
//  Created by Nobuyuki Umetani on 2020-03-11.
//

#ifndef smplio_hpp
#define smplio_hpp

#include <stdio.h>

namespace delfem2{
namespace cnpy{

void LoadSmpl(std::vector<double>& aXYZ0,
              std::vector<double>& aW,
              std::vector<unsigned int>& aTri,
              std::vector<int>& aIndBoneParent,
              std::vector<double>& aJntRgrs,
              const std::string& fpath);

  
}
}

#endif /* smplio_hpp */
