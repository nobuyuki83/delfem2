/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_CUBEGRID_H
#define DFM2_CUBEGRID_H

#include <vector>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

class CCubeGrid
{
public:
  CCubeGrid(){
    this->ivx = 0;
    this->ivy = 0;
    this->ivz = 0;
    is_active = true;
  }
  CCubeGrid(int i, int j, int k){
    this->ivx = i;
    this->ivy = j;
    this->ivz = k;
    is_active = true;
  }
  bool operator<(const CCubeGrid& rhs) const
  {
    if( this->ivx != rhs.ivx ){ return this->ivx < rhs.ivx; }
    if( this->ivy != rhs.ivy ){ return this->ivy < rhs.ivy; }
    if( this->ivz != rhs.ivz ){ return this->ivz < rhs.ivz; }
    return false;
  }
public:
  int ivx, ivy, ivz;
  bool is_active;
};


void Pick_CubeGrid(
    unsigned int& icube_pic,
    int& iface_pic,
    const double src_pic[3],
    const double dir_pic_[3],
    double elen,
    const double org[3],
    const std::vector<CCubeGrid>& aCube);

void Adj_CubeGrid(
    int& ivx,
    int& ivy,
    int& ivz,
    unsigned int ivox,
    int iface,
    std::vector<CCubeGrid>& aCube);

void Add_CubeGrid(
    std::vector<CCubeGrid>& aVox,
    int ivx1,
    int ivy1,
    int ivz1);

void Del_CubeGrid(
    std::vector<CCubeGrid>& aCube,
    int i1,
    int j1,
    int k1);

void AABB_CubeGrid(
    int aabb[6],
    const std::vector<CCubeGrid>& aCube);

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/gridcube.cpp"
#endif

#endif /* cubegrid_h */
