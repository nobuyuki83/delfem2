#ifndef VOXEL_H
#define VOXEL_H

#include <vector>

#include "vec3.h"

int Adj_Grid
(int ivox_picked, int iface_picked,
 int ndivx, int ndivy, int ndivz);

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

void Draw_CubeGrid
(bool is_picked, int iface_picked,
 double elen, const CVector3& org,
 const CCubeGrid& cube);

void Pick_CubeGrid
(int& icube_pic, int& iface_pic,
 const CVector3& src_pic, const CVector3& dir_pic,
 double elen,
 const CVector3& org,
 const std::vector<CCubeGrid>& aCube);

void Adj_CubeGrid
(int& ivx, int& ivy, int& ivz,
 int ivox, int iface,
 std::vector<CCubeGrid>& aCube);

void Add_CubeGrid
(std::vector<CCubeGrid>& aVox,
 int ivx1, int ivy1, int ivz1);

void Del_CubeGrid
(std::vector<CCubeGrid>& aCube,
 int i1, int j1, int k1);

#endif
