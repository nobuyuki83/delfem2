#ifndef VOXEL_H
#define VOXEL_H

#include <vector>

// TODO: remove dependency to vec3
#include "vec3.h"

int Adj_Grid
(int ivox_picked, int iface_picked,
 int ndivx, int ndivy, int ndivz);

void GetQuad_VoxelGrid
(std::vector<double>& aXYZ, std::vector<int>& aQuad,
 int ndivx, int ndivy, int ndivz,
 int iorgx, int iorgy, int iorgz,
 const std::vector<int>& aIsVox);

////////////////////////////////////////////////////////////////////

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

void Draw_CubeGrid(bool is_picked, int iface_picked,
                   double elen, const CVector3& org,
                   const CCubeGrid& cube);
void Pick_CubeGrid(int& icube_pic, int& iface_pic,
                   const CVector3& src_pic, const CVector3& dir_pic,
                   double elen,
                   const CVector3& org,
                   const std::vector<CCubeGrid>& aCube);
void Adj_CubeGrid(int& ivx, int& ivy, int& ivz,
                  int ivox, int iface,
                  std::vector<CCubeGrid>& aCube);
void Add_CubeGrid(std::vector<CCubeGrid>& aVox,
                  int ivx1, int ivy1, int ivz1);
void Del_CubeGrid(std::vector<CCubeGrid>& aCube,
                  int i1, int j1, int k1);
void AABB_CubeGrid(int aabb[6],
                   const std::vector<CCubeGrid>& aCube);

////////////////////////////////////////////////////////////////////////////////////


class CVoxelGrid
{
public:
  CVoxelGrid(){
    ndivx = ndivy = ndivz = 0;
    iorgx = iorgy = iorgz = 0;
    aIsVox.clear();
  }
  void Init_AABB(const int aabb[6]){
    ndivx = aabb[1]-aabb[0];
    ndivy = aabb[3]-aabb[2];
    ndivz = aabb[5]-aabb[4];
    iorgx = aabb[0];
    iorgy = aabb[2];
    iorgz = aabb[4];
    const int nvoxel = ndivx*ndivy*ndivz;
    aIsVox.assign(nvoxel,0);
  }
  void AABB(int aabb[8]) const {
    bool is_initial = true;
    aabb[0] = +1;
    aabb[1] = -1;
    for(int igvx=0;igvx<ndivx;++igvx){
      for(int igvy=0;igvy<ndivy;++igvy){
        for(int igvz=0;igvz<ndivz;++igvz){
          const int ivoxel = igvx*(ndivy*ndivz)+igvy*ndivz+igvz;
          if( aIsVox [ivoxel] ==0 ){ continue; }
          const int igpx0 = igvx+iorgx+0;  const int igpx1 = igvx+iorgx+1;
          const int igpy0 = igvy+iorgy+0;  const int igpy1 = igvy+iorgy+1;
          const int igpz0 = igvz+iorgz+0;  const int igpz1 = igvz+iorgz+1;
          if( is_initial ){
            aabb[0] = igpx0;  aabb[1] = igpx1;
            aabb[2] = igpy0;  aabb[3] = igpy1;
            aabb[4] = igpz0;  aabb[5] = igpz1;
          }
          else{
            if( igpx0 < aabb[0] ){ aabb[0] = igpx0; }
            if( igpx1 > aabb[1] ){ aabb[1] = igpx1; }
            if( igpy0 < aabb[2] ){ aabb[2] = igpy0; }
            if( igpy1 > aabb[3] ){ aabb[3] = igpy1; }
            if( igpz0 < aabb[4] ){ aabb[4] = igpz0; }
            if( igpz1 > aabb[5] ){ aabb[5] = igpz1; }
          }
        }
      }
    }
  }
  /*
  void Add(int ivx, int ivy, int ivz){
    int aabb[6];
    
  }
   */
  void Set(int ivx, int ivy, int ivz, int isVox){
    int igvx = ivx-iorgx;
    int igvy = ivy-iorgy;
    int igvz = ivz-iorgz;
    if( igvx<0 || igvx>ndivx ){ return; }
    if( igvy<0 || igvy>ndivy ){ return; }
    if( igvz<0 || igvz>ndivz ){ return; }
    const int ivoxel = igvx*(ndivy*ndivz)+igvy*ndivz+igvz;
    aIsVox[ivoxel] = isVox;
  }
  void GetQuad(std::vector<double>& aXYZ, std::vector<int>& aQuad){
    GetQuad_VoxelGrid(aXYZ, aQuad,
                      ndivx, ndivy, ndivz,
                      iorgx, iorgy, iorgz,
                      aIsVox);
  }
public:
  int ndivx, ndivy, ndivz;
  int iorgx, iorgy, iorgz;
  std::vector<int> aIsVox;
};

#endif
