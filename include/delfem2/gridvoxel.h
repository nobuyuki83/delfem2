/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_GRIDVOXEL_H
#define DFM2_GRIDVOXEL_H

#include <vector>

#include "delfem2/dfm2_inline.h"
#include "delfem2/mat4.h"
#include "delfem2/vec3.h"

namespace delfem2 {

int Adj_Grid(
    int igridvox, int iface,
    int ndivx, int ndivy, int ndivz);

bool IsInclude_AABB(
    const int aabb[8], int igvx, int igvy, int igvz);

void Add_AABB(
    int aabb[8],
    int ivx, int ivy, int ivz);

// -----------------------------------------------

void MeshQuad3D_VoxelGrid(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aQuad,
    unsigned int ndivx,
    unsigned int ndivy,
    unsigned int ndivz,
    const std::vector<int> &aIsVox);

void MeshHex3D_VoxelGrid(
    std::vector<double> &aXYZ,
    std::vector<int> &aQuad,
    unsigned int ndivx,
    unsigned int ndivy,
    unsigned int ndivz,
    const std::vector<int> &aIsVox);

void MeshTet3D_VoxelGrid(
    std::vector<double> &aXYZ,
    std::vector<int> &aTet,
    unsigned int ndivx,
    unsigned int ndivy,
    unsigned int ndivz,
    const std::vector<int> &aIsVox);

template<typename VAL>
class CGrid3 {
 public:
  CGrid3() {
    ndivx = ndivy = ndivz = 0;
    aVal.clear();
  }
  void Initialize(const unsigned int ndivx0,
                  const unsigned int ndivy0,
                  const unsigned int ndivz0,
                  const VAL v) {
    ndivx = ndivx0;
    ndivy = ndivy0;
    ndivz = ndivz0;
    const unsigned int nvoxel = ndivx * ndivy * ndivz;
    aVal.assign(nvoxel, v);
  }
  bool IsInclude(unsigned int ivx, unsigned int ivy, unsigned int ivz) {
    if (ivx >= ndivx) { return false; }
    if (ivy >= ndivy) { return false; }
    if (ivz >= ndivz) { return false; }
    return true;
  }
  void Set(int ivx, int ivy, int ivz, int isVox) {
    if (!this->IsInclude(ivx, ivy, ivz)) { return; }
    const int ivoxel = ivx * (ndivy * ndivz) + ivy * ndivz + ivz;
    aVal[ivoxel] = isVox;
  }
 public:
  unsigned int ndivx, ndivy, ndivz;
  std::vector<VAL> aVal;
  CMat4d am; // affine matrix
};

void Grid3Voxel_Dilation(CGrid3<int> &grid);

void Grid3Voxel_Erosion(CGrid3<int> &grid);

/**
 * @brief comute voxel geodesic distance from one seed voxel
 * @param aDist (out) geodesic distance at the center of the voxle
 * @param el (in) edge length of the voxel
 */
void VoxelGeodesic(
    std::vector<double> &aDist,
    const std::vector<std::pair<unsigned int, double> > &aIdvoxDist,
    double el,
    const CGrid3<int> &grid);

/**
 * @details this is the implementation of the following paper
 * Amanatides, John, and Andrew Woo. "A fast voxel traversal algorithm for ray tracing." In Eurographics, vol. 87, no. 3, pp. 3-10. 1987.
 * TODO: zero division might occur when (ps - pe ) is aligned to axis.
 */
void Intersection_VoxelGrid_LinSeg(
    std::vector<unsigned int> &aIndVox,
    const CGrid3<int> &grid,
    const CVec3d &ps,
    const CVec3d &pe);


/*
void Add(int ivx, int ivy, int ivz){
  if( this->IsInclude(ivx, ivy, ivz) ){
    Set(ivx,ivy,ivz,1);
  }
  else{ // size change
    int aabb[8]; this->AABB(aabb);
    Add_AABB(aabb,ivx,ivy,ivz);
    CVoxelGrid3D vg0 = (*this);
    this->Init_AABB(aabb);
    Set(ivx,ivy,ivz,1);
    for(int igvx0=0;igvx0<vg0.ndivx;++igvx0){
      for(int igvy0=0;igvy0<vg0.ndivy;++igvy0){
        for(int igvz0=0;igvz0<vg0.ndivz;++igvz0){
          const int ivoxel0 = igvx0*(vg0.ndivy*vg0.ndivz)+igvy0*vg0.ndivz+igvz0;
          if( vg0.aIsVox[ivoxel0] == 0 ) continue;
          int ivx1 = igvx0+vg0.iorgx;
          int ivy1 = igvy0+vg0.iorgy;
          int ivz1 = igvz0+vg0.iorgz;
          assert( IsInclude(ivx1,ivy1,ivz1) );
          Set(ivx1,ivy1,ivz1,1);
        }
      }
    }
  }
}
 */

}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/gridvoxel.cpp"
#endif

#endif
