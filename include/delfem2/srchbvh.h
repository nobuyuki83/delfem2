/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @details this file is independent from any other delfem codes
 */

#ifndef DFM2_BVH_H
#define DFM2_BVH_H

#include <stack>
#include <vector>
#include <set>
#include <cassert>
#include <climits>
#include <iostream>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

/**
 * @details this class is used from cuda. So don't add default constructor. Cuda needs constructor with __device__  flag.
 */
class CNodeBVH2
{
public:
  unsigned int iparent; // UINT_MAX: root, else: parent
  unsigned int ichild[2]; // if ichild[1] == UINTMAX, then this is leaf and ichild[0] is the stored number
};

/**
 * @brief compute number of leading zeros
 * @function compute number of leading zeros
 * @param x input
 * @details clz(0) needs to be 32 to run BVH
 */
DFM2_INLINE unsigned int nbits_leading_zero(uint32_t x);

/**
 * @details make BVH topology in a top-down manner
 */
int BVHTopology_TopDown_MeshElem(
    std::vector<CNodeBVH2>& aNodeBVH,
    const unsigned int nfael,
    const std::vector<unsigned int>& aElSuEl,
    const std::vector<double>& aElemCenter);

/**
 * @details check if the leaf is visited once
 */
void Check_BVH(
    const std::vector<CNodeBVH2>& aNodeBVH,
    size_t nLeaf);


// --------------------------------------------
// code related to Morten code

/**
 * @returns return -1 if start == last
 * @details find split in BVH construction
 * https://devblogs.nvidia.com/thinking-parallel-part-iii-tree-construction-gpu/
 */
unsigned int MortonCode_FindSplit(
    const std::uint32_t* sortedMC,
    unsigned int start,
    unsigned int last);

/**
 * @details find range in parallel BVH construction
 * https://devblogs.nvidia.com/thinking-parallel-part-iii-tree-construction-gpu/
 */
std::pair<unsigned int,unsigned int> MortonCode_DeterminRange(
    const std::uint32_t* sortedMC,
    size_t nMC,
    unsigned int i);

/**
 * @brief compute morton code for 3d coordinates of a point. Each coordinate must be within the range of [0,1]
 * @details defined for "float" and "double"
 * https://devblogs.nvidia.com/thinking-parallel-part-iii-tree-construction-gpu/
 */
template <typename REAL>
DFM2_INLINE std::uint32_t MortonCode(REAL x, REAL y, REAL z);


/**
 * @details defined for "float" and "double"
 */
template <typename REAL>
void SortedMortenCode_Points3(
    std::vector<unsigned int> &aSortedId,
    std::vector<unsigned int> &aSortedMc,
    const std::vector<REAL> &aXYZ,
    const REAL min_xyz[3],
    const REAL max_xyz[3]);

void BVHTopology_Morton(
    std::vector<CNodeBVH2>& aNodeBVH,
    const std::vector<unsigned int>& aSortedId,
    const std::vector<std::uint32_t>& aSortedMc);

void Check_MortonCode_RangeSplit(
    const std::vector<std::uint32_t>& aSortedMc);

void Check_MortonCode_Sort(
    const std::vector<unsigned int>& aSortedId,
    const std::vector<std::uint32_t>& aSortedMc,
    const std::vector<double>& aXYZ,
    const double bbmin[3],
    const double bbmax[3]);

// above: code related to morton code
// -------------------------------------------------------------------
// below: template functions from here

/**
 * @brief build Bounding Box for AABB
 */
template <typename BBOX, typename LEAF_VOLUME_MAKER>
void BVH_BuildBVHGeometry(
    std::vector<BBOX>& aBB,
    unsigned int ibvh,
    const std::vector<CNodeBVH2>& aNodeBVH,
    const LEAF_VOLUME_MAKER& leafvolume);

template <typename BBOX, typename REAL>
class CLeafVolumeMaker_Mesh{
public:
  CLeafVolumeMaker_Mesh(
      REAL margin,
      const REAL* aXYZ,
      size_t nXYZ,
      const unsigned int* aElem,
      size_t nElem,
      unsigned int nnoel)
      : margin(margin),
        aXYZ(aXYZ), nXYZ(nXYZ),
        aElem(aElem), nElem(nElem), nnoel(nnoel)
  {}
  void SetVolume(BBOX& bb,
                 unsigned ielem) const {
    assert( ielem < nElem );
    bb.Set_Inactive();
    for(unsigned int inoel=0;inoel<nnoel;++inoel){
      const unsigned int ino0 = aElem[ielem*nnoel+inoel];
      bb.AddPoint(aXYZ+ino0*3, margin);
    }
  }
public:
  REAL margin;
  const REAL* aXYZ;
  size_t nXYZ;
  const unsigned int* aElem;
  size_t nElem;
  unsigned int nnoel;
};

template <typename BBOX, typename REAL>
class CLeafVolumeMaker_Point {
public:
  CLeafVolumeMaker_Point(
      const REAL* aXYZ,
      size_t nXYZ) : 
	  aXYZ(aXYZ), 
	  nXYZ(nXYZ) {
  }
  void SetVolume(BBOX& bb,
	  unsigned int ielem) const {
    assert( ielem < nXYZ );
    bb.AddPoint(aXYZ+ielem*3, 0.0);
    return;
  }
public:
  const REAL* aXYZ;
  size_t nXYZ;
};

template <typename BBOX, typename REAL>
class CLeafVolumeMaker_DynamicTriangle {
public:
  CLeafVolumeMaker_DynamicTriangle(
      double dt,
      const std::vector<double>& aXYZ,
      const std::vector<double>& aUVW,
      const std::vector<unsigned int>& aTri,
      double eps)
      : dt(dt), aXYZ(aXYZ), aUVW(aUVW), aTri(aTri), eps(eps) {}
  void SetVolume(BBOX& bb,
                 unsigned int itri) const {
    assert( itri < aTri.size() );
    const int ino0 = aTri[itri*3+0];
    const int ino1 = aTri[itri*3+1];
    const int ino2 = aTri[itri*3+2];
    // bb.bbmin[0] = +1;
    // bb.bbmax[0] = -1;
    bb.Set_Inactive();
    bb.AddPoint(aXYZ.data()+ino0*3, eps);
    bb.AddPoint(aXYZ.data()+ino1*3, eps);
    bb.AddPoint(aXYZ.data()+ino2*3, eps);
    const double p0[3] = {aXYZ[ino0*3+0]+dt*aUVW[ino0*3+0], aXYZ[ino0*3+1]+dt*aUVW[ino0*3+1], aXYZ[ino0*3+2]+dt*aUVW[ino0*3+2]};
    const double p1[3] = {aXYZ[ino1*3+0]+dt*aUVW[ino1*3+0], aXYZ[ino1*3+1]+dt*aUVW[ino1*3+1], aXYZ[ino1*3+2]+dt*aUVW[ino1*3+2]};
    const double p2[3] = {aXYZ[ino2*3+0]+dt*aUVW[ino2*3+0], aXYZ[ino2*3+1]+dt*aUVW[ino2*3+1], aXYZ[ino2*3+2]+dt*aUVW[ino2*3+2]};
    bb.AddPoint(p0, eps);
    bb.AddPoint(p1, eps);
    bb.AddPoint(p2, eps);
  }
public:
  double dt;
  const std::vector<double>& aXYZ;
  const std::vector<double>& aUVW;
  const std::vector<unsigned int>& aTri;
  double eps;
};

template <typename BBOX, typename PREDICATE>
void BVH_GetIndElem_Predicate(
    std::vector<unsigned int>& aIndElem,
    //
    PREDICATE pred,
    unsigned int ibvh,
    const std::vector<delfem2::CNodeBVH2>& aBVH,
    const std::vector<BBOX>& aBB)
{
  assert( ibvh < aBVH.size() );
  const bool is_intersect = pred.IsTrue(ibvh,aBB);
  if( !is_intersect ) return;
  //
  const unsigned int ichild0 = aBVH[ibvh].ichild[0];
  const unsigned int ichild1 = aBVH[ibvh].ichild[1];
  if( ichild1 == UINT_MAX ){ // leaf
    assert( ichild0 < aBB.size() );
    aIndElem.push_back(ichild0);
    return;
  }
  BVH_GetIndElem_Predicate(aIndElem, pred, ichild0,aBVH,aBB);
  BVH_GetIndElem_Predicate(aIndElem, pred, ichild1,aBVH,aBB);
}

template <class BV,typename REAL>
class CIsBV_IntersectLine
{
public:
  CIsBV_IntersectLine(
      const REAL src_[3],
      const REAL dir_[3]) :
      src{src_[0],src_[1],src_[2]},
      dir{dir_[0],dir_[1],dir_[2]}
  {}
  bool IsTrue(unsigned int ibv,
              const std::vector<BV>& aBV)
  {
    return aBV[ibv].IsIntersectLine(src,dir);
  }
public:
  const REAL src[3];
  const REAL dir[3];
};


template <typename BV>
class CIsBV_IncludePoint
{
public:
  CIsBV_IncludePoint(const double pos_[3]) :
    pos{pos_[0],pos_[1],pos_[2]} {}
  bool IsTrue(unsigned int ibvh, const std::vector<BV>& aBB) {
    return aBB[ibvh].isInclude_Point(pos[0],pos[1],pos[2]);
  }
public:
  const double pos[3];
};

template <typename BV>
class CIsBV_IntersectRay
{
public:
  CIsBV_IntersectRay(const double src_[3], const double dir_[3]) :
      src{src_[0],src_[1],src_[2]},
      dir{dir_[0],dir_[1],dir_[2]} {}
  bool IsTrue(unsigned int ibvh, const std::vector<BV>& aBB){
    return aBB[ibvh].IsIntersectRay(src,dir);
  }
public:
  const double src[3];
  const double dir[3];
};


template <typename BV>
class CIsBV_InsideRange
{
public:
  CIsBV_InsideRange(
      const double pos_[3],
      double min_,
      double max_) :
      pos{pos_[0],pos_[1],pos_[2]},
      min(min_),
      max(max_) {}
  bool IsTrue(unsigned int ibvh, const std::vector<BV>& aBB){
    double min0, max0;
    aBB[ibvh].Range_DistToPoint(min0,max0, pos[0],pos[1],pos[2]);
    if( max0<min || min0>max ){ return false; }
    return true;
  }
public:
  double pos[3];
  double min, max;
};

/**
 * @brief potential maximum distance of the nearest point
 * @details set some value with min > max for input  e.g,. min=+1, max=-1
 */
template <typename BBOX>
void BVH_Range_DistToNearestPoint(
    double& min, double& max,
    //
    const double p[3],
    unsigned int ibvh,
    const std::vector<delfem2::CNodeBVH2>& aBVH,
    const std::vector<BBOX>& aBB);

/**
 * @brief find nearest point
 * @param min
 * @param aBB a bounding box of nodes
 */
template <typename BBOX, typename REAL>
void BVH_IndPoint_NearestPoint(
    unsigned int& ip,
    REAL& dist_cur,
    //
    const REAL p[3],
    unsigned int ibvh,
    const std::vector<delfem2::CNodeBVH2>& aBVH,
    const std::vector<BBOX>& aBB);

} // end namespace delfem2


// -------------------------------------------------------------------------
// below: building geometry


/**
 * @brief build Bounding Box for BVH
 */
template <typename BBOX, typename LEAF_VOLUME_MAKER>
void delfem2::BVH_BuildBVHGeometry(
    std::vector<BBOX>& aBB,
    unsigned int ibvh,
    const std::vector<delfem2::CNodeBVH2>& aNodeBVH,
    const LEAF_VOLUME_MAKER& lvm)
{
  aBB.resize( aNodeBVH.size() );
  assert( ibvh < aNodeBVH.size() );
  const unsigned int ichild0 = aNodeBVH[ibvh].ichild[0];
  const unsigned int ichild1 = aNodeBVH[ibvh].ichild[1];
  if( ichild1 == UINT_MAX ){ // leaf node    
    const unsigned int ielem = ichild0;
    lvm.SetVolume(aBB[ibvh],
                  ielem);
    return;
  }
  // branch node is the bounding volume of child nodes
  assert( aNodeBVH[ichild0].iparent == ibvh );
  assert( aNodeBVH[ichild1].iparent == ibvh );
  BVH_BuildBVHGeometry(aBB, ichild0,aNodeBVH, lvm);
  BVH_BuildBVHGeometry(aBB, ichild1,aNodeBVH, lvm);
  BBOX& bb = aBB[ibvh];
  bb  = aBB[ichild0];
  bb += aBB[ichild1];
  return;
}

// ------------------------------------------------------------------------



/**
 * @brief potential maximum distance of the nearest point
 * @details set some value with min > max for input  e.g,. min=+1, max=-1
 */
template <typename BBOX>
void delfem2::BVH_Range_DistToNearestPoint(
    double& min, double& max,
    //
    const double p[3],
    unsigned int ibvh,
    const std::vector<delfem2::CNodeBVH2>& aBVH,
    const std::vector<BBOX>& aBB)
{
  double min0=+1.0, max0=-1.0;
  aBB[ibvh].Range_DistToPoint(min0,max0, p[0],p[1],p[2]);
  if( max0 < min0 ){ return; } // ibvh is a inactive bvh the children should be inactive too
  //
  if( max>=min && min0>max ){ return; } // current range [min,max] is valid and nearer than [min0,min0].
  const unsigned int ichild0 = aBVH[ibvh].ichild[0];
  const unsigned int ichild1 = aBVH[ibvh].ichild[1];
  if( ichild1 == UINT_MAX ){ // leaf
    if( max<min ){ // current range is inactive
      max = max0;
      min = min0;
      return;
    }
    if( max0 < max ){ max = max0; }
    if( min0 < min ){ min = min0; }
    return;
  }
  //
  BVH_Range_DistToNearestPoint(min,max, p, ichild0,aBVH,aBB);
  BVH_Range_DistToNearestPoint(min,max, p, ichild1,aBVH,aBB);
}

/**
 * @brief index of the point nearest to the given point
 * @details the cur_dist should be input as a negative value (e.g., cur_dist=-1)
 */
template <typename BBOX, typename REAL>
void delfem2::BVH_IndPoint_NearestPoint(
    unsigned int& ip,
    REAL& cur_dist,
    //
    const REAL p[3],
    unsigned int ibvh,
    const std::vector<delfem2::CNodeBVH2>& aBVH,
    const std::vector<BBOX>& aBB)
{
  assert( aBVH.size() == aBB.size() );
  REAL min0=+1.0, max0=-1.0;
  aBB[ibvh].Range_DistToPoint(min0,max0, p[0],p[1],p[2]);
  if( max0 < min0 ){ return; } // ibvh is a inactive bvh the children should be inactive too
  if( cur_dist > 0 && min0>cur_dist ){ return; } // current range [min,max] is valid and nearer than [min0,min0].
  const unsigned int ichild0 = aBVH[ibvh].ichild[0];
  const unsigned int ichild1 = aBVH[ibvh].ichild[1];
  if( ichild1 == UINT_MAX ){ // leaf
    assert( min0 == max0 ); // because this is point
    if( cur_dist < 0 || max0 < cur_dist ){ // current range is inactive
      cur_dist = max0;
      ip = ichild0;
    }
    return;
  }
  //
  BVH_IndPoint_NearestPoint(ip,cur_dist, p, ichild0,aBVH,aBB);
  BVH_IndPoint_NearestPoint(ip,cur_dist, p, ichild1,aBVH,aBB);
}







#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/srchbvh.cpp"
#endif


#endif
