#ifndef BVH_H
#define BVH_H

#include <stack>
#include <vector>
#include <set>
#include <assert.h>

class CNodeBVH
{
public:
  int iroot; // -1: root, else: parent
  int ichild[2]; // if ichild[1] == -1, then this is leaf and ichild[0] is the stored number
};

// make BVH topology
int MakeTreeTopologyBVH_TopDown
(std::vector<CNodeBVH>& aNodeBVH,
 const std::vector<int>& aElemSurInd,
 const std::vector<int>& aElemSur,
 const std::vector<double>& aElemCenter);

///////////////////////////////////////////////////////////////////////////////////////////
// following is the template function to define bounding volume

// build Bounding Box for AABB
template <typename T>
void BuildBoundingBoxesBVH
(int ibvh,
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aElem,
 int nnoel,
 const std::vector<CNodeBVH>& aNodeBVH,
 std::vector<T>& aBB)
{
  aBB.resize( aNodeBVH.size() );
  assert( ibvh < (int)aNodeBVH.size() );
  int ichild0 = aNodeBVH[ibvh].ichild[0];
  int ichild1 = aNodeBVH[ibvh].ichild[1];
  if( ichild1 == -1 ){ // leaf node
    const int ielem = ichild0;
    assert( ielem < (int)aElem.size()/nnoel );
    T& bb = aBB[ibvh];
    bb.is_active = false;
    for(int inoel=0;inoel<nnoel;++inoel){
      const int ino0 = aElem[ielem*nnoel+inoel];
      bb.AddPoint(aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2], delta*0.5);
    }
    return;
  }
  // branch node is the bounding volume of child nodes
  assert( aNodeBVH[ichild0].iroot == ibvh );
  assert( aNodeBVH[ichild1].iroot == ibvh );
  BuildBoundingBoxesBVH(ichild0,delta, aXYZ,aElem,nnoel,aNodeBVH,aBB);
  BuildBoundingBoxesBVH(ichild1,delta, aXYZ,aElem,nnoel,aNodeBVH,aBB);
  T& bb = aBB[ibvh];
  bb.is_active = false;
  bb  = aBB[ichild0];
  bb += aBB[ichild1];
  return;
}

// build Bounding Box for AABB
template <typename T>
void BuildBoundingBoxesBVH
(int ibvh,
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aElemInd,
 const std::vector<unsigned int>& aElem,
 const std::vector<CNodeBVH>& aNodeBVH,
 std::vector<T>& aBB)
{
  aBB.resize( aNodeBVH.size() );
  assert( ibvh < aNodeBVH.size() );
  int ichild0 = aNodeBVH[ibvh].ichild[0];
  int ichild1 = aNodeBVH[ibvh].ichild[1];
  if( ichild1 == -1 ){ // leaf node
    const int ielem = ichild0;
    assert( ielem < aElemInd.size()-1 );
    T& bb = aBB[ibvh];
    bb.is_active = false;
    for(int iip=aElemInd[ielem];iip<aElemInd[ielem+1];++iip){
      const int ino0 = aElem[iip];
      bb.AddPoint(aXYZ[ino0*3+0], aXYZ[ino0*3+1], aXYZ[ino0*3+2], delta*0.5);
    }
    return;
  }
  // branch node is the bounding volume of child nodes
  assert( aNodeBVH[ichild0].iroot == ibvh );
  assert( aNodeBVH[ichild1].iroot == ibvh );
  BuildBoundingBoxesBVH(ichild0,delta, aXYZ,aElemInd,aElem,aNodeBVH,aBB);
  BuildBoundingBoxesBVH(ichild1,delta, aXYZ,aElemInd,aElem,aNodeBVH,aBB);
  T& bb = aBB[ibvh];
  bb.is_active = false;
  bb  = aBB[ichild0];
  bb += aBB[ichild1];
  return;
}

template <typename T>
void BuildBoundingBoxesBVH_Dynamic
(int ibvh,
 double dt,
 const std::vector<double>& aXYZ,
 const std::vector<double>& aUVW,
 const std::vector<unsigned int>& aTri,
 const std::vector<CNodeBVH>& aNodeBVH,
 std::vector<T>& aBB)
{
  double eps = 1.0e-10;
  assert( ibvh < aNodeBVH.size() );
  int ichild0 = aNodeBVH[ibvh].ichild[0];
  int ichild1 = aNodeBVH[ibvh].ichild[1];
  if( ichild1 == -1 ){ // leaf
    const int itri = ichild0;
    assert( itri < aTri.size() );
    const int ino0 = aTri[itri*3+0];
    const int ino1 = aTri[itri*3+1];
    const int ino2 = aTri[itri*3+2];
    T& bb = aBB[ibvh];
    bb.is_active = false; // initialize
    bb.AddPoint(aXYZ[ino0*3+0],aXYZ[ino0*3+1],aXYZ[ino0*3+2], eps);
    bb.AddPoint(aXYZ[ino1*3+0],aXYZ[ino1*3+1],aXYZ[ino1*3+2], eps);
    bb.AddPoint(aXYZ[ino2*3+0],aXYZ[ino2*3+1],aXYZ[ino2*3+2], eps);
    bb.AddPoint(aXYZ[ino0*3+0]+dt*aUVW[ino0*3+0], aXYZ[ino0*3+1]+dt*aUVW[ino0*3+1], aXYZ[ino0*3+2]+dt*aUVW[ino0*3+2], eps);
    bb.AddPoint(aXYZ[ino1*3+0]+dt*aUVW[ino1*3+0], aXYZ[ino1*3+1]+dt*aUVW[ino1*3+1], aXYZ[ino1*3+2]+dt*aUVW[ino1*3+2], eps);
    bb.AddPoint(aXYZ[ino2*3+0]+dt*aUVW[ino2*3+0], aXYZ[ino2*3+1]+dt*aUVW[ino2*3+1], aXYZ[ino2*3+2]+dt*aUVW[ino2*3+2], eps);
    return;
  }
  // internal node,内部ノードは子ノードのBounding Volume
  assert( aNodeBVH[ichild0].iroot == ibvh );
  assert( aNodeBVH[ichild1].iroot == ibvh );
  BuildBoundingBoxesBVH_Dynamic(ichild0,dt, aXYZ,aUVW,aTri,aNodeBVH,aBB);
  BuildBoundingBoxesBVH_Dynamic(ichild1,dt, aXYZ,aUVW,aTri,aNodeBVH,aBB);
  T& bb = aBB[ibvh];
  bb.is_active = false;
  bb  = aBB[ichild0];
  bb += aBB[ichild1];
  return;
}

//////////////////////////////////////////////////////////////////////////////////


template <typename T>
void getBVH_IncludePoint
(std::vector<int>& aIndElem,
 /////
 double px, double py, double pz,
 int ibvh,
 const std::vector<CNodeBVH>& aBVH,
 const std::vector<T>& aBB)
{
  if( !aBB[ibvh].isInclude_Point(px,py,pz) ){ return; }
  const int ichild0 = aBVH[ibvh].ichild[0];
  const int ichild1 = aBVH[ibvh].ichild[1];
  if( ichild1 == -1 ){ // leaf
    aIndElem.push_back(ichild0);
    return;
  }
  /////
  getBVH_IncludePoint(aIndElem, px,py,pz, ichild0,  aBVH,aBB);
  getBVH_IncludePoint(aIndElem, px,py,pz, ichild1,  aBVH,aBB);
}

// potential maximum distance of the nearest point
template <typename T>
void getMinMaxDist_nearPoint
(double& max,
 /////
 double px, double py, double pz,
 int ibvh,
 const std::vector<CNodeBVH>& aBVH,
 const std::vector<T>& aBB)
{
  double min0, max0;
  aBB[ibvh].getRange_Point(min0,max0, px,py,pz);
  ////
  if( max>=0 && min0>max ){ return; }
  const int ichild0 = aBVH[ibvh].ichild[0];
  const int ichild1 = aBVH[ibvh].ichild[1];
  if( ichild1 == -1 ){ // leaf
    assert( aBB[ibvh].is_active );
    if( max<0 || max0<max ){ max=max0; }
    return;
  }
  /////
  getMinMaxDist_nearPoint(max, px,py,pz, ichild0,aBVH,aBB);
  getMinMaxDist_nearPoint(max, px,py,pz, ichild1,aBVH,aBB);
}

template <typename T>
void getBVH_NearPoint_OverlapDist
(std::vector<int>& aIndElem,
 /////
 double dist,
 double px, double py, double pz,
 int ibvh,
 const std::vector<CNodeBVH>& aBVH,
 const std::vector<T>& aBB)
{
  {
    double min0, max0;
    aBB[ibvh].getRange_Point(min0,max0, px,py,pz);
    if( dist<min0 ){ return; }
  }
  const int ichild0 = aBVH[ibvh].ichild[0];
  const int ichild1 = aBVH[ibvh].ichild[1];
  if( ichild1 == -1 ){ // leaf
    aIndElem.push_back(ichild0);
    return;
  }
  /////
  getBVH_NearPoint_OverlapDist(aIndElem, dist,px,py,pz, ichild0,aBVH,aBB);
  getBVH_NearPoint_OverlapDist(aIndElem, dist,px,py,pz, ichild1,aBVH,aBB);
}


#endif
