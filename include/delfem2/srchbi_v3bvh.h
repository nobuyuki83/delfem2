/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @details bilatial search (elemnet-to-element) using bvh
 */

#ifndef DFM2_SRCHBI_V3BVH_H
#define DFM2_SRCHBI_V3BVH_H

#include "delfem2/srchbvh.h"
#include "delfem2/geoproximity3_v3.h"
#include "delfem2/geosolidelm_v3.h"
#include "delfem2/vec3.h"
#include <stdio.h>


namespace delfem2 {

template <typename BBOX>
bool IsContact_FV_Proximity(
    int ino0, int ino1, int ino2, int ino3,
    const CVec3d& p0, const CVec3d& p1, const CVec3d& p2, const CVec3d& p3,
    const BBOX& bb,
    const double delta);

/**
 * @brief check two edge elements collide or not in the continuous time
 * @tparam BBOX bounding box class must equipp with AddPoint and IsIntersect
 * @return whether this element collide or not
 */
template <typename BBOX>
bool IsContact_EE_CCD(
    int ino0,         int ino1,         int jno0,         int jno1,
    const CVec3d& p0s, const CVec3d& p1s, const CVec3d& q0s, const CVec3d& q1s,
    const CVec3d& p0e, const CVec3d& p1e, const CVec3d& q0e, const CVec3d& q1e);

template <typename BBOX>
bool IsContact_FV_CCD(
    int ino0,        int ino1,        int ino2,        int ino3,
    const CVec3d& p0, const CVec3d& p1, const CVec3d& p2, const CVec3d& p3,
    const CVec3d& q0, const CVec3d& q1, const CVec3d& q2, const CVec3d& q3,
    const BBOX& bb);

// -------------
class CContactElement;

template <typename BBOX>
void GetContactElement_Proximity(
    std::set<CContactElement>& aContactElem,
    // ----------
    double delta,
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri,
    int ibvh0, int ibvh1,
    const std::vector<CNodeBVH2>& aBVH,
    const std::vector<BBOX>& aBB);

template <typename BBOX>
void GetContactElement_Proximity(
    std::set<CContactElement>& aContactElem,
    // ----------
    double delta,
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri,
    int ibvh,
    const std::vector<CNodeBVH2>& aBVH,
    const std::vector<BBOX>& aBB);

template <typename BBOX>
void GetContactElement_CCD(
    std::set<CContactElement>& aContactElem,
    // ------------
    double dt,
    double delta,
    const std::vector<double>& aXYZ,
    const std::vector<double>& aUVW,
    const std::vector<unsigned int>& aTri,
    int ibvh,
    const std::vector<CNodeBVH2>& aBVH,
    const std::vector<BBOX>& aBB);

template <typename BBOX>
void GetContactElement_CCD(
    std::set<CContactElement>& aContactElem,
    // --------------
    double dt,
    double delta,
    const std::vector<double>& aXYZ,
    const std::vector<double>& aUVW,
    const std::vector<unsigned int>& aTri,
    int ibvh0, int ibvh1,
    const std::vector<CNodeBVH2>& aBVH,
    const std::vector<BBOX>& aBB);

// ------------------------
template <typename REAL>
class CIntersectTriPair;

template <typename BBOX>
void GetIntersectTriPairs(
    std::vector<CIntersectTriPair<double>>& aIntersectTriPair,
    // --------------
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri,
    int ibvh,
    const std::vector<CNodeBVH2>& aBVH,
    const std::vector<BBOX>& aBB);

template <typename BBOX>
void GetIntersectTriPairs(
    std::vector<CIntersectTriPair<double>>& aIntersectTriPair,
    // ------------------
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri,
    int ibvh0, int ibvh1,
    const std::vector<CNodeBVH2>& aBVH,
    const std::vector<BBOX>& aBB);

} // end namespace delfem2


// --------------------------------------------

template <typename BBOX>
bool delfem2::IsContact_FV_Proximity
(int ino0, int ino1, int ino2, int ino3,
 const CVec3d& p0, const CVec3d& p1, const CVec3d& p2, const CVec3d& p3,
 const BBOX& bb,
 const double delta)
{
  if( ino3 == ino0 || ino3 == ino1 || ino3 == ino2 ){ return false; }
  if( !bb.isInclude_Point(p3.x,p3.y,p3.z) ) return false;
  double height = fabs( Height(p0,p1,p2,p3) );
  if( height > delta ) return false;
  double w0,w1;
  const double dist = DistanceFaceVertex(p0,p1,p2,p3,w0,w1);
  const double w2 = 1-w0-w1;
  if( dist > delta ) return false;
  double mgn = ( Distance(p0, p1) + Distance(p1, p2) + Distance(p2, p3) ) / 3.0;
  mgn = 0;
  if( w0 < -mgn || w0 > 1+mgn ) return false;
  if( w1 < -mgn || w1 > 1+mgn ) return false;
  if( w2 < -mgn || w2 > 1+mgn ) return false;
  return true;
}

// check if two edge elements collide or not
template <typename BBOX>
bool delfem2::IsContact_EE_CCD
(int ino0,         int ino1,         int jno0,         int jno1,
 const CVec3d& p0s, const CVec3d& p1s, const CVec3d& q0s, const CVec3d& q1s,
 const CVec3d& p0e, const CVec3d& p1e, const CVec3d& q0e, const CVec3d& q1e)
{
  double eps = 1.0e-10;
  if( ino0 == jno0 || ino0 == jno1 || ino1 == jno0 || ino1 == jno1 ) return false;
  BBOX bbq;
  bbq.AddPoint(q0s.data(), eps);
  bbq.AddPoint(q1s.data(), eps);
  bbq.AddPoint(q0e.data(), eps);
  bbq.AddPoint(q1e.data(), eps);
  //
  BBOX bbp;
  bbp.AddPoint(p0s.data(), eps);
  bbp.AddPoint(p1s.data(), eps);
  bbp.AddPoint(p0e.data(), eps);
  bbp.AddPoint(p1e.data(), eps);
  if( !bbp.IsIntersect(bbq) ) return false;
  //
  double t;
  {
    const bool res = FindCoplanerInterp(t,
                                        p0s,p1s,q0s,q1s, p0e,p1e,q0e,q1e);
    if( !res ) return false;
    assert( t >= 0 && t <= 1 );
  }
  CVec3d p0m = (1-t)*p0s + t*p0e;
  CVec3d p1m = (1-t)*p1s + t*p1e;
  CVec3d q0m = (1-t)*q0s + t*q0e;
  CVec3d q1m = (1-t)*q1s + t*q1e;
  double w0,w1;
  double dist = DistanceEdgeEdge(p0m, p1m, q0m, q1m, w0,w1);
  if( w0 < 0 || w0 > 1 ) return false;
  if( w1 < 0 || w1 > 1 ) return false;
  if( dist > 1.0e-2 ) return false;
  return true;
}


// ---------------------

namespace delfem2 {

// smallest element of contact (vertex-face or edge-edge)
class CContactElement
{
public:
  CContactElement(bool is_fv,int j0, int j1, int j2, int j3)
  {
    this->is_fv = is_fv;
    if( is_fv ){
      this->is_fv = true;  //  re-order such that coputation is fast.
      if(      j0 < j1 && j0 < j2 && j1 < j2 ){ ino0=j0;  ino1=j1;  ino2=j2;  ino3=j3; }
      else if( j0 < j1 && j0 < j2 && j2 < j1 ){ ino0=j0;  ino1=j2;  ino2=j1;  ino3=j3; }
      else if( j1 < j0 && j1 < j2 && j0 < j2 ){ ino0=j1;  ino1=j0;  ino2=j2;  ino3=j3; }
      else if( j1 < j0 && j1 < j2 && j2 < j0 ){ ino0=j1;  ino1=j2;  ino2=j0;  ino3=j3; }
      else if( j2 < j0 && j2 < j1 && j0 < j1 ){ ino0=j2;  ino1=j0;  ino2=j1;  ino3=j3; }
      else if( j2 < j0 && j2 < j1 && j1 < j0 ){ ino0=j2;  ino1=j1;  ino2=j0;  ino3=j3; }
      else { assert(0); }
    }
    else{
      this->is_fv = false; //  re-order such that coputation is fast.
      if(      j0 < j1 && j0 < j2 && j0 < j3 && j2 < j3 ){ ino0=j0;  ino1=j1;  ino2=j2;  ino3=j3; }
      else if( j0 < j1 && j0 < j2 && j0 < j3 && j3 < j2 ){ ino0=j0;  ino1=j1;  ino2=j3;  ino3=j2; }
      else if( j1 < j0 && j1 < j2 && j1 < j3 && j2 < j3 ){ ino0=j1;  ino1=j0;  ino2=j2;  ino3=j3; }
      else if( j1 < j0 && j1 < j2 && j1 < j3 && j3 < j2 ){ ino0=j1;  ino1=j0;  ino2=j3;  ino3=j2; }
      else if( j2 < j0 && j2 < j1 && j2 < j3 && j0 < j1 ){ ino0=j2;  ino1=j3;  ino2=j0;  ino3=j1; }
      else if( j2 < j0 && j2 < j1 && j2 < j3 && j1 < j0 ){ ino0=j2;  ino1=j3;  ino2=j1;  ino3=j0; }
      else if( j3 < j0 && j3 < j1 && j3 < j2 && j0 < j1 ){ ino0=j3;  ino1=j2;  ino2=j0;  ino3=j1; }
      else if( j3 < j0 && j3 < j1 && j3 < j2 && j1 < j0 ){ ino0=j3;  ino1=j2;  ino2=j1;  ino3=j0; }
      else { assert(0); }
    }
  }
  bool operator < (const CContactElement& p2) const
  {
    if( ino0 != p2.ino0 ){ return ino0 < p2.ino0; }
    if( ino1 != p2.ino1 ){ return ino1 < p2.ino1; }
    if( ino2 != p2.ino2 ){ return ino2 < p2.ino2; }
    return ino3 < p2.ino3;
  }
public:
  bool is_fv; // true: ee contact, false: vf contact
  int ino0, ino1, ino2, ino3; // four points in the contact
};
  
}

template <typename BBOX>
void delfem2::GetContactElement_Proximity
(std::set<CContactElement>& aContactElem,
 // -------------------
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 int ibvh0, int ibvh1,
 const std::vector<delfem2::CNodeBVH2>& aBVH,
 const std::vector<BBOX>& aBB)
{
  assert( ibvh0 < (int)aBB.size() );
  assert( ibvh1 < (int)aBB.size() );
  if( !aBB[ibvh0].IsIntersect(aBB[ibvh1]) ) return;
  const int ichild0_0 = aBVH[ibvh0].ichild[0];
  const int ichild0_1 = aBVH[ibvh0].ichild[1];
  const int ichild1_0 = aBVH[ibvh1].ichild[0];
  const int ichild1_1 = aBVH[ibvh1].ichild[1];
  const bool is_leaf0 = (ichild0_1 == -1);
  const bool is_leaf1 = (ichild1_1 == -1);
  if(      !is_leaf0 && !is_leaf1 ){
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_0,ichild1_0, aBVH,aBB);
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_1,ichild1_0, aBVH,aBB);
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_0,ichild1_1, aBVH,aBB);
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_1,ichild1_1, aBVH,aBB);
  }
  else if( !is_leaf0 &&  is_leaf1 ){
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_0,ibvh1,aBVH,aBB);
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0_1,ibvh1,aBVH,aBB);
  }
  else if(  is_leaf0 && !is_leaf1 ){
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ibvh0,ichild1_0,aBVH,aBB);
    GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ibvh0,ichild1_1,aBVH,aBB);
  }
  else if(  is_leaf0 &&  is_leaf1 ){
    const int itri = ichild0_0;
    const int jtri = ichild1_0;
    const int in0 = aTri[itri*3+0];
    const int in1 = aTri[itri*3+1];
    const int in2 = aTri[itri*3+2];
    const int jn0 = aTri[jtri*3+0];
    const int jn1 = aTri[jtri*3+1];
    const int jn2 = aTri[jtri*3+2];
    const CVec3d p0(aXYZ[in0*3+0], aXYZ[in0*3+1], aXYZ[in0*3+2]);
    const CVec3d p1(aXYZ[in1*3+0], aXYZ[in1*3+1], aXYZ[in1*3+2]);
    const CVec3d p2(aXYZ[in2*3+0], aXYZ[in2*3+1], aXYZ[in2*3+2]);
    const CVec3d q0(aXYZ[jn0*3+0], aXYZ[jn0*3+1], aXYZ[jn0*3+2]);
    const CVec3d q1(aXYZ[jn1*3+0], aXYZ[jn1*3+1], aXYZ[jn1*3+2]);
    const CVec3d q2(aXYZ[jn2*3+0], aXYZ[jn2*3+1], aXYZ[jn2*3+2]);
    if( IsContact_FV_Proximity(   in0,in1,in2,jn0, p0,p1,p2,q0, aBB[ichild0_0], delta) ){
      aContactElem.insert( CContactElement(true,    in0,in1,in2,jn0) );
    }
    if( IsContact_FV_Proximity(   in0,in1,in2,jn1, p0,p1,p2,q1, aBB[ichild0_0], delta) ){
      aContactElem.insert( CContactElement(true,    in0,in1,in2,jn1) );
    }
    if( IsContact_FV_Proximity(   in0,in1,in2,jn2, p0,p1,p2,q2, aBB[ichild0_0], delta) ){
      aContactElem.insert( CContactElement(true,    in0,in1,in2,jn2) );
    }
    if( IsContact_FV_Proximity(   jn0,jn1,jn2,in0, q0,q1,q2,p0, aBB[ichild1_0], delta) ){
      aContactElem.insert( CContactElement(true,    jn0,jn1,jn2,in0) );
    }
    if( IsContact_FV_Proximity(   jn0,jn1,jn2,in1, q0,q1,q2,p1, aBB[ichild1_0], delta) ){
      aContactElem.insert( CContactElement(true,    jn0,jn1,jn2,in1) );
    }
    if( IsContact_FV_Proximity(   jn0,jn1,jn2,in2, q0,q1,q2,p2, aBB[ichild1_0], delta) ){
      aContactElem.insert( CContactElement(true,    jn0,jn1,jn2,in2) );
    }
    ////
    if( IsContact_EE_Proximity(      in0,in1,jn0,jn1, p0,p1,q0,q1, delta) ){
      aContactElem.insert( CContactElement(false,    in0,in1,jn0,jn1) );
    }
    if( IsContact_EE_Proximity(      in0,in1,jn1,jn2, p0,p1,q1,q2, delta) ){
      aContactElem.insert( CContactElement(false,    in0,in1,jn1,jn2) );
    }
    if( IsContact_EE_Proximity(      in0,in1,jn2,jn0, p0,p1,q2,q0, delta) ){
      aContactElem.insert( CContactElement(false,    in0,in1,jn2,jn0) );
    }
    if( IsContact_EE_Proximity(      in1,in2,jn0,jn1, p1,p2,q0,q1, delta) ){
      aContactElem.insert( CContactElement(false,    in1,in2,jn0,jn1) );
    }
    if( IsContact_EE_Proximity(      in1,in2,jn1,jn2, p1,p2,q1,q2, delta) ){
      aContactElem.insert( CContactElement(false,    in1,in2,jn1,jn2) );
    }
    if( IsContact_EE_Proximity(      in1,in2,jn2,jn0, p1,p2,q2,q0, delta) ){
      aContactElem.insert( CContactElement(false,    in1,in2,jn2,jn0) );
    }
    if( IsContact_EE_Proximity(      in2,in0,jn0,jn1, p2,p0,q0,q1, delta) ){
      aContactElem.insert( CContactElement(false,    in2,in0,jn0,jn1) );
    }
    if( IsContact_EE_Proximity(      in2,in0,jn1,jn2, p2,p0,q1,q2, delta) ){
      aContactElem.insert( CContactElement(false,    in2,in0,jn1,jn2) );
    }
    if( IsContact_EE_Proximity(      in2,in0,jn2,jn0, p2,p0,q2,q0, delta) ){
      aContactElem.insert( CContactElement(false,    in2,in0,jn2,jn0) );
    }
  }
}

template <typename BBOX>
void delfem2::GetContactElement_Proximity
(std::set<delfem2::CContactElement>& aContactElem,
 ////
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 int ibvh,
 const std::vector<delfem2::CNodeBVH2>& aBVH,
 const std::vector<BBOX>& aBB)
{
  const int ichild0 = aBVH[ibvh].ichild[0];
  const int ichild1 = aBVH[ibvh].ichild[1];
  const bool is_leaf = (ichild1 == -1);
  if( is_leaf ) return;
  GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0,ichild1,aBVH,aBB);
  GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild0,        aBVH,aBB);
  GetContactElement_Proximity(aContactElem, delta,aXYZ,aTri, ichild1,        aBVH,aBB);
}


// CCDのFVで接触する要素を検出
template <typename T>
bool delfem2::IsContact_FV_CCD
(int ino0,        int ino1,        int ino2,        int ino3,
 const CVec3d& p0, const CVec3d& p1, const CVec3d& p2, const CVec3d& p3,
 const CVec3d& q0, const CVec3d& q1, const CVec3d& q2, const CVec3d& q3,
 const T& bb)
{
  double eps = 1.0e-10;
  if( ino3 == ino0 || ino3 == ino1 || ino3 == ino2 ){ return false; }
  { // culling
    T bbp;
    bbp.AddPoint(p3.data(), eps);
    bbp.AddPoint(q3.data(), eps);
    if( !bb.IsIntersect(bbp) ) return false;
  }
  return IsContact_FV_CCD2(ino0, ino1,ino2,ino3, p0,p1,p2,p3, q0, q1,q2,q3);
}

// detect contact element with Continous Collision Detection (CCD)
template <typename BBOX>
void delfem2::GetContactElement_CCD
(std::set<CContactElement>& aContactElem,
 // --------------------
 double dt,
 double delta,
 const std::vector<double>& aXYZ,
 const std::vector<double>& aUVW,
 const std::vector<unsigned int>& aTri,
 int ibvh0, int ibvh1,
 const std::vector<CNodeBVH2>& aBVH,
 const std::vector<BBOX>& aBB)
{
  assert( ibvh0 < (int)aBB.size() );
  assert( ibvh1 < (int)aBB.size() );
  if( !aBB[ibvh0].IsIntersect(aBB[ibvh1]) ) return;
  const int ichild0_0 = aBVH[ibvh0].ichild[0];
  const int ichild0_1 = aBVH[ibvh0].ichild[1];
  const int ichild1_0 = aBVH[ibvh1].ichild[0];
  const int ichild1_1 = aBVH[ibvh1].ichild[1];
  const bool is_leaf0 = (ichild0_1 == -1);
  const bool is_leaf1 = (ichild1_1 == -1);
  if(      !is_leaf0 && !is_leaf1 ){
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_0,ichild1_0, aBVH,aBB);
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_1,ichild1_0, aBVH,aBB);
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_0,ichild1_1, aBVH,aBB);
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_1,ichild1_1, aBVH,aBB);
  }
  else if( !is_leaf0 &&  is_leaf1 ){
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_0,ibvh1,     aBVH,aBB);
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0_1,ibvh1,     aBVH,aBB);
  }
  else if(  is_leaf0 && !is_leaf1 ){
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ibvh0,    ichild1_0, aBVH,aBB);
    GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ibvh0,    ichild1_1, aBVH,aBB);
  }
  else if(  is_leaf0 &&  is_leaf1 ){
    const int itri = ichild0_0;
    const int jtri = ichild1_0;
    int in0 = aTri[itri*3+0];
    int in1 = aTri[itri*3+1];
    int in2 = aTri[itri*3+2];
    int jn0 = aTri[jtri*3+0];
    int jn1 = aTri[jtri*3+1];
    int jn2 = aTri[jtri*3+2];
    const CVec3d p0s(aXYZ[in0*3+0],                  aXYZ[in0*3+1],                  aXYZ[in0*3+2]);
    const CVec3d p1s(aXYZ[in1*3+0],                  aXYZ[in1*3+1],                  aXYZ[in1*3+2]);
    const CVec3d p2s(aXYZ[in2*3+0],                  aXYZ[in2*3+1],                  aXYZ[in2*3+2]);
    const CVec3d q0s(aXYZ[jn0*3+0],                  aXYZ[jn0*3+1],                  aXYZ[jn0*3+2]);
    const CVec3d q1s(aXYZ[jn1*3+0],                  aXYZ[jn1*3+1],                  aXYZ[jn1*3+2]);
    const CVec3d q2s(aXYZ[jn2*3+0],                  aXYZ[jn2*3+1],                  aXYZ[jn2*3+2]);
    const CVec3d p0e(aXYZ[in0*3+0]+dt*aUVW[in0*3+0], aXYZ[in0*3+1]+dt*aUVW[in0*3+1], aXYZ[in0*3+2]+dt*aUVW[in0*3+2]);
    const CVec3d p1e(aXYZ[in1*3+0]+dt*aUVW[in1*3+0], aXYZ[in1*3+1]+dt*aUVW[in1*3+1], aXYZ[in1*3+2]+dt*aUVW[in1*3+2]);
    const CVec3d p2e(aXYZ[in2*3+0]+dt*aUVW[in2*3+0], aXYZ[in2*3+1]+dt*aUVW[in2*3+1], aXYZ[in2*3+2]+dt*aUVW[in2*3+2]);
    const CVec3d q0e(aXYZ[jn0*3+0]+dt*aUVW[jn0*3+0], aXYZ[jn0*3+1]+dt*aUVW[jn0*3+1], aXYZ[jn0*3+2]+dt*aUVW[jn0*3+2]);
    const CVec3d q1e(aXYZ[jn1*3+0]+dt*aUVW[jn1*3+0], aXYZ[jn1*3+1]+dt*aUVW[jn1*3+1], aXYZ[jn1*3+2]+dt*aUVW[jn1*3+2]);
    const CVec3d q2e(aXYZ[jn2*3+0]+dt*aUVW[jn2*3+0], aXYZ[jn2*3+1]+dt*aUVW[jn2*3+1], aXYZ[jn2*3+2]+dt*aUVW[jn2*3+2]);
    
    if( IsContact_FV_CCD(      in0,in1,in2,jn0, p0s,p1s,p2s,q0s, p0e,p1e,p2e,q0e, aBB[ibvh0]) ){
      aContactElem.insert( CContactElement(true, in0,in1,in2,jn0) );
    }
    if( IsContact_FV_CCD(      in0,in1,in2,jn1, p0s,p1s,p2s,q1s, p0e,p1e,p2e,q1e, aBB[ibvh0]) ){
      aContactElem.insert( CContactElement(true, in0,in1,in2,jn1) );
    }
    if( IsContact_FV_CCD(      in0,in1,in2,jn2, p0s,p1s,p2s,q2s, p0e,p1e,p2e,q2e, aBB[ibvh0]) ){
      aContactElem.insert( CContactElement(true, in0,in1,in2,jn2) );
    }
    if( IsContact_FV_CCD(      jn0,jn1,jn2,in0, q0s,q1s,q2s,p0s, q0e,q1e,q2e,p0e, aBB[ibvh1]) ){
      aContactElem.insert( CContactElement(true, jn0,jn1,jn2,in0) );
    }
    if( IsContact_FV_CCD(      jn0,jn1,jn2,in1, q0s,q1s,q2s,p1s, q0e,q1e,q2e,p1e, aBB[ibvh1]) ){
      aContactElem.insert( CContactElement(true, jn0,jn1,jn2,in1) );
    }
    if( IsContact_FV_CCD(      jn0,jn1,jn2,in2, q0s,q1s,q2s,p2s, q0e,q1e,q2e,p2e, aBB[ibvh1]) ){
      aContactElem.insert( CContactElement(true, jn0,jn1,jn2,in2) );
    }
    ////
    if( IsContact_EE_CCD<BBOX>(          in0,in1,jn0,jn1, p0s,p1s,q0s,q1s,  p0e,p1e,q0e,q1e) ){
      aContactElem.insert( CContactElement(false,  in0,in1,jn0,jn1) );
    }
    if( IsContact_EE_CCD<BBOX>(          in0,in1,jn1,jn2, p0s,p1s,q1s,q2s,  p0e,p1e,q1e,q2e) ){
      aContactElem.insert( CContactElement(false,  in0,in1,jn1,jn2) );
    }
    if( IsContact_EE_CCD<BBOX>(          in0,in1,jn2,jn0, p0s,p1s,q2s,q0s,  p0e,p1e,q2e,q0e) ){
      aContactElem.insert( CContactElement(false,  in0,in1,jn2,jn0) );
    }
    if( IsContact_EE_CCD<BBOX>(          in1,in2,jn0,jn1, p1s,p2s,q0s,q1s,  p1e,p2e,q0e,q1e) ){
      aContactElem.insert( CContactElement(false,  in1,in2,jn0,jn1) );
    }
    if( IsContact_EE_CCD<BBOX>(          in1,in2,jn1,jn2, p1s,p2s,q1s,q2s,  p1e,p2e,q1e,q2e) ){
      aContactElem.insert( CContactElement(false,  in1,in2,jn1,jn2) );
    }
    if( IsContact_EE_CCD<BBOX>(          in1,in2,jn2,jn0, p1s,p2s,q2s,q0s,  p1e,p2e,q2e,q0e) ){
      aContactElem.insert( CContactElement(false,  in1,in2,jn2,jn0) );
    }
    if( IsContact_EE_CCD<BBOX>(          in2,in0,jn0,jn1, p2s,p0s,q0s,q1s,  p2e,p0e,q0e,q1e) ){
      aContactElem.insert( CContactElement(false,  in2,in0,jn0,jn1) );
    }
    if( IsContact_EE_CCD<BBOX>(          in2,in0,jn1,jn2, p2s,p0s,q1s,q2s,  p2e,p0e,q1e,q2e) ){
      aContactElem.insert( CContactElement(false,  in2,in0,jn1,jn2) );
    }
    if( IsContact_EE_CCD<BBOX>(          in2,in0,jn2,jn0, p2s,p0s,q2s,q0s,  p2e,p0e,q2e,q0e) ){
      aContactElem.insert( CContactElement(false,  in2,in0,jn2,jn0) );
    }
  }
}

template <typename BBOX>
void delfem2::GetContactElement_CCD(
	std::set<CContactElement>& aContactElem,
	//
	double dt,
	double delta,
	const std::vector<double>& aXYZ,
	const std::vector<double>& aUVW,
	const std::vector<unsigned int>& aTri,
	int ibvh,
	const std::vector<delfem2::CNodeBVH2>& aBVH,
	const std::vector<BBOX>& aBB)
{
  const int ichild0 = aBVH[ibvh].ichild[0];
  const int ichild1 = aBVH[ibvh].ichild[1];
  const bool is_leaf = (ichild1 == -1);
  if( is_leaf ) return;
  GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0,ichild1,aBVH,aBB);
  GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild0,        aBVH,aBB);
  GetContactElement_CCD(aContactElem, dt,delta, aXYZ,aUVW,aTri, ichild1,        aBVH,aBB);
}

// ---------------------------------------------------------------------------

namespace delfem2 {

template <typename REAL>
class CIntersectTriPair
{
public:
  int itri, jtri;
  CVec3<REAL> P[2];
};
  
}

// detect contact element with Continous Collision Detection (CCD)
template <typename BBOX>
void delfem2::GetIntersectTriPairs
(std::vector<delfem2::CIntersectTriPair<double>>& aIntersectTriPair,
 ////
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 int ibvh0, int ibvh1,
 const std::vector<delfem2::CNodeBVH2>& aBVH,
 const std::vector<BBOX>& aBB)
{
  assert( ibvh0 < (int)aBB.size() );
  assert( ibvh1 < (int)aBB.size() );
  if( !aBB[ibvh0].IsIntersect(aBB[ibvh1]) ) return;
  const int ichild0_0 = aBVH[ibvh0].ichild[0];
  const int ichild0_1 = aBVH[ibvh0].ichild[1];
  const int ichild1_0 = aBVH[ibvh1].ichild[0];
  const int ichild1_1 = aBVH[ibvh1].ichild[1];
  const bool is_leaf0 = (ichild0_1 == -1);
  const bool is_leaf1 = (ichild1_1 == -1);
  if(      !is_leaf0 && !is_leaf1 ){
    GetIntersectTriPairs(aIntersectTriPair, aXYZ,aTri, ichild0_0,ichild1_0, aBVH,aBB);
    GetIntersectTriPairs(aIntersectTriPair, aXYZ,aTri, ichild0_1,ichild1_0, aBVH,aBB);
    GetIntersectTriPairs(aIntersectTriPair, aXYZ,aTri, ichild0_0,ichild1_1, aBVH,aBB);
    GetIntersectTriPairs(aIntersectTriPair, aXYZ,aTri, ichild0_1,ichild1_1, aBVH,aBB);
  }
  else if( !is_leaf0 &&  is_leaf1 ){
    GetIntersectTriPairs(aIntersectTriPair, aXYZ,aTri, ichild0_0,ibvh1,     aBVH,aBB);
    GetIntersectTriPairs(aIntersectTriPair, aXYZ,aTri, ichild0_1,ibvh1,     aBVH,aBB);
  }
  else if(  is_leaf0 && !is_leaf1 ){
    GetIntersectTriPairs(aIntersectTriPair, aXYZ,aTri, ibvh0,    ichild1_0, aBVH,aBB);
    GetIntersectTriPairs(aIntersectTriPair, aXYZ,aTri, ibvh0,    ichild1_1, aBVH,aBB);
  }
  else if(  is_leaf0 &&  is_leaf1 ){
    const int itri = ichild0_0;
    const int jtri = ichild1_0;
    CVec3d P0,P1;
    bool res = isIntersectTriPair(P0,P1,
                                  itri, jtri, aTri, aXYZ);
    if( !res ){ return; }
    delfem2::CIntersectTriPair<double> itp;
    itp.itri = itri;
    itp.jtri = jtri;
    itp.P[0] = P0;
    itp.P[1] = P1;
    aIntersectTriPair.push_back(itp);
  }
}

template <typename BBOX>
void delfem2::GetIntersectTriPairs
(std::vector<delfem2::CIntersectTriPair<double>>& aIntersectTriPair,
 // ------------------------------------
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 int ibvh,
 const std::vector<delfem2::CNodeBVH2>& aBVH,
 const std::vector<BBOX>& aBB)
{
  const int ichild0 = aBVH[ibvh].ichild[0];
  const int ichild1 = aBVH[ibvh].ichild[1];
  const bool is_leaf = (ichild1 == -1);
  if( is_leaf ) return;
  GetIntersectTriPairs(aIntersectTriPair, aXYZ,aTri, ichild0,ichild1,aBVH,aBB);
  GetIntersectTriPairs(aIntersectTriPair, aXYZ,aTri, ichild0,        aBVH,aBB);
  GetIntersectTriPairs(aIntersectTriPair, aXYZ,aTri, ichild1,        aBVH,aBB);
}




#endif /* collisionTri_BvhVec3_hpp */
