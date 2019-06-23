//
//  self_collision_cloth.cpp
//  cloth_simple_contact
//
//  Created by Nobuyuki Umetani on 10/28/12.
//  Copyright (c) 2012 Nobuyuki Umetani. All rights reserved.
//


#include <stack>
#include <map>

#include "delfem2/cloth_selfcollision.h"
#include "delfem2/vec3.h"
#include "delfem2/bv.h"
#include "delfem2/bvh.h"
#include "delfem2/collisiontri_v3bvh.h"

// 撃力を計算
void SelfCollisionImpulse_Proximity
(std::vector<double>& aUVWm, // (in,out)velocity
 ////
 double delta,
 double stiffness,
 double dt,
 double mass,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 const std::vector<CContactElement>& aContactElem)
{
  for(int ice=0;ice<aContactElem.size();ice++){
    const CContactElement& ce = aContactElem[ice];
    const int ino0 = ce.ino0;
    const int ino1 = ce.ino1;
    const int ino2 = ce.ino2;
    const int ino3 = ce.ino3;
    CVector3 p0( aXYZ[ ino0*3+0], aXYZ[ ino0*3+1], aXYZ[ ino0*3+2] );
    CVector3 p1( aXYZ[ ino1*3+0], aXYZ[ ino1*3+1], aXYZ[ ino1*3+2] );
    CVector3 p2( aXYZ[ ino2*3+0], aXYZ[ ino2*3+1], aXYZ[ ino2*3+2] );
    CVector3 p3( aXYZ[ ino3*3+0], aXYZ[ ino3*3+1], aXYZ[ ino3*3+2] );
    CVector3 v0( aUVWm[ino0*3+0], aUVWm[ino0*3+1], aUVWm[ino0*3+2] );
    CVector3 v1( aUVWm[ino1*3+0], aUVWm[ino1*3+1], aUVWm[ino1*3+2] );
    CVector3 v2( aUVWm[ino2*3+0], aUVWm[ino2*3+1], aUVWm[ino2*3+2] );
    CVector3 v3( aUVWm[ino3*3+0], aUVWm[ino3*3+1], aUVWm[ino3*3+2] );
    if( ce.is_fv ){ // face-vtx      
      double w0,w1;
      {
        double dist = DistanceFaceVertex(p0, p1, p2, p3, w0,w1);
        if( w0 < 0 || w0 > 1 ) continue;
        if( w1 < 0 || w1 > 1 ) continue;
        if( dist > delta ) continue;
      }
      double w2 = 1.0 - w0 - w1;
      CVector3 pc = w0*p0 + w1*p1 + w2*p2;
      CVector3 norm = p3-pc; norm.SetNormalizedVector();
      double p_depth = delta - Dot(p3-pc,norm); // penetration depth 
      double rel_v = Dot(v3-w0*v0-w1*v1-w2*v2,norm);
      if( rel_v > 0.1*p_depth/dt ) continue;
      double imp_el = dt*stiffness*p_depth;
      double imp_ie = mass*(0.1*p_depth/dt-rel_v);
      double imp_min = ( imp_el < imp_ie ) ? imp_el : imp_ie;
      double imp_mod = 2*imp_min / (1+w0*w0+w1*w1+w2*w2);
      imp_mod /= mass;
      imp_mod *= 0.25;
      aUVWm[ino0*3+0] += -norm.x*imp_mod*w0;
      aUVWm[ino0*3+1] += -norm.y*imp_mod*w0;
      aUVWm[ino0*3+2] += -norm.z*imp_mod*w0;
      aUVWm[ino1*3+0] += -norm.x*imp_mod*w1;
      aUVWm[ino1*3+1] += -norm.y*imp_mod*w1;
      aUVWm[ino1*3+2] += -norm.z*imp_mod*w1;
      aUVWm[ino2*3+0] += -norm.x*imp_mod*w2;
      aUVWm[ino2*3+1] += -norm.y*imp_mod*w2;
      aUVWm[ino2*3+2] += -norm.z*imp_mod*w2;
      aUVWm[ino3*3+0] += +norm.x*imp_mod;
      aUVWm[ino3*3+1] += +norm.y*imp_mod;
      aUVWm[ino3*3+2] += +norm.z*imp_mod;
    }
    else{ // edge-edge
      double w01,w23;
      {
        double dist = DistanceEdgeEdge(p0, p1, p2, p3, w01,w23);
        if( w01 < 0 || w01 > 1 ) continue;
        if( w23 < 0 || w23 > 1 ) continue;
        if( dist > delta ) continue;
      }
      CVector3 c01 = (1-w01)*p0 + w01*p1;
      CVector3 c23 = (1-w23)*p2 + w23*p3;
      CVector3 norm = (c23-c01); norm.SetNormalizedVector();
      double p_depth = delta - (c23-c01).Length();
      double rel_v = Dot((1-w23)*v2+w23*v3-(1-w01)*v0-w01*v1,norm);
      if( rel_v > 0.1*p_depth/dt ) continue;
      double imp_el = dt*stiffness*p_depth;
      double imp_ie = mass*(0.1*p_depth/dt-rel_v);
      double imp_min = ( imp_el < imp_ie ) ? imp_el : imp_ie;
      double imp_mod = 2*imp_min / ( w01*w01+(1-w01)*(1-w01) + w23*w23+(1-w23)*(1-w23) );
      imp_mod /= mass;
      imp_mod *= 0.25;      
      aUVWm[ino0*3+0] += -norm.x*imp_mod*(1-w01);
      aUVWm[ino0*3+1] += -norm.y*imp_mod*(1-w01);
      aUVWm[ino0*3+2] += -norm.z*imp_mod*(1-w01);
      aUVWm[ino1*3+0] += -norm.x*imp_mod*w01;
      aUVWm[ino1*3+1] += -norm.y*imp_mod*w01;
      aUVWm[ino1*3+2] += -norm.z*imp_mod*w01;
      aUVWm[ino2*3+0] += +norm.x*imp_mod*(1-w23);
      aUVWm[ino2*3+1] += +norm.y*imp_mod*(1-w23);
      aUVWm[ino2*3+2] += +norm.z*imp_mod*(1-w23);
      aUVWm[ino3*3+0] += +norm.x*imp_mod*w23;
      aUVWm[ino3*3+1] += +norm.y*imp_mod*w23;
      aUVWm[ino3*3+2] += +norm.z*imp_mod*w23;
    }
  }
}


// Impulseの計算
void SelfCollisionImpulse_CCD
(std::vector<double>& aUVWm, // (in,out)velocity
 ////
 double delta,
 double stiffness,
 double dt,
 double mass,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
 const std::vector<CContactElement>& aContactElem)
{
  for(int ice=0;ice<aContactElem.size();ice++){
    const CContactElement& ce = aContactElem[ice];
    const int ino0 = ce.ino0;
    const int ino1 = ce.ino1;
    const int ino2 = ce.ino2;
    const int ino3 = ce.ino3;
    CVector3 p0( aXYZ[ ino0*3+0], aXYZ[ ino0*3+1], aXYZ[ ino0*3+2] );
    CVector3 p1( aXYZ[ ino1*3+0], aXYZ[ ino1*3+1], aXYZ[ ino1*3+2] );
    CVector3 p2( aXYZ[ ino2*3+0], aXYZ[ ino2*3+1], aXYZ[ ino2*3+2] );
    CVector3 p3( aXYZ[ ino3*3+0], aXYZ[ ino3*3+1], aXYZ[ ino3*3+2] );
    CVector3 v0( aUVWm[ino0*3+0], aUVWm[ino0*3+1], aUVWm[ino0*3+2] );
    CVector3 v1( aUVWm[ino1*3+0], aUVWm[ino1*3+1], aUVWm[ino1*3+2] );
    CVector3 v2( aUVWm[ino2*3+0], aUVWm[ino2*3+1], aUVWm[ino2*3+2] );
    CVector3 v3( aUVWm[ino3*3+0], aUVWm[ino3*3+1], aUVWm[ino3*3+2] );    
    double t = FindCoplanerInterp(p0,p1,p2,p3, p0+v0,p1+v1,p2+v2,p3+v3);
    if( t < 0 || t > 1 ) continue;
    if( ce.is_fv ){ // face-vtx
      double w0,w1;
      {        
        CVector3 p0m = p0 + t*v0;
        CVector3 p1m = p1 + t*v1;
        CVector3 p2m = p2 + t*v2;
        CVector3 p3m = p3 + t*v3;
        double dist = DistanceFaceVertex(p0m, p1m, p2m, p3m, w0,w1);
        if( w0 < 0 || w0 > 1 ) continue;
        if( w1 < 0 || w1 > 1 ) continue;
        if( dist > delta ) continue;
      }
      double w2 = 1.0 - w0 - w1;
      CVector3 pc = w0*p0 + w1*p1 + w2*p2;
      CVector3 norm = p3 - pc; norm.SetNormalizedVector();
      double rel_v = Dot(v3-w0*v0-w1*v1-w2*v2,norm); // relative velocity (positive if separating)
      if( rel_v > 0.1*delta/dt ) continue; // separating
      double imp = mass*(0.1*delta/dt-rel_v);
      double imp_mod = 2*imp/(1.0+w0*w0+w1*w1+w2*w2);
      imp_mod /= mass;
      imp_mod *= 0.1;
      aUVWm[ino0*3+0] += -norm.x*imp_mod*w0;
      aUVWm[ino0*3+1] += -norm.y*imp_mod*w0;
      aUVWm[ino0*3+2] += -norm.z*imp_mod*w0;
      aUVWm[ino1*3+0] += -norm.x*imp_mod*w1;
      aUVWm[ino1*3+1] += -norm.y*imp_mod*w1;
      aUVWm[ino1*3+2] += -norm.z*imp_mod*w1;
      aUVWm[ino2*3+0] += -norm.x*imp_mod*w2;
      aUVWm[ino2*3+1] += -norm.y*imp_mod*w2;
      aUVWm[ino2*3+2] += -norm.z*imp_mod*w2;
      aUVWm[ino3*3+0] += +norm.x*imp_mod;
      aUVWm[ino3*3+1] += +norm.y*imp_mod;
      aUVWm[ino3*3+2] += +norm.z*imp_mod;
    }
    else{ // edge-edge
      double w01,w23;
      {
        CVector3 p0m = p0 + t*v0;
        CVector3 p1m = p1 + t*v1;
        CVector3 p2m = p2 + t*v2;
        CVector3 p3m = p3 + t*v3;
        double dist = DistanceEdgeEdge(p0m, p1m, p2m, p3m, w01,w23);
        if( w01 < 0 || w01 > 1 ) continue;
        if( w23 < 0 || w23 > 1 ) continue;
        if( dist > delta ) continue;
      }      
      CVector3 c01 = (1-w01)*p0 + w01*p1;
      CVector3 c23 = (1-w23)*p2 + w23*p3;
      CVector3 norm = (c23-c01); norm.SetNormalizedVector();
      double rel_v = Dot((1-w23)*v2+w23*v3-(1-w01)*v0-w01*v1,norm);
      if( rel_v > 0.1*delta/dt ) continue; // separating
      double imp = mass*(0.1*delta/dt-rel_v); // reasonable
      double imp_mod = 2*imp/( w01*w01+(1-w01)*(1-w01) + w23*w23+(1-w23)*(1-w23) );
      imp_mod /= mass;
      imp_mod *= 0.1;
      aUVWm[ino0*3+0] += -norm.x*imp_mod*(1-w01);
      aUVWm[ino0*3+1] += -norm.y*imp_mod*(1-w01);
      aUVWm[ino0*3+2] += -norm.z*imp_mod*(1-w01);
      aUVWm[ino1*3+0] += -norm.x*imp_mod*w01;
      aUVWm[ino1*3+1] += -norm.y*imp_mod*w01;
      aUVWm[ino1*3+2] += -norm.z*imp_mod*w01;
      aUVWm[ino2*3+0] += +norm.x*imp_mod*(1-w23);
      aUVWm[ino2*3+1] += +norm.y*imp_mod*(1-w23);
      aUVWm[ino2*3+2] += +norm.z*imp_mod*(1-w23);
      aUVWm[ino3*3+0] += +norm.x*imp_mod*w23;
      aUVWm[ino3*3+1] += +norm.y*imp_mod*w23;
      aUVWm[ino3*3+2] += +norm.z*imp_mod*w23;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// RIZを更新する
void MakeRigidImpactZone
(std::vector< std::set<int> >& aRIZ, // (in,ou)RIZに属する節点のインデックスの集合の配列
 const std::vector<CContactElement>& aContactElem, // 自己交差する接触要素の配列
// const CJaggedArray& aEdge
 const std::vector<int>& psup_ind,
 const std::vector<int>& psup) // 三角形メッシュの辺の配列    
{
  for(int ice=0;ice<aContactElem.size();ice++){
    const CContactElement& ce = aContactElem[ice];
    const int n[4] = {ce.ino0, ce.ino1, ce.ino2, ce.ino3};
    std::set<int> ind_inc; // 接触要素が接するRIZの集合
    for(int i=0;i<4;i++){
      const int ino = n[i];
      for(int iriz=0;iriz<aRIZ.size();iriz++){
        if( aRIZ[iriz].find(ino) != aRIZ[iriz].end() ){
          ind_inc.insert(iriz);
        }
        else{
          for(int iedge=psup_ind[ino];iedge<psup_ind[ino+1];iedge++){
            int jno = psup[iedge];
            if( aRIZ[iriz].find(jno) != aRIZ[iriz].end() ){
              ind_inc.insert(iriz);  break;
            }
          }
        }
      }
    }
    if( ind_inc.size() == 0 ){ // 接触要素はどのRIZにも属していない
      int ind0 = (int)aRIZ.size();
      aRIZ.resize(ind0+1);
      aRIZ[ind0].insert(n[0]); aRIZ[ind0].insert(n[1]); aRIZ[ind0].insert(n[2]); aRIZ[ind0].insert(n[3]);
    }
    else if( ind_inc.size() == 1 ){ // 接触要素は一つのRIZに接する
      int ind0 = *(ind_inc.begin());
      aRIZ[ind0].insert(n[0]); aRIZ[ind0].insert(n[1]); aRIZ[ind0].insert(n[2]); aRIZ[ind0].insert(n[3]);
    }
    else{ // overlapping two reagion，接触要素が複数のRIZに接するー＞接する複数のRIZをマージする
      std::vector< std::set<int> > aRIZ1; // マージされた後のRIZの配列
      for(int iriz=0;iriz<aRIZ.size();iriz++){ // 接さないRIZをコピー
        if( ind_inc.find(iriz) != ind_inc.end() ) continue;
        aRIZ1.push_back( aRIZ[iriz] );
      }
      // マージしたRIZを，aRIZ1の最後に追加
      int ind0 = (int)aRIZ1.size();
      aRIZ1.resize(ind0+1);
      for(std::set<int>::iterator itr=ind_inc.begin();itr!=ind_inc.end();itr++){
        int ind1 = *itr;
        assert( ind1 < aRIZ.size() );
        for(std::set<int>::iterator jtr=aRIZ[ind1].begin();jtr!=aRIZ[ind1].end();jtr++){
          aRIZ1[ind0].insert(*jtr);
        }
      }
      aRIZ1[ind0].insert(n[0]); aRIZ1[ind0].insert(n[1]); aRIZ1[ind0].insert(n[2]); aRIZ1[ind0].insert(n[3]);
      aRIZ = aRIZ1;
    }
  }
}


// t is a tmporary buffer size of 9
static inline void CalcInvMat3(double ainv[], const double a[])
{
	const double det =
  + a[0]*a[4]*a[8] + a[3]*a[7]*a[2] + a[6]*a[1]*a[5]
  - a[0]*a[7]*a[5] - a[6]*a[4]*a[2] - a[3]*a[1]*a[8];
	const double inv_det = 1.0/det;
  
	ainv[0] = inv_det*(a[4]*a[8]-a[5]*a[7]);
	ainv[1] = inv_det*(a[2]*a[7]-a[1]*a[8]);
	ainv[2] = inv_det*(a[1]*a[5]-a[2]*a[4]);
  
	ainv[3] = inv_det*(a[5]*a[6]-a[3]*a[8]);
	ainv[4] = inv_det*(a[0]*a[8]-a[2]*a[6]);
	ainv[5] = inv_det*(a[2]*a[3]-a[0]*a[5]);
  
	ainv[6] = inv_det*(a[3]*a[7]-a[4]*a[6]);
	ainv[7] = inv_det*(a[1]*a[6]-a[0]*a[7]);
	ainv[8] = inv_det*(a[0]*a[4]-a[1]*a[3]);
}

void ApplyRigidImpactZone
(std::vector<double>& aUVWm, // (in,out)RIZで更新された中間速度
 ////
 const std::vector< std::set<int> >& aRIZ,  // (in)各RIZに属する節点の集合(set)の配列
 const std::vector<double>& aXYZ, // (in) 前ステップの節点の位置の配列
 const std::vector<double>& aUVWm0) // (in) RIZを使う前の中間速度
{
  for(int iriz=0;iriz<aRIZ.size();iriz++){
    std::vector<int> aInd; // index of points belong to this RIZ
    for(std::set<int>::iterator jtr=aRIZ[iriz].begin();jtr!=aRIZ[iriz].end();jtr++){
      aInd.push_back(*jtr);
    }
    CVector3 gc(0,0,0); // 重心位置
    CVector3 av(0,0,0); // 平均速度
    for(int iv=0;iv<aInd.size();iv++){
      int ino = aInd[iv];
      gc += CVector3(aXYZ[  ino*3+0],aXYZ[  ino*3+1],aXYZ[  ino*3+2]);
      av += CVector3(aUVWm0[ino*3+0],aUVWm0[ino*3+1],aUVWm0[ino*3+2]);
    }
    gc /= (double)aInd.size();
    av /= (double)aInd.size();
    CVector3 L(0,0,0); // 角運動量
    double I[9] = {0,0,0, 0,0,0, 0,0,0}; // 慣性テンソル
    for(int iv=0;iv<aInd.size();iv++){
      int ino = aInd[iv];
      CVector3 p(aXYZ[  ino*3+0],aXYZ[  ino*3+1],aXYZ[  ino*3+2]);
      CVector3 v(aUVWm0[ino*3+0],aUVWm0[ino*3+1],aUVWm0[ino*3+2]);
      L += Cross(p-gc,v-av);
      CVector3 q = p-gc;
      I[0] += v*v - q[0]*q[0];  I[1] +=     - q[0]*q[1];  I[2] +=     - q[0]*q[2];
      I[3] +=     - q[1]*q[0];  I[4] += v*v - q[1]*q[1];  I[5] +=     - q[1]*q[2];
      I[6] +=     - q[2]*q[0];  I[7] +=     - q[2]*q[1];  I[8] += v*v - q[2]*q[2];
    }
    // 角速度を求める
    double Iinv[9];
    CalcInvMat3(Iinv,I);
    CVector3 omg;
    omg.x = Iinv[0]*L.x + Iinv[1]*L.y + Iinv[2]*L.z;
    omg.y = Iinv[3]*L.x + Iinv[4]*L.y + Iinv[5]*L.z;
    omg.z = Iinv[6]*L.x + Iinv[7]*L.y + Iinv[8]*L.z;
    // 中間速度の更新
    for(int iv=0;iv<aInd.size();iv++){
      int ino = aInd[iv];
      CVector3 p(aXYZ[  ino*3+0],aXYZ[  ino*3+1],aXYZ[  ino*3+2]);
      CVector3 rot = -Cross(p-gc,omg);
      aUVWm[ino*3+0] = av.x + rot.x;
      aUVWm[ino*3+1] = av.y + rot.y;
      aUVWm[ino*3+2] = av.z + rot.z;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 衝突が解消された中間速度を返す
void GetIntermidiateVelocityContactResolved
(std::vector<double>& aUVWm,
 bool& is_impulse_applied,
 ////
 double dt,
 double contact_clearance,
 double mass_point,
 double cloth_contact_stiffness,
 const std::vector<double>& aXYZ,
 const std::vector<unsigned int>& aTri,
// const CJaggedArray& aEdge,
 const std::vector<int>& psup_ind,
 const std::vector<int>& psup,
 int iroot_bvh,
 const std::vector<CNodeBVH>& aNodeBVH,
 std::vector<CBV3D_AABB>& aBB)
{
  {
    std::vector<CContactElement> aContactElem;
    {
      BuildBoundingBoxesBVH(iroot_bvh,
                            contact_clearance,
                            aXYZ,aTri,3,aNodeBVH,aBB);
      std::set<CContactElement> setCE;
      GetContactElement_Proximity(setCE,
                                  contact_clearance,
                                  aXYZ,aTri,
                                  iroot_bvh,
                                  aNodeBVH,aBB); // output
      aContactElem.assign(setCE.begin(),setCE.end());
//      aContactElem.clear();
//      for(std::set<CContactElement>::iterator itr=setCE.begin();itr!=setCE.end();itr++){
//        aContactElem.push_back(*itr);
//      }
      std::cout << "  Proximity      Contact Elem Size: " << aContactElem.size() << std::endl;
    }
    is_impulse_applied = aContactElem.size() > 0;
    SelfCollisionImpulse_Proximity(aUVWm,
                              contact_clearance,
                              cloth_contact_stiffness,
                              dt,
                              mass_point,
                              aXYZ,aTri,
                              aContactElem);
  }
  ////////
  for(int itr=0;itr<5;itr++){
    std::vector<CContactElement> aContactElem;
    {
      BuildBoundingBoxesBVH_Dynamic(iroot_bvh,
                                    dt,
                                    aXYZ,aUVWm,aTri,aNodeBVH,aBB);
      std::set<CContactElement> setCE;
      GetContactElement_CCD(setCE,
                            dt,contact_clearance,
                            aXYZ,aUVWm,aTri,
                            iroot_bvh,
                            aNodeBVH,aBB); // output
      aContactElem.assign(setCE.begin(),setCE.end());
//      aContactElem.clear();
//      for(std::set<CContactElement>::iterator itr=setCE.begin();itr!=setCE.end();itr++){
//        aContactElem.push_back(*itr);
//      }
    }
      std::cout << "  CCD iter: " << itr << "    Contact Elem Size: " << aContactElem.size() << std::endl;    
    if( aContactElem.size() == 0 ){ return; }
    is_impulse_applied = is_impulse_applied || (aContactElem.size() > 0);    
    SelfCollisionImpulse_CCD(aUVWm,
                              contact_clearance,
                              cloth_contact_stiffness,
                              dt,
                              mass_point,
                              aXYZ,aTri,
                              aContactElem);
  }
  std::vector<double> aUVWm0 = aUVWm;
  std::vector< std::set<int> > aRIZ;
  for(int itr=0;itr<100;itr++){
    std::vector<CContactElement> aContactElem;    
    {
      BuildBoundingBoxesBVH_Dynamic(iroot_bvh,
                                    dt,
                                    aXYZ,aUVWm,aTri,aNodeBVH,aBB);
      std::set<CContactElement> setCE;
      GetContactElement_CCD(setCE,
                            dt,contact_clearance,
                            aXYZ,aUVWm,aTri,
                            iroot_bvh,
                            aNodeBVH,aBB); // output
      aContactElem.assign(setCE.begin(),setCE.end());
//      for(std::set<CContactElement>::iterator itr=setCE.begin();itr!=setCE.end();itr++){
//        aContactElem.push_back(*itr);
//      }
    }
    int nnode_riz = 0;
    for(int iriz=0;iriz<aRIZ.size();iriz++){
      nnode_riz += aRIZ[iriz].size();
    }
    std::cout << "  RIZ iter: " << itr << "    Contact Elem Size: " << aContactElem.size() << "   NNode In RIZ: " << nnode_riz << std::endl;
    if( aContactElem.size() == 0 ){
      std::cout << "Resolved All Collisions : " << std::endl;
      break;
    }
    MakeRigidImpactZone(aRIZ, aContactElem, psup_ind,psup);
    ApplyRigidImpactZone(aUVWm, aRIZ,aXYZ,aUVWm0);
  }
}

