/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/geo3_v23m34q.h"
#include "delfem2/defarapenergy_geo3.h"

// =======================================

DFM2_INLINE double delfem2::W_ArapEnergy(
	const std::vector<double>& aXYZ0,
	const std::vector<double>& aXYZ1,
	const std::vector<double>& aQuat1,
	const std::vector<unsigned int>& psup_ind,
	const std::vector<unsigned int>& psup)
{
  const size_t np = aXYZ0.size()/3;
  assert( aXYZ1.size() == np*3 );
  assert( aQuat1.size() == np*4 );
  assert( psup_ind.size() == np+1 );
  double w = 0.0;
  for(unsigned int ip=0;ip<np;++ip){
    CVec3d Pi(aXYZ0.data()+ip*3);
    CVec3d pi(aXYZ1.data()+ip*3);
    CQuatd Qi(aQuat1.data()+4*ip);
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const unsigned int jp = psup[ipsup];
      const CVec3d v0 = Qi*(CVec3d(aXYZ0.data()+jp*3)-Pi);
      CVec3d pj(aXYZ1.data()+jp*3);
      const CVec3d v1 = pj-pi;
      CVec3d v = v0-v1;
      w += v.dot(v);
      //      w += v1*v1;
    }
  }
  return 0.5*w;
}

DFM2_INLINE void delfem2::dW_ArapEnergy
(std::vector<double>& aRes,
 const std::vector<double>& aXYZ0,
 const std::vector<double>& aXYZ1,
 const std::vector<double>& aQuat1,
 const std::vector<unsigned int>& psup_ind,
 const std::vector<unsigned int>& psup)
{
  const size_t np = aXYZ0.size()/3;
  assert( aXYZ1.size() == np*3 );
  aRes.assign(np*3, 0.0);
  for(unsigned int ip=0;ip<np;++ip){
    const CVec3d Pi(aXYZ0.data()+ip*3);
    const CVec3d pi(aXYZ1.data()+ip*3);
    const CQuatd Qi(aQuat1.data()+ip*4);
    CMat3d LM; LM.setZero();
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const unsigned int jp = psup[ipsup];
      const CVec3d v0 = Qi*(CVec3d(aXYZ0.data()+jp*3)-Pi);
      CVec3d pj(aXYZ1.data()+jp*3);
      const CVec3d v1 = pj-pi;
      const CVec3d r = -(v1-v0);
      //      const CVec3d r = -v1;
      r.AddToScale(aRes.data()+ip*3, +1.);
      r.AddToScale(aRes.data()+jp*3, -1.);
    }
  }
}

DFM2_INLINE void delfem2::ddW_ArapEnergy(
	std::vector<double>& eM,
	const std::vector<unsigned int>& aIP,
	const std::vector<double>& aXYZ0,
	const std::vector<double>& aQuat1)
{
  const size_t nIP = aIP.size();
  const size_t nNg = nIP-1; // number of neighbor
  unsigned int ip = aIP[nNg];
  const CVec3d Pi(aXYZ0.data()+ip*3);
  CMat3d LM; LM.setZero();
  for(unsigned int iip=0;iip<nNg;++iip){
    const unsigned int jp = aIP[iip];
    const CVec3d v0 = (CVec3d(aXYZ0.data()+jp*3)-Pi);
    LM += Mat3_CrossCross(v0);
  }
  CMat3d LMi = LM.Inverse();
  CMat3d Mrot = CMat3d::Quat(aQuat1.data()+ip*4);
  //    LMi = R*LMi*R.Trans();
  eM.assign(nIP*nIP*9, 0.0);
  for(unsigned int jjp=0;jjp<nNg;++jjp){
    for(unsigned int kkp=0;kkp<nNg;++kkp){
      const CVec3d vj = (CVec3d(aXYZ0.data()+aIP[jjp]*3)-Pi);
      const CVec3d vk = (CVec3d(aXYZ0.data()+aIP[kkp]*3)-Pi);
      CMat3d L1 = Mrot*CMat3d::Spin(vk.p)*LMi*CMat3d::Spin(vj.p)*Mrot.transpose();
      //        L1 = CMat3d::Spin(vk.p)*LMi*CMat3d::Spin(vj.p);
      L1.AddToScale(eM.data()+(kkp*nIP+jjp)*9, -1.0);
      L1.AddToScale(eM.data()+(nNg*nIP+nNg)*9, -1.0);
      L1.AddToScale(eM.data()+(nNg*nIP+jjp)*9, +1.0);
      L1.AddToScale(eM.data()+(kkp*nIP+nNg)*9, +1.0);
    }
    {
      CMat3d L1 = CMat3d::Identity();
      L1.AddToScale(eM.data()+(jjp*nIP+jjp)*9, +1.0);
      L1.AddToScale(eM.data()+(nNg*nIP+nNg)*9, +1.0);
      L1.AddToScale(eM.data()+(nNg*nIP+jjp)*9, -1.0);
      L1.AddToScale(eM.data()+(jjp*nIP+nNg)*9, -1.0);
    }
  }
}
