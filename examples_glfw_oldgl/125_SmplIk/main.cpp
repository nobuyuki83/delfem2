/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief SMPL model
 * @details skinning
 */

#include <cstdlib>
#include <random>
#include <time.h>
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mat4.h"
#include "delfem2/quat.h"
#include "delfem2/mats.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vecxitrsol.h"
//
#include "delfem2/rig_v3q.h"
#include "delfem2/v23m3q.h"
//
#include "delfem2/cnpy/smpl_cnpy.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_v23.h"
#include "delfem2/opengl/glold_rig_v23q.h"

namespace dfm2 = delfem2;

// -------------------

void UpdateRotationsByClusterMatching
(std::vector<double>& aQuat1,
 const std::vector<double>& aXYZ0,
 const std::vector<double>& aXYZ1,
 const std::vector<unsigned int>& psup_ind,
 const std::vector<unsigned int>& psup)
{
  const unsigned int np = aXYZ0.size()/3;
  for(unsigned int ip=0;ip<np;++ip){
    const dfm2::CVec3d Pi(aXYZ0.data()+ip*3);
    const dfm2::CVec3d pi(aXYZ1.data()+ip*3);
    const dfm2::CQuatd Qi(aQuat1.data()+ip*4);
    dfm2::CMat3d Mat; Mat.SetZero();
    dfm2::CVec3d rhs; rhs.SetZero();
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const unsigned int jp = psup[ipsup];
      const dfm2::CVec3d v0 = Qi*(dfm2::CVec3d(aXYZ0.data()+jp*3)-Pi);
      const dfm2::CVec3d v1 = dfm2::CVec3d(aXYZ1.data()+jp*3)-pi-v0;
      Mat += dfm2::Mat3_CrossCross(v0);
      rhs += v1^v0;
    }
    dfm2::CVec3d sol = Mat.Inverse()*rhs;
    dfm2::CQuatd q0 = dfm2::Quat_CartesianAngle(sol);
    dfm2::CQuatd q1 = q0*dfm2::CQuatd(aQuat1.data()+ip*4);
    q1.CopyTo(aQuat1.data()+ip*4);
  }
}

void SensitivityRigSkin
 (std::vector<double>& dXYZ12,
  unsigned int ib_s,
  unsigned int idim_s,
  const std::vector<dfm2::CRigBone> aBone1,
  const std::vector<double>& aXYZ0,
  const std::vector<double>& aW)
{
  const unsigned int nBone = aBone1.size();
  std::vector<double> aM(nBone*16);
  {
    for(std::size_t ibone=0;ibone<aBone1.size();++ibone){
      dfm2::CMat4d m01 = dfm2::CMat4d::Scale(aBone1[ibone].scale);
      m01 = dfm2::CMat4d::Quat(aBone1[ibone].quatRelativeRot) * m01;
      if( ibone == ib_s ){
        dfm2::CMat3d dn1 = dfm2::CMat3d::Spin(dfm2::CVec3d::Axis(idim_s).p) + dfm2::CMat3d::Identity();
        dfm2::CMat4d dm1 = dfm2::CMat4d::Mat3(dn1.mat);
        m01 = dm1 * m01;
      }
      m01 = dfm2::CMat4d::Translate(aBone1[ibone].transRelative) * m01;
      const int ibone_p = aBone1[ibone].ibone_parent;
      if( ibone_p < 0 || ibone_p >= (int)aBone1.size() ){ // root bone
        dfm2::Copy_Mat4( aM.data()+ibone*16, m01.mat );
        continue;
      }
      dfm2::MatMat4(aM.data()+ibone*16,
                    aM.data()+ibone_p*16, m01.mat);
    }
    for(std::size_t ibone=0;ibone<aBone1.size();++ibone){
      for(int i=0;i<16;++i){
        aM[ibone*16+i] -= aBone1[ibone].affmat3Global[i];
      }
    }
  }

  const unsigned int np = aXYZ0.size()/3;
  dXYZ12.assign(np*3,0.0);
  for(int ip=0;ip<np;++ip){
    double p0a[4] = {aXYZ0[ip*3+0], aXYZ0[ip*3+1], aXYZ0[ip*3+2], 1.0};
    for(int ib=0;ib<nBone;++ib){
      double p0b[4]; dfm2::MatVec4(p0b,
                                   aBone1[ib].invBindMat, p0a);
      double p0c[4]; dfm2::MatVec4(p0c,
                                   aM.data()+ib*16, p0b);
      dXYZ12[ip*3+0] += aW[ip*nBone+ib]*p0c[0];
      dXYZ12[ip*3+1] += aW[ip*nBone+ib]*p0c[1];
      dXYZ12[ip*3+2] += aW[ip*nBone+ib]*p0c[2];
    }
  }
}


void Check_SensitivityRigSkin(
    std::vector<dfm2::CRigBone> aBone1,
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aW)
{
  const unsigned int nBone = aBone1.size();
  unsigned int ib_s = (unsigned int)(nBone*(rand()/(RAND_MAX+1.0)));
  unsigned int idim_s = (unsigned int)(3*(rand()/(RAND_MAX+1.0)));

  const double eps = 1.0e-4;
  std::vector<dfm2::CRigBone> aBone2 = aBone1;
  {
    dfm2::CQuatd dq = dfm2::Quat_CartesianAngle(eps*dfm2::CVec3d::Axis(idim_s));
    dfm2::CQuatd q0 = dq*dfm2::CQuatd(aBone2[ib_s].quatRelativeRot);
    q0.CopyTo(aBone2[ib_s].quatRelativeRot);
  }
  std::vector<double> dXYZ12;
  SensitivityRigSkin(dXYZ12,
                     ib_s, idim_s,
                     aBone1, aXYZ0, aW);
  {
    std::vector<double> aXYZ1;
    dfm2::UpdateBoneRotTrans(aBone1);
    dfm2::Skinning_LBS(aXYZ1,
                       aXYZ0, aBone1, aW);
    // ----------------
    std::vector<double> aXYZ2;
    dfm2::UpdateBoneRotTrans(aBone2);
    dfm2::Skinning_LBS(aXYZ2,
                       aXYZ0, aBone2, aW);
    // ----------------
    double max_ratio = 0.0;
    for(int i=0;i<aXYZ0.size();++i){
      double val0 = (aXYZ2[i] - aXYZ1[i])/eps;
      double val1 =  dXYZ12[i];
      double ratio = fabs(val0-val1)/(fabs(val1)+1.0);
      max_ratio = (ratio >max_ratio ) ? ratio : max_ratio;
    }
    std::cout << "max ratio: " << max_ratio << std::endl;
  }
}


void W_ArapEdgeDiff(
    std::vector<double>& aDiff,
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aXYZ1,
    const std::vector<double>& aQuat1,
    const std::vector<unsigned int>& psup_ind,
    const std::vector<unsigned int>& psup)
{
  const unsigned int np = aXYZ0.size()/3;
  assert( aXYZ1.size() == np*3 );
  assert( aQuat1.size() == np*4 );
  assert( psup_ind.size() == np+1 );
  aDiff.resize(psup.size()*3);
  for(unsigned int ip=0;ip<np;++ip){
    dfm2::CVec3d Pi(aXYZ0.data()+ip*3);
    dfm2::CVec3d pi(aXYZ1.data()+ip*3);
    dfm2::CQuatd Qi(aQuat1.data()+4*ip);
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const unsigned int jp = psup[ipsup];
      dfm2::CVec3d Pj(aXYZ0.data()+jp*3);
      dfm2::CVec3d pj(aXYZ1.data()+jp*3);
      dfm2::CVec3d v = Qi*(Pj-Pi)-(pj-pi);
      v.CopyTo(aDiff.data()+ipsup*3);
    }
  }
}



void dW_ArapEdgeDiff(
    std::vector<double>& dDiff,
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aXYZ1,
    const std::vector<double>& aQuat1,
    const std::vector<double>& dXYZ12,
    const std::vector<unsigned int>& psup_ind,
    const std::vector<unsigned int>& psup)
{
  const unsigned int np = aXYZ0.size()/3;
  assert( aXYZ1.size() == np*3 );
  assert( aQuat1.size() == np*4 );
  assert( psup_ind.size() == np+1 );
  dDiff.resize(psup.size()*3);
  for(unsigned int ip=0;ip<np;++ip){
    const dfm2::CVec3d Pi(aXYZ0.data()+ip*3);
    const dfm2::CQuatd Qi(aQuat1.data()+ip*4);
    dfm2::CMat3d Mat; Mat.SetZero();
    dfm2::CVec3d rhs; rhs.SetZero();
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const unsigned int jp = psup[ipsup];
      const dfm2::CVec3d v0 = Qi*(dfm2::CVec3d(aXYZ0.data()+jp*3)-Pi);
      const dfm2::CVec3d v1 = dfm2::CVec3d(dXYZ12.data()+jp*3) - dfm2::CVec3d(dXYZ12.data()+ip*3);
      Mat += dfm2::Mat3_CrossCross(v0);
      rhs += (v0^v1);
    }
    const dfm2::CVec3d sol = Mat.Inverse()*rhs;
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const unsigned int jp = psup[ipsup];
      const dfm2::CVec3d v0 = Qi*(dfm2::CVec3d(aXYZ0.data()+jp*3)-Pi);
      const dfm2::CVec3d v1 = dfm2::CVec3d(dXYZ12.data()+jp*3) - dfm2::CVec3d(dXYZ12.data()+ip*3);
      dfm2::CVec3d d = +(v0^sol)-v1;
      d.CopyTo(dDiff.data()+ipsup*3);
    }
  }
}

void Check_ArapEdgeDiff(const std::vector<double>& aXYZ0,
                        const std::vector<double>& aXYZ1,
                        const std::vector<double>& aQuat0,
                        const std::vector<double>& dXYZ12,
                        const std::vector<unsigned int>& psup_ind,
                        const std::vector<unsigned int>& psup)
{
  std::vector<double> aQuat1 = aQuat0;
  for(int itr=0;itr<20;++itr){
    UpdateRotationsByClusterMatching(aQuat1,
                                     aXYZ0,aXYZ1,psup_ind,psup);
  }

  std::vector<double> aDiff1;
  W_ArapEdgeDiff(aDiff1,
                 aXYZ0, aXYZ1, aQuat1, psup_ind, psup);

  const double eps = 1.0e-5;
  std::vector<double> aXYZ2 = aXYZ1;
  for(int i=0;i<aXYZ2.size();++i){ aXYZ2[i] += eps*dXYZ12[i]; }
  
  std::vector<double> aQuat2 = aQuat1;
  UpdateRotationsByClusterMatching(aQuat2,
                                   aXYZ0, aXYZ2, psup_ind, psup);
    
  std::vector<double> aDiff2;
  W_ArapEdgeDiff(aDiff2,
                 aXYZ0, aXYZ2, aQuat2, psup_ind, psup);
  
  std::vector<double> dDiff12;
  dW_ArapEdgeDiff(dDiff12,
                  aXYZ0, aXYZ1, aQuat1, dXYZ12, psup_ind, psup);
  
  assert( aDiff1.size() == dDiff12.size() );
  assert( aDiff2.size() == dDiff12.size() );
  
  double max_ratio = 0.0;
  for(int i=0;i<dDiff12.size();++i){
    double val0 = dDiff12[i];
    double val1 = (aDiff2[i]-aDiff1[i])/eps;
    double ratio = fabs(val1-val0)/(fabs(val0)+1.0);
    max_ratio = (ratio > max_ratio) ? ratio : max_ratio;
  }
  std::cout << "max_diff_ratio: " << max_ratio << std::endl;
}

double W_ArapEnergy(const std::vector<double>& aXYZ0,
                    const std::vector<double>& aXYZ1,
                    const std::vector<double>& aQuat1,
                    const std::vector<unsigned int>& psup_ind,
                    const std::vector<unsigned int>& psup)
{
  const unsigned int np = aXYZ0.size()/3;
  assert( aXYZ1.size() == np*3 );
  assert( aQuat1.size() == np*4 );
  assert( psup_ind.size() == np+1 );
  double w = 0.0;
  for(unsigned int ip=0;ip<np;++ip){
    dfm2::CVec3d Pi(aXYZ0.data()+ip*3);
    dfm2::CVec3d pi(aXYZ1.data()+ip*3);
    dfm2::CQuatd Qi(aQuat1.data()+4*ip);
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const unsigned int jp = psup[ipsup];
      const dfm2::CVec3d v0 = Qi*(dfm2::CVec3d(aXYZ0.data()+jp*3)-Pi);
      dfm2::CVec3d pj(aXYZ1.data()+jp*3);
      const dfm2::CVec3d v1 = pj-pi;
      dfm2::CVec3d v = v0-v1;
      w += v*v;
//      w += v1*v1;
    }
  }
  return 0.5*w;
}

void dW_ArapEnergy
 (std::vector<double>& aRes,
  const std::vector<double>& aXYZ0,
  const std::vector<double>& aXYZ1,
  const std::vector<double>& aQuat1,
  const std::vector<unsigned int>& psup_ind,
  const std::vector<unsigned int>& psup)
{
  const unsigned int np = aXYZ0.size()/3;
  assert( aXYZ1.size() == np*3 );
  aRes.assign(np*3, 0.0);
  for(unsigned int ip=0;ip<np;++ip){
    const dfm2::CVec3d Pi(aXYZ0.data()+ip*3);
    const dfm2::CVec3d pi(aXYZ1.data()+ip*3);
    const dfm2::CQuatd Qi(aQuat1.data()+ip*4);
    dfm2::CMat3d LM; LM.SetZero();
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const unsigned int jp = psup[ipsup];
      const dfm2::CVec3d v0 = Qi*(dfm2::CVec3d(aXYZ0.data()+jp*3)-Pi);
      dfm2::CVec3d pj(aXYZ1.data()+jp*3);
      const dfm2::CVec3d v1 = pj-pi;
//      const dfm2::CVec3d r = -v1;
      const dfm2::CVec3d r = -(v1-v0);
      aRes[ip*3+0] += r.x();
      aRes[ip*3+1] += r.y();
      aRes[ip*3+2] += r.z();
      aRes[jp*3+0] -= r.x();
      aRes[jp*3+1] -= r.y();
      aRes[jp*3+2] -= r.z();
    }
  }
}

void JArray_Extend
(std::vector<unsigned int>& psup_ind1,
 std::vector<unsigned int>& psup1,
 const std::vector<unsigned int>& psup_ind0,
 const std::vector<unsigned int>& psup0)
{
  
}

void ddW_ArapEnergy
 (dfm2::CMatrixSparse<double>& Mat,
  const std::vector<unsigned int>& psup_ind,
  const std::vector<unsigned int>& psup)
{
  const unsigned int np = psup_ind.size()-1;
  Mat.Initialize(np, 3, true);
  Mat.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
  Mat.SetZero();
  std::vector<int> tmp_buffer;
  for(unsigned int ip=0;ip<np;++ip){
    /*
    const dfm2::CVec3d Pi(aXYZ0.data()+ip*3);
    const dfm2::CQuatd Qi(aQuat1.data()+ip*4);
    dfm2::CMat3d LM; LM.SetZero();
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      const unsigned int jp = psup[ipsup];
      const dfm2::CVec3d v0 = Qi*(dfm2::CVec3d(aXYZ0.data()+jp*3)-Pi);
      LM += dfm2::Mat3_CrossCross(v0);
    }
    dfm2::CMat3d LMi = LM.Inverse();
     */
    for(unsigned int ipsup=psup_ind[ip];ipsup<psup_ind[ip+1];++ipsup){
      unsigned int jp = psup[ipsup];
//      const dfm2::CVec3d v0 = Qi*(dfm2::CVec3d(aXYZ0.data()+jp*3)-Pi);
      //        const dfm2::CVec3d v1 = dfm2::CVec3d(aXYZ1.data()+jp*3) - dfm2::CVec3d(aXYZ1.data()+ip*3);
      //         dfm2::CMat3d L1 = -dfm2::CMat3d::Spin(v0.p)*LMi*dfm2::CMat3d::Spin(v0.p);
      //        dfm2::CMat3d L1 = -dfm2::CMat3d::Spin(v0.p)*LMi;
      dfm2::CMat3d L1 = dfm2::CMat3d::Identity();
      double eMat[2][2][9];
      L1.CopyToScale(eMat[0][0],+1);
      L1.CopyToScale(eMat[1][1],+1);
      L1.CopyToScale(eMat[0][1],-1);
      L1.CopyToScale(eMat[1][0],-1);
      const unsigned int aIP[2] = {ip,jp};
      Mat.Mearge(2, aIP, 2, aIP, 9, &eMat[0][0][0], tmp_buffer);
    }
  }
}

void Check_ArapEnergy(const std::vector<double>& aXYZ0,
                      const std::vector<double>& aXYZ1,
                      const std::vector<double>& aQuat0,
                      const std::vector<double>& dXYZ12,
                      const std::vector<unsigned int>& psup_ind,
                      const std::vector<unsigned int>& psup)
{
  std::vector<double> aQuat1 = aQuat0;
  for(int itr=0;itr<20;++itr){
    UpdateRotationsByClusterMatching(aQuat1,
                                     aXYZ0,aXYZ1,psup_ind,psup);
  }
  
  double w1 = W_ArapEnergy(aXYZ0, aXYZ1, aQuat1, psup_ind, psup);
  std::vector<double> aRes1; dW_ArapEnergy(aRes1,
                                           aXYZ0, aXYZ1, aQuat1, psup_ind, psup);
  
  const double eps = 1.0e-5;
  std::vector<double> aXYZ2 = aXYZ1;
  for(int i=0;i<aXYZ2.size();++i){ aXYZ2[i] += eps*dXYZ12[i]; }
  
  std::vector<double> aQuat2 = aQuat1;
  UpdateRotationsByClusterMatching(aQuat2,
                                   aXYZ0, aXYZ2, psup_ind, psup);

  double w2 = W_ArapEnergy(aXYZ0, aXYZ2, aQuat2, psup_ind, psup);
  std::vector<double> aRes2; dW_ArapEnergy(aRes2,
                                           aXYZ0, aXYZ2, aQuat2, psup_ind, psup);

  {
    double dw = dfm2::Dot(aRes1,dXYZ12);
    double val2 = (w2-w1)/eps;
    std::cout << "       diff_sns: " << fabs(dw-val2)/(fabs(dw)+1.0) << std::endl;
  }

  /*
  {
    const unsigned int np = psup_ind.size()-1;
    std::vector<double> aRes12(np*3);
    Mat.MatVec(aRes12.data(),
               1.0, dXYZ12.data(), 0.0);
    for(int i=0;i<np*3;++i){
      std::cout << i << " " << (aRes2[i]-aRes1[i])/eps << " " << aRes12[i] << std::endl;
    }
    std::cout << "sym: " << dfm2::CheckSymmetry(Mat) << std::endl;
  }
   */
        
}

// --------------------

int main()
{
  std::vector<double> aXYZ0;
  std::vector<double> aW;
  std::vector<unsigned int> aTri;
  std::vector<dfm2::CRigBone> aBone;
  {
    std::vector<int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl(aXYZ0,
                         aW,
                         aTri,
                         aIndBoneParent,
                         aJntRgrs,
                         std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    Smpl2Rig(aBone,
             aIndBoneParent, aXYZ0, aJntRgrs);
    
  }
  std::vector<unsigned int> psup_ind, psup;
  {
    dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                               aTri.data(), aTri.size()/3, 3,
                               (int)aXYZ0.size()/3);
    dfm2::JArray_Sort(psup_ind, psup);
  }
  std::cout << psup.size() << std::endl;
  
  std::vector<double> aXYZ1 = aXYZ0;
  for(int ibone=0;ibone<aBone.size();++ibone){
    dfm2::CQuatd::Random(0.2).CopyTo(aBone[ibone].quatRelativeRot);
  }
  dfm2::UpdateBoneRotTrans(aBone);
  
  dfm2::Skinning_LBS(aXYZ1,
                     aXYZ0, aBone, aW);
  
  std::vector<double> aQuat1;
  {
    const unsigned int np = aXYZ1.size()/3;
    aQuat1.resize(np*4);
    for(int ip=0;ip<np;++ip){
      dfm2::Quat_Identity(aQuat1.data()+4*ip);
    }
    for(int itr=0;itr<20;++itr){
      UpdateRotationsByClusterMatching(aQuat1,
                                       aXYZ0,aXYZ1,psup_ind,psup);
    }
  }
  
  Check_SensitivityRigSkin(aBone,aXYZ0,aW);
  Check_SensitivityRigSkin(aBone,aXYZ0,aW);
  
  for(int itr=0;itr<5;++itr){
    unsigned int ib_s = aBone.size()*(rand()/(RAND_MAX+1.0));
    std::cout << "check sensitivity: " << ib_s << std::endl;
    unsigned int idim_s = 0;
    std::vector<double> dXYZ12;
    SensitivityRigSkin(dXYZ12,
                       ib_s, idim_s,
                       aBone, aXYZ0, aW);
    Check_ArapEnergy(aXYZ0,aXYZ1,aQuat1,
                     dXYZ12,psup_ind,psup);
  }
  

  
  /*
  {
    std::cout << psup.size() << " " << aXYZ0.size()/3 << std::endl;
    std::vector<double> dDiffdJoint; // nEdge*3*nJoint*3
    dDiffdJoint.resize(psup.size()*3*aBone.size()*3);
    time_t t0 = clock();
    std::cout << "hoge0" << std::endl;
    for(int ibone=0;ibone<aBone.size();++ibone){
      for(int irotbone=0;irotbone<3;++irotbone){
        std::vector<double> dXYZ12;
        SensitivityRigSkin(dXYZ12,
                           ibone, irotbone,
                           aBone, aXYZ0, aW);
        std::vector<double> dDiff12;
        dW_ArapEdgeDiff(dDiff12,
                        aXYZ0, aXYZ1, aQuat1, dXYZ12, psup_ind, psup);
      }
    }
    time_t t1 = clock();
    std::cout << "hoge1: " << (t1-t0)/1000000.0 << std::endl;
  }
   */
    
  // -----------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  dfm2::opengl::setSomeLighting();

  int iframe = 0;
  while (true)
  {
    iframe = (iframe+1)%50;
    if( iframe ==0 ){
      dfm2::UpdateBoneRotTrans(aBone);
      Skinning_LBS(aXYZ1,
                   aXYZ0, aBone, aW);
    }
    
    // -------------------
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    ::glEnable(GL_DEPTH_TEST);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), aTri.data(), aTri.size()/3);
    ::glDisable(GL_DEPTH_TEST);
    delfem2::opengl::DrawBone(aBone,
                              -1, -1,
                              0.01, 1.0);
    ::glEnable(GL_DEPTH_TEST);
    dfm2::opengl::Draw_QuaternionsCoordinateAxes(aXYZ1,aQuat1,0.02);
//    dfm2::opengl::DrawJoints(aJntPos1, aIndBoneParent);
//    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
//    dfm2::opengl::DrawJoints(aJntPos0, aIndBoneParent);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
