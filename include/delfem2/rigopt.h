//
//  rigoptim.h
//  examples_glfwold_static
//
//  Created by Nobuyuki Umetani on 2020-04-06.
//

#ifndef DFM2_RIGOPT_H
#define DFM2_RIGOPT_H

#include "delfem2/rig_geo3.h"
#include "delfem2/quat.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/lsvecx.h"

namespace delfem2 {


/**
 * How bone's affine matrix will change due to the imput of the articulated bones
 * @param[out,in] add sensitivity of bone's 0-th row of afffine matrix [ *, nBone, 12 ]
 * @param[in] ib_s changing bone
 * @param[in] idim_s changing bone
 * @param[in] is_rot is bone change is rotation or translation
 * @param[in] aBone1 current sleleton
 */
DFM2_INLINE void Rig_SensitivityBoneTransform(
    std::vector<double>& L, // [ nsns, nBone*12 ]
    unsigned int ib_s,
    unsigned int idim_s,
    bool is_rot,
    const std::vector<CRigBone>& aBone1)
{
  const size_t nb = aBone1.size();
  const unsigned int istat0 = static_cast<unsigned int>(L.size());
  assert( istat0%12 == 0 );
  L.resize(istat0+nb*12);
  std::vector<double> aM(nb*16);

  for(std::size_t ibone=0;ibone<aBone1.size();++ibone){
    CMat4d m01 = CMat4d::Scale(aBone1[ibone].scale);
    m01 = CMat4d::Quat(aBone1[ibone].quatRelativeRot) * m01;
    if( ibone == ib_s && is_rot ){
      CMat3d dn1 = CMat3d::Spin(CVec3d::Axis(idim_s).p) + CMat3d::Identity();
      CMat4d dm1 = CMat4d::Mat3(dn1.data());
      m01 = dm1 * m01;
    }
    m01 = CMat4d::Translate(aBone1[ibone].transRelative) * m01;
    if( ibone == ib_s && !is_rot ){
      m01 = CMat4d::Translate(CVec3d::Axis(idim_s).p) * m01;
    }
    const int ibone_p = aBone1[ibone].ibone_parent;
    if( ibone_p  == -1 ){ // root bone
      Copy_Mat4( aM.data()+ibone*16, m01.data() );
      continue;
    }
    MatMat4(aM.data()+ibone*16,
            aM.data()+ibone_p*16, m01.data() );
  }
  for(std::size_t ib=0;ib<nb;++ib){
    for(int jdim=0;jdim<4;++jdim){
      L[istat0+ib*12+4*0+jdim] = aM[ib*16+0*4+jdim] - aBone1[ib].affmat3Global[0*4+jdim];
      L[istat0+ib*12+4*1+jdim] = aM[ib*16+1*4+jdim] - aBone1[ib].affmat3Global[1*4+jdim];
      L[istat0+ib*12+4*2+jdim] = aM[ib*16+2*4+jdim] - aBone1[ib].affmat3Global[2*4+jdim];
    }
  }

}

DFM2_INLINE void
Rig_Sensitivity_Skeleton(
    std::vector<double>& L,
    const std::vector<CRigBone>& aBone)
{
  L.resize(0);
  for(unsigned int ibs=0;ibs<aBone.size();++ibs){ // add rotation of bone
    for(int idims=0;idims<3;++idims){
      Rig_SensitivityBoneTransform(
          L,
          ibs,idims,true,
          aBone);
    }
  }
  for(int idims=0;idims<3;++idims){ // add translation of root bone
    Rig_SensitivityBoneTransform(
        L,
        0,idims,false,
        aBone);
  }
}

void Sensitivity_RigSkinPoint(
    std::vector<double>& dP,
    unsigned int ip,
    const double* p0,
    const std::vector<CRigBone>& aBone,
    const std::vector<double>& L,  // [ nsns, nb*12 ]
    unsigned int nb_par_pt, // number of bone par pointk
    const std::vector<double>& aSkinSparseW,
    const std::vector<unsigned int>& aSkinSparseI)
{
  assert( aSkinSparseI.size() == aSkinSparseW.size() );
  const size_t nb = aBone.size();
  const size_t nsns = L.size()/(nb*12);
  dP.assign(3*nsns, 0.0);
  const double p0a[4] = { p0[0], p0[1], p0[2], 1.0};
  for (unsigned int jjb = 0; jjb < nb_par_pt; ++jjb) {
    const unsigned int jb0 = aSkinSparseI[ip*nb_par_pt+jjb];
    const double wb0 = aSkinSparseW[ip*nb_par_pt+jjb];
    double p0b[4];
    delfem2::MatVec4(p0b,
                     aBone[jb0].invBindMat, p0a);
    for(unsigned int isns=0;isns<nsns;++isns) {
      for (unsigned int jdim = 0; jdim < 4; ++jdim) {
        dP[0*nsns+isns] += wb0 * p0b[jdim] * L[isns * (nb * 12) + jb0 * 12 + 4 * 0 + jdim];
        dP[1*nsns+isns] += wb0 * p0b[jdim] * L[isns * (nb * 12) + jb0 * 12 + 4 * 1 + jdim];
        dP[2*nsns+isns] += wb0 * p0b[jdim] * L[isns * (nb * 12) + jb0 * 12 + 4 * 2 + jdim];
      }
    }
  }
}

class CTarget {
public:
  unsigned int ib;
  CVec3d pos;
public:
};

DFM2_INLINE void Rig_WdW_Target(
    std::vector<double>& aW,
    std::vector<double>& adW,
    const std::vector<CRigBone>& aBone,
    const CTarget& target,
    const std::vector<double>& L) // [ nsns, nBone, 12]
{
  const size_t nb = aBone.size();
  const size_t nsns = L.size()/(nb*12);
  assert( L.size() == nsns*nb*12 );
  // --------
  unsigned int ib = target.ib;
  const CVec3d pos = target.pos;
  // ---------------
  const CVec3d p0 = aBone[ib].Pos();
  const unsigned int ncnst = 2;
  {
    double sqx = pos.x-p0.x;
    double sqy = pos.y-p0.y;
    aW.push_back(sqx);
    aW.push_back(sqy);
  }
  const unsigned int istat = static_cast<unsigned int>(adW.size());
  adW.resize(istat+ncnst*nsns);
  for(unsigned int isns=0;isns<nsns;++isns){
    double dx = L[isns*(nb*12) + ib*12+3];
    double dy = L[isns*(nb*12) + ib*12+7];
    adW[istat+0*nsns+isns] = -dx;
    adW[istat+1*nsns+isns] = -dy;
  }
}



DFM2_INLINE void
Solve_MinRigging(
    std::vector<CRigBone>& aBone,
    const std::vector<CTarget>& aTarget)
{
  std::vector<double> L; // [nsns, nb*12]
  Rig_Sensitivity_Skeleton(
      L,
      aBone);
  assert( L.size() == (aBone.size()+1)*3*aBone.size()*12 );
  
  // -----------------------
  
  std::vector<double> aC0; // [nC]
  std::vector<double> adC0; // [nC, nb*3 ]
  for(unsigned int it=0;it<aTarget.size();++it){
     Rig_WdW_Target(aC0,adC0,
         aBone,aTarget[it],L);
  }
  
  const size_t nsns = L.size()/(aBone.size()*12);
  const size_t nC = aC0.size();
  
  class CSystemMatrix{
  public:
    CSystemMatrix(const std::vector<double>& adC_,
                  size_t nC_,
                  size_t nsns_) :
    adC(adC_), nC(nC_), nsns(nsns_)
    {
      //      std::cout << "constructor reduced system matrix " << std::endl;
      assert(adC.size()==nsns*nC);
      tmpC0.resize(nC_);
    }
  public:
    void MatVec(double* y,
                double alpha, const double* x, double beta) const {
       ::delfem2::MatVec(tmpC0.data(),                         
		   adC.data(), 
		   static_cast<unsigned int>(nC), 
		   static_cast<unsigned int>(nsns), 
		   x);
       MatTVec(y,
		   adC.data(),
		   static_cast<unsigned int>(nC), 
		   static_cast<unsigned int>(nsns), 
		   tmpC0.data());
      for(unsigned int i=0;i<nsns;++i){ y[i] += (beta+0.01)*x[i]; }
    }
  public:
    const std::vector<double>& adC;
    size_t nC;
    size_t nsns;
    mutable std::vector<double> tmpC0;
  } mat(adC0, nC, nsns);
  
  std::vector<double> r(nsns,0.0);
   MatTVec(r.data(),
	   adC0.data(), 
	   static_cast<unsigned int>(nC), 
	   static_cast<unsigned int>(nsns), 
	   aC0.data());
  
  std::vector<double> u(nsns,0.0);
  std::vector<double> reshist;
  {
    const std::size_t n = nsns;
    std::vector<double> tmp0(n), tmp1(n);
    auto vr = delfem2::CVecXd(r);
    auto vu = delfem2::CVecXd(u);
    auto vs = delfem2::CVecXd(tmp0);
    auto vt = delfem2::CVecXd(tmp1);
    reshist = Solve_CG(
        vr, vu, vs, vt,
        1.0e-3, 100, mat);
  }
  //  std::cout << "convergence" << reshist.size() << std::endl;
  for(unsigned int ib=0;ib<aBone.size();++ib){
     CVec3d vec_rot(u.data()+ib*3);
      CQuatd dq; Quat_CartesianAngle(dq.p, (-vec_rot).p);
     CQuatd q0 = dq* CQuatd(aBone[ib].quatRelativeRot);
    q0.CopyTo(aBone[ib].quatRelativeRot);
  }
  {
     CVec3d vec_trans(u.data()+aBone.size()*3);
    aBone[0].transRelative[0] -= vec_trans.x;
    aBone[0].transRelative[1] -= vec_trans.y;
    aBone[0].transRelative[2] -= vec_trans.z;
  }
  UpdateBoneRotTrans(aBone);
}

}


#endif /* rigoptim_h */
