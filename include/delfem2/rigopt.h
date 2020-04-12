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

namespace delfem2 {

void Solve_MinRigging
 (std::vector<CRigBone>& aBone,
  const std::vector<CTarget>& aTarget)
{
  std::vector<double> Lx, Ly, Lz; // [nsns, nb*4]
  for(int ibs=0;ibs<aBone.size();++ibs){
    for(int idims=0;idims<3;++idims){
       Rig_SensitivityBoneTransform_Eigen(Lx,Ly,Lz,
                                               ibs,idims,true,
                                               aBone);
    }
  }
  for(int idims=0;idims<3;++idims){
     Rig_SensitivityBoneTransform_Eigen(Lx,Ly,Lz,
                                             0,idims,false,
                                             aBone);
  }
  
  // -----------------------
  
  std::vector<double> aC0; // [nC]
  std::vector<double> adC0; // [nC, nb*3 ]
  for(int it=0;it<aTarget.size();++it){
     Rig_WdW_Target_Eigen(aC0,adC0,
                               aBone,aTarget[it],Lx,Ly,Lz);
  }
  
  const unsigned int nsns = Lx.size()/(aBone.size()*4);
  const unsigned int nC = aC0.size();
  
  class CSystemMatrix{
  public:
    CSystemMatrix(const std::vector<double>& adC_,
                  unsigned int nC_,
                  unsigned int nsns_) :
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
                         adC.data(), nC, nsns, x);
       MatTVec(y,
                    adC.data(), nC, nsns, tmpC0.data());
      for(int i=0;i<nsns;++i){ y[i] += (beta+0.01)*x[i]; }
    }
  public:
    const std::vector<double>& adC;
    unsigned int nC;
    unsigned int nsns;
    mutable std::vector<double> tmpC0;
  } mat(adC0, nC, nsns);
  
  std::vector<double> r(nsns,0.0);
   MatTVec(r.data(),
                adC0.data(), nC, nsns, aC0.data());
  
  std::vector<double> u(nsns,0.0);
  std::vector<double> reshist =  Solve_CG(r.data(), u.data(),
                                               nsns, 1.0e-3, 100, mat);
  //  std::cout << "convergence" << reshist.size() << std::endl;
  for(int ib=0;ib<aBone.size();++ib){
     CVec3d vec_rot(u.data()+ib*3);
      CQuatd dq; Quat_CartesianAngle(dq.q, (-vec_rot).p);
     CQuatd q0 = dq* CQuatd(aBone[ib].quatRelativeRot);
    q0.CopyTo(aBone[ib].quatRelativeRot);
  }
  {
     CVec3d vec_trans(u.data()+aBone.size()*3);
    aBone[0].transRelative[0] -= vec_trans.x();
    aBone[0].transRelative[1] -= vec_trans.y();
    aBone[0].transRelative[2] -= vec_trans.z();
  }
   UpdateBoneRotTrans(aBone);
}

}


#endif /* rigoptim_h */
