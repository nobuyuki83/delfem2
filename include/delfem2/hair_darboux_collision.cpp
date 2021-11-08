/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/hair_darboux_collision.h"

#include "delfem2/hair_darboux.h"
#include "delfem2/hair_sparse.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/lsvecx.h"
#include "delfem2/mat3.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/mshuni.h"

namespace delfem2 {
namespace hair_darboux_collision {

class CMatContact {
 public:
  CMatContact(
      const CMatrixSparse<double> &mats_,
      const std::vector<CContactHair> &aContact_,
      const double stiffness_)
      : mats(mats_), aContact(aContact_), stiffness(stiffness_) {}
  void MatVec(double *y,
              double alpha, const double *x,
              double beta) const {
    mats.MatVec(y,
                alpha, x, beta);
    for (unsigned int ic = 0; ic < aContact.size(); ++ic) {
      const CContactHair &ch = aContact[ic];
      const unsigned int aIP[4] = {ch.ip0, ch.ip1, ch.iq0, ch.iq1};
      const double aW[4] = {1 - ch.s, ch.s, -(1 - ch.t), -ch.t};
      const CVec3d &nrm = ch.norm;
//        std::cout << "lennorm: " << nrm.Length() << std::endl;
      const auto NN = Mat3_OuterProduct(nrm, nrm);
      for (int iip = 0; iip < 4; ++iip) {
        unsigned int ip0 = aIP[iip];
        for (int jjp = 0; jjp < 4; ++jjp) {
          unsigned int jp0 = aIP[jjp];
          double v0[3];
          MatVec3(v0,
                  NN.p_, x + jp0 * 4);
          y[ip0 * 4 + 0] += alpha * stiffness * aW[iip] * aW[jjp] * v0[0];
          y[ip0 * 4 + 1] += alpha * stiffness * aW[iip] * aW[jjp] * v0[1];
          y[ip0 * 4 + 2] += alpha * stiffness * aW[iip] * aW[jjp] * v0[2];
        }
      }
    }
  }
 public:
  const CMatrixSparse<double> &mats;
  const std::vector<CContactHair> &aContact;
  const double stiffness;
};
}
}

DFM2_INLINE void delfem2::Solve_RodHairContact(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    CMatrixSparse<double> &mats,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    double mdtt,
    const std::vector<CVec3d> &aPt0,
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aS0,
    const std::vector<int> &aBCFlag,
    const std::vector<unsigned int> &aIP_HairRoot,
    const double clearance,
    const double stiff_contact,
    const std::vector<CContactHair> &aContact) {
  namespace lcl = ::delfem2::hair_darboux_collision;
  assert(mats.nrowdim_ == 4);
  assert(mats.ncoldim_ == 4);
  const size_t np = aP.size();
  assert(aP0.size() == np);
  assert(aS0.size() == np);
  assert(aP.size() == np);
  assert(aS.size() == np);
  mats.setZero();
  std::vector<double> vec_r;
  vec_r.assign(np * 4, 0.0);
  double W = MergeLinSys_Hair(
      vec_r, mats,
      stiff_stretch, stiff_bendtwist,
      aIP_HairRoot, aP, aS, aP0, aS0);
  for (unsigned int ip = 0; ip < aP.size(); ++ip) {
    mats.val_dia_[ip * 16 + 0 * 4 + 0] += mdtt;
    mats.val_dia_[ip * 16 + 1 * 4 + 1] += mdtt;
    mats.val_dia_[ip * 16 + 2 * 4 + 2] += mdtt;
  }
  for (unsigned int ip = 0; ip < aP.size(); ++ip) { // this term has effect only in the second nonliear iteration
    const CVec3d dp = aPt0[ip] - aP[ip];
    vec_r[ip * 4 + 0] += dp.x * mdtt;
    vec_r[ip * 4 + 1] += dp.y * mdtt;
    vec_r[ip * 4 + 2] += dp.z * mdtt;
  }
//  std::cout << "energy:" << W << std::endl;
  //    std::cout << "sym: " << CheckSymmetry(mats) << std::endl;
  for (const auto &ch : aContact) {
    const unsigned int aIP[4] = {ch.ip0, ch.ip1, ch.iq0, ch.iq1};
    const double aW[4] = {1 - ch.s, ch.s, -(1 - ch.t), -ch.t};
    CVec3d a = aW[0] * aP[aIP[0]] + aW[1] * aP[aIP[1]] + aW[2] * aP[aIP[2]] + aW[3] * aP[aIP[3]];
    double r0 = a.dot(ch.norm) - clearance;
//    std::cout << "    contact: " << r0 << std::endl;
    for (int iip = 0; iip < 4; ++iip) {
      const unsigned int ip0 = aIP[iip];
      vec_r[ip0 * 4 + 0] -= stiff_contact * r0 * aW[iip] * ch.norm.x;
      vec_r[ip0 * 4 + 1] -= stiff_contact * r0 * aW[iip] * ch.norm.y;
      vec_r[ip0 * 4 + 2] -= stiff_contact * r0 * aW[iip] * ch.norm.z;
    }
  }
  // --------------
  W = 0.0; // to remove warning
  assert(aBCFlag.size() == np * 4);
  mats.SetFixedBC(aBCFlag.data());
  std::vector<double> vec_x;
  vec_x.assign(np * 4, 0.0);
  setRHS_Zero(vec_r, aBCFlag, 0);
  lcl::CMatContact mc(mats, aContact, stiff_contact);
  {
    const std::size_t n = vec_r.size();
    std::vector<double> tmp0(n), tmp1(n);
    auto vr = CVecXd(vec_r);
    auto vu = CVecXd(vec_x);
    auto vs = CVecXd(tmp0);
    auto vt = CVecXd(tmp1);
    auto aConvHist = Solve_CG(
        vr, vu, vs, vt,
        1.0e-6, 3000, mc);
    /*
    if( aConvHist.size() > 0 ){
      std::cout << "            conv: " << aConvHist.size() << " " << aConvHist[0] << " " << aConvHist[aConvHist.size()-1] << std::endl;
    }
     */
//    std::cout << "  updates: " << DotX(vec_x.data(), vec_x.data(), vec_x.size()) << std::endl;
  }
  UpdateSolutionHair(aP, aS,
                     vec_x, aIP_HairRoot, aBCFlag);
  /*
  std::cout << "hogehoge " << aContact.size() << std::endl;
  for(const auto& ch : aContact){
    const unsigned int aIP[4] = {ch.ip0, ch.ip1, ch.iq0, ch.iq1};
    const double aW[4] = {1-ch.s, ch.s, -(1-ch.t), -ch.t};
    CVec3d a = aW[0]*aP[aIP[0]] + aW[1]*aP[aIP[1]] + aW[2]*aP[aIP[2]] + aW[3]*aP[aIP[3]];
    double r0 = a*ch.norm - clearance;
    std::cout << "       contact: " << r0 << " " << clearance << std::endl;
  }
   */
}

