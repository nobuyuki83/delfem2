/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/hair_darboux_collision.h"

#include "delfem2/hair_darboux_util.h"
#include "delfem2/hair_darboux_solver.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/view_vectorx.h"
#include "delfem2/mat3.h"
#include "delfem2/mat3_funcs.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/mshuni.h"

namespace delfem2::hair_darboux_collision {

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
      const std::array<double,9> NN = Mat3_OuterProduct(nrm, nrm);
      for (int iip = 0; iip < 4; ++iip) {
        unsigned int ip0 = aIP[iip];
        for (int jjp = 0; jjp < 4; ++jjp) {
          unsigned int jp0 = aIP[jjp];
          double v0[3];
          MatVec3(v0, NN.data(), x + jp0 * 4);
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

DFM2_INLINE void delfem2::Solve_RodHairContact(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    delfem2::LinearSystemSolver_BlockSparse &mats,
    double stiff_stretch,
    const double stiff_bendtwist[3],
    double mdtt,
    const std::vector<CVec3d> &aPt0,
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aS0,
    const std::vector<unsigned int> &aIP_HairRoot,
    double clearance,
    double stiff_contact,
    const std::vector<CContactHair> &aContact){
  namespace lcl = delfem2::hair_darboux_collision;
  assert(mats.ndim() == 4);
  const size_t np = aP.size();
  assert(aP0.size() == np);
  assert(aS0.size() == np);
  assert(aP.size() == np);
  assert(aS.size() == np);
  mats.BeginMerge();
  double W = 0.0;
  W += Merge_HairDarboux(
      mats,
      stiff_bendtwist,
      aIP_HairRoot, aP, aS, aP0, aS0);
  W += Merge_HairStretch(
      mats,
      stiff_stretch,
      aIP_HairRoot, aP, aP0);
  for (unsigned int ip = 0; ip < aP.size(); ++ip) {
    mats.AddValueToDiagonal(ip,0,mdtt);
    mats.AddValueToDiagonal(ip,1,mdtt);
    mats.AddValueToDiagonal(ip,2,mdtt);
  }
  for (unsigned int ip = 0; ip < aP.size(); ++ip) { // this term has effect only in the second nonliear iteration
    const CVec3d dp = aPt0[ip] - aP[ip];
    mats.vec_r[ip*4+0] += dp.x * mdtt;
    mats.vec_r[ip*4+1] += dp.y * mdtt;
    mats.vec_r[ip*4+2] += dp.z * mdtt;
  }
  for (const auto &ch : aContact) {
    const unsigned int aIP[4] = {ch.ip0, ch.ip1, ch.iq0, ch.iq1};
    const double aW[4] = {1 - ch.s, ch.s, -(1 - ch.t), -ch.t};
    CVec3d a = aW[0] * aP[aIP[0]] + aW[1] * aP[aIP[1]] + aW[2] * aP[aIP[2]] + aW[3] * aP[aIP[3]];
    double r0 = a.dot(ch.norm) - clearance;
    for (int iip = 0; iip < 4; ++iip) {
      const unsigned int ip0 = aIP[iip];
      mats.vec_r[ip0*4+0] -= stiff_contact * r0 * aW[iip] * ch.norm.x;
      mats.vec_r[ip0*4+1] -= stiff_contact * r0 * aW[iip] * ch.norm.y;
      mats.vec_r[ip0*4+2] -= stiff_contact * r0 * aW[iip] * ch.norm.z;
    }
  }
  // --------------
  W = 0.0; // to remove warning
  assert(mats.dof_bcflag.size() == np * 4);
  mats.matrix.SetFixedBC(mats.dof_bcflag.data());
  setRHS_Zero(mats.vec_r, mats.dof_bcflag, 0);
  lcl::CMatContact mc(mats.matrix, aContact, stiff_contact);
  {
    const std::size_t n = mats.ndof();
    mats.vec_x.assign(n, 0.0);
    mats.tmp0.resize(n);
    mats.tmp1.resize(n);
    auto vr = ViewAsVectorXd(mats.vec_r);
    auto vu = ViewAsVectorXd(mats.vec_x);
    auto vs = ViewAsVectorXd(mats.tmp0);
    auto vt = ViewAsVectorXd(mats.tmp1);
    auto aConvHist = Solve_CG(
        vr, vu, vs, vt,
        1.0e-6, 3000, mc);
  }
  UpdateSolutionHair(
      aP, aS,
      mats.vec_x, aIP_HairRoot, mats.dof_bcflag);
}
