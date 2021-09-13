/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/femrod.h"

#include "delfem2/fem_distance3.h"
#include "delfem2/fem_rod3_darboux.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/lsvecx.h"
#include "delfem2/mat3.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/mshuni.h"
#include "delfem2/jagarray.h"

// ========================================
// RodHair

DFM2_INLINE void delfem2::ParallelTransport_RodHair(
    std::vector<CVec3d> &aP0,
    std::vector<CVec3d> &aS0,
    const std::vector<unsigned int> &aIP_HairRoot) {
  assert(aP0.size() == aS0.size());
  assert(!aIP_HairRoot.empty() && aIP_HairRoot[0] == 0);
  assert(aP0.size() == aIP_HairRoot[aIP_HairRoot.size() - 1]);
  for (unsigned int ih = 0; ih < aIP_HairRoot.size() - 1; ++ih) {
    const unsigned int ip_r = aIP_HairRoot[ih];
    const unsigned int np = aIP_HairRoot[ih + 1] - ip_r;
    for (unsigned int ir = 0; ir < np - 2; ++ir) {
      const unsigned int ip0 = ip_r + ir + 0;
      const unsigned int ip1 = ip_r + ir + 1;
      const unsigned int ip2 = ip_r + ir + 2;
      const unsigned int is0 = ip_r + ir + 0;
      const unsigned int is1 = ip_r + ir + 1;
      const CMat3d CMat3 = Mat3_MinimumRotation(aP0[ip1] - aP0[ip0], aP0[ip2] - aP0[ip1]);
      CVec3d s1 = CMat3 * aS0[is0] + aS0[is1];
      const CVec3d v = (aP0[ip2] - aP0[ip1]).normalized();
      aS0[is1] = (s1 - (s1.dot(v)) * v).normalized();
    }
  }
}

DFM2_INLINE void delfem2::MakeBCFlag_RodHair(
    std::vector<int> &aBCFlag,
    const std::vector<unsigned int> &aIP_HairRoot) {
  assert(!aIP_HairRoot.empty() && aIP_HairRoot[0] == 0);
  const unsigned int np = aIP_HairRoot[aIP_HairRoot.size() - 1];
  aBCFlag.assign(np * 4, 0);
  for (unsigned int ihair = 0; ihair < aIP_HairRoot.size() - 1; ++ihair) {
    assert(aIP_HairRoot[ihair + 1] > aIP_HairRoot[ihair]);
    unsigned int ips0 = aIP_HairRoot[ihair] + 0;
    unsigned int ips1 = aIP_HairRoot[ihair] + 1;
    unsigned int ipe1 = aIP_HairRoot[ihair + 1] - 1;
    aBCFlag[ips0 * 4 + 0] = 1;
    aBCFlag[ips0 * 4 + 1] = 1;
    aBCFlag[ips0 * 4 + 2] = 1;
    aBCFlag[ips1 * 4 + 0] = 1;
    aBCFlag[ips1 * 4 + 1] = 1;
    aBCFlag[ips1 * 4 + 2] = 1;
    aBCFlag[ips0 * 4 + 3] = 1;
    aBCFlag[ipe1 * 4 + 3] = 1;
  }
}

DFM2_INLINE void delfem2::MakeSparseMatrix_RodHair(
    CMatrixSparse<double> &mats,
    const std::vector<unsigned int> &aIP_HairRoot,
    unsigned int ndof_par_node) {
  assert(!aIP_HairRoot.empty() && aIP_HairRoot[0] == 0);
  const unsigned int np = aIP_HairRoot[aIP_HairRoot.size() - 1];
  std::vector<unsigned int> psup_ind, psup;
  {
    std::vector<unsigned int> aElemRod;
    for (unsigned int ihair = 0; ihair < aIP_HairRoot.size() - 1; ++ihair) {
      unsigned int ip0 = aIP_HairRoot[ihair];
      unsigned int nr = aIP_HairRoot[ihair + 1] - ip0 - 2; // number of rod elements
      for (unsigned int ir = 0; ir < nr; ++ir) {
        aElemRod.push_back(ip0 + ir + 0);
        aElemRod.push_back(ip0 + ir + 1);
        aElemRod.push_back(ip0 + ir + 2);
      }
    }
    JArray_PSuP_MeshElem(
        psup_ind, psup,
        aElemRod.data(), aElemRod.size() / 3, 3, np);
  }
  JArray_Sort(psup_ind, psup);
  mats.Initialize(np, ndof_par_node, true);
  mats.SetPattern(psup_ind.data(), psup_ind.size(), psup.data(), psup.size());
}

DFM2_INLINE void delfem2::MakeDirectorOrthogonal_RodHair(
    std::vector<CVec3d> &aS,
    const std::vector<CVec3d> &aP) {
  for (unsigned int is = 0; is < aP.size() - 1; ++is) {
    assert(is < aS.size());
    const unsigned int ip0 = is + 0;
    const unsigned int ip1 = is + 1;
    const CVec3d &p0 = aP[ip0];
    const CVec3d &p1 = aP[ip1];
    const CVec3d e01 = (p1 - p0).normalized();
    aS[is] -= (aS[is].dot(e01)) * e01;
    aS[is].normalize();
  }
  /*
  for(unsigned int ihair=0;ihair<aIP_HairRoot.size()-1;++ihair){
    unsigned int ips = aIP_HairRoot[ihair];
    unsigned int ns = aIP_HairRoot[ihair+1] - aIP_HairRoot[ihair] -1;
    for(unsigned int is=0;is<ns;++is){
      const unsigned int ip0 = ips+is+0;
      const unsigned int ip1 = ips+is+1;
      const CVec3d& p0 = aP[ip0];
      const CVec3d& p1 = aP[ip1];
      const CVec3d e01 = (p1-p0).Normalize();
      aS[ip0] -= (aS[ip0]*e01)*e01;
      aS[ip0].SetNormalizedVector();
    }
  }
  */
}

double delfem2::MergeLinSys_Hair(
    std::vector<double> &vec_r,
    delfem2::CMatrixSparse<double> &mats,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    const std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<delfem2::CVec3d> &aP,
    const std::vector<delfem2::CVec3d> &aS,
    const std::vector<delfem2::CVec3d> &aP0,
    const std::vector<delfem2::CVec3d> &aS0) {
  using namespace delfem2;
  std::vector<unsigned int> tmp_buffer;
  double W = 0;
  for (unsigned int ihair = 0; ihair < aIP_HairRoot.size() - 1; ++ihair) {
    unsigned int ips = aIP_HairRoot[ihair];
    unsigned int ns = aIP_HairRoot[ihair + 1] - aIP_HairRoot[ihair] - 1;
    for (unsigned int is = 0; is < ns; ++is) {
      const unsigned int ip0 = ips + is + 0;
      const unsigned int ip1 = ips + is + 1;
      const unsigned int aINoel[2] = {ip0, ip1};
      const double L0 = (aP0[ip0] - aP0[ip1]).norm();
      const CVec3d aPE[2] = {aP[ip0], aP[ip1]};
      // --------------
      CVec3d dW_dP[2];
      CMat3d ddW_ddP[2][2];
      W += WdWddW_SquareLengthLineseg3D(
          dW_dP, ddW_ddP,
          stiff_stretch, aPE, L0);
      {
        double eM[2][2][4][4];
        for (unsigned int i = 0; i < 2 * 2 * 4 * 4; ++i) { (&eM[0][0][0][0])[i] = 0.0; }
        for (int in = 0; in < 2; ++in) {
          for (int jn = 0; jn < 2; ++jn) {
            ddW_ddP[in][jn].CopyToMat4(&eM[in][jn][0][0]);
          }
        }
        Merge<2, 2, 4, 4, double>(mats, aINoel, aINoel, eM, tmp_buffer);
      }
      for (int in = 0; in < 2; in++) {
        const unsigned int ip = aINoel[in];
        vec_r[ip * 4 + 0] -= dW_dP[in].x;
        vec_r[ip * 4 + 1] -= dW_dP[in].y;
        vec_r[ip * 4 + 2] -= dW_dP[in].z;
      }
    }
  }
  // --------------------------
  for (unsigned int ihair = 0; ihair < aIP_HairRoot.size() - 1; ++ihair) {
    const unsigned int ips = aIP_HairRoot[ihair];
    const unsigned int nr = aIP_HairRoot[ihair + 1] - aIP_HairRoot[ihair] - 2;
    for (unsigned int ir = 0; ir < nr; ++ir) {
      const unsigned int ip0 = ips + ir + 0;
      const unsigned int ip1 = ips + ir + 1;
      const unsigned int ip2 = ips + ir + 2;
      const unsigned int aINoel[3] = {ip0, ip1, ip2};
      const CVec3d aPE[3] = {aP[ip0], aP[ip1], aP[ip2]};
      const CVec3d aSE[2] = {aS[ip0], aS[ip1]};
      CVec3d Darboux0;
      { // Daboux vector of initial configuration
        const CVec3d aPE0[3] = {aP0[ip0], aP0[ip1], aP0[ip2]};
        const CVec3d aSE0[2] = {aS0[ip0], aS0[ip1]};
        Darboux0 = Darboux_Rod(aPE0, aSE0);
      }
      // ------
      CVec3d dW_dP[3];
      double dW_dt[2];
      CMat3d ddW_ddP[3][3];
      CVec3d ddW_dtdP[2][3];
      double ddW_ddt[2][2];
      W += WdWddW_Rod(
          dW_dP, dW_dt, ddW_ddP, ddW_dtdP, ddW_ddt,
          stiff_bendtwist,
          aPE, aSE, Darboux0, false);
      {
        double eM[3][3][4][4];
        for (unsigned int i = 0; i < 3 * 3 * 4 * 4; ++i) { (&eM[0][0][0][0])[i] = 0.0; }
        for (int in = 0; in < 3; ++in) {
          for (int jn = 0; jn < 3; ++jn) {
            ddW_ddP[in][jn].CopyToMat4(&eM[in][jn][0][0]);
          }
        }
        for (int in = 0; in < 3; ++in) {
          for (int jn = 0; jn < 2; ++jn) {
            eM[in][jn][0][3] = eM[jn][in][3][0] = ddW_dtdP[jn][in].x;
            eM[in][jn][1][3] = eM[jn][in][3][1] = ddW_dtdP[jn][in].y;
            eM[in][jn][2][3] = eM[jn][in][3][2] = ddW_dtdP[jn][in].z;
          }
        }
        for (int in = 0; in < 2; ++in) {
          for (int jn = 0; jn < 2; ++jn) {
            eM[in][jn][3][3] = ddW_ddt[in][jn];
          }
        }
        Merge<3, 3, 4, 4, double>(mats, aINoel, aINoel, eM, tmp_buffer);
      }
      {
        for (int ino = 0; ino < 3; ino++) {
          const unsigned int ip = aINoel[ino];
          vec_r[ip * 4 + 0] -= dW_dP[ino].x;
          vec_r[ip * 4 + 1] -= dW_dP[ino].y;
          vec_r[ip * 4 + 2] -= dW_dP[ino].z;
        }
        for (int in = 0; in < 2; in++) {
          const unsigned int in0 = aINoel[in];
          vec_r[in0 * 4 + 3] -= dW_dt[in];
        }
      }
    }
  }
  return W;
}

DFM2_INLINE void delfem2::UpdateSolutionHair(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    const std::vector<double> &vec_x,
    const std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<int> &aBCFlag) {
  for (unsigned int ihair = 0; ihair < aIP_HairRoot.size() - 1; ++ihair) {
    unsigned int ips = aIP_HairRoot[ihair];
    unsigned int ns = aIP_HairRoot[ihair + 1] - aIP_HairRoot[ihair] - 1;
    for (unsigned int is = 0; is < ns; ++is) {
      const unsigned int ip0 = ips + is + 0;
      const unsigned int ip1 = ips + is + 1;
      CVec3d V01 = aP[ip1] - aP[ip0];
      CVec3d du(vec_x[ip1 * 4 + 0] - vec_x[ip0 * 4 + 0],
                vec_x[ip1 * 4 + 1] - vec_x[ip0 * 4 + 1],
                vec_x[ip1 * 4 + 2] - vec_x[ip0 * 4 + 2]);
      const double dtheta = vec_x[ip0 * 4 + 3];
      CVec3d frm[3];
      RodFrameTrans(frm,
                    aS[ip0], V01, du, dtheta);
      aS[ip0] = frm[0];
    }
  }
  for (unsigned int ip = 0; ip < aP.size(); ++ip) {
    if (aBCFlag[ip * 4 + 0] != 0) continue;
    aP[ip].p[0] += vec_x[ip * 4 + 0];
    aP[ip].p[1] += vec_x[ip * 4 + 1];
    aP[ip].p[2] += vec_x[ip * 4 + 2];
  }
}

DFM2_INLINE void delfem2::Solve_RodHair(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    CMatrixSparse<double> &mats,
    const double stiff_stretch,
    const double stiff_bendtwist[3],
    double mdtt,
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aS0,
    const std::vector<int> &aBCFlag,
    const std::vector<unsigned int> &aIP_HairRoot) {
  assert(mats.nrowdim_ == 4 && mats.ncoldim_ == 4);
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
  std::cout << "energy:" << W << std::endl;
  //    std::cout << "sym: " << CheckSymmetry(mats) << std::endl;
  assert(aBCFlag.size() == np * 4);
  mats.SetFixedBC(aBCFlag.data());
  setRHS_Zero(vec_r, aBCFlag, 0);
  std::vector<double> vec_x;
  vec_x.assign(np * 4, 0.0);
  {
    const std::size_t n = vec_r.size();
    std::vector<double> tmp0(n), tmp1(n);
    auto aConvHist = Solve_CG(
        CVecXd(vec_r), CVecXd(vec_x), CVecXd(tmp0), CVecXd(tmp1),
        1.0e-4, 300, mats);
    if (aConvHist.size() > 0) {
      std::cout << "            conv: " << aConvHist.size() << " " << aConvHist[0] << " "
                << aConvHist[aConvHist.size() - 1] << std::endl;
    }
  }
  UpdateSolutionHair(aP, aS,
                     vec_x, aIP_HairRoot, aBCFlag);
}

namespace delfem2 {
namespace femrod {

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
  femrod::CMatContact mc(mats, aContact, stiff_contact);
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

