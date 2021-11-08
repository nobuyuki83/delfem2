/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/hair_darboux.h"

#include "delfem2/fem_rod3_darboux.h"
#include "delfem2/geo3_v23m34q.h"

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


void delfem2::MakeProblemSetting_Spiral(
    std::vector<delfem2::CVec3d> &aP0,
    std::vector<delfem2::CVec3d> &aS0,
    std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<CHairShape> &aHairShape) {
  aIP_HairRoot.assign(1, 0);
  aP0.clear();
  aS0.clear();
  for (unsigned int ihair = 0; ihair < aHairShape.size(); ++ihair) {
    const unsigned int np = aHairShape[ihair].np;
    const double pitch = aHairShape[ihair].pitch;
    const double dangle = aHairShape[ihair].dangle;
    const double rad0 = aHairShape[ihair].rad0;
    const double *p0 = aHairShape[ihair].p0;
    for (unsigned int ip = 0; ip < np; ++ip) {
      CVec3d p = CVec3d(
          p0[0] + ip * pitch,
          p0[1] + rad0 * cos(dangle * ip),
          p0[2] + rad0 * sin(dangle * ip));
      aP0.push_back(p);
    }
    const unsigned int np0 = aIP_HairRoot[ihair];
    for (unsigned int is = 0; is < np - 1; ++is) {
      const CVec3d v = (aP0[np0 + is + 1] - aP0[np0 + is + 0]).normalized();
      CVec3d s(1.3, 1.5, 1.7);
      s = (s - (s.dot(v)) * v).normalized();
      aS0.push_back(s);
    }
    aS0.emplace_back(1, 0, 0);
    aIP_HairRoot.push_back(static_cast<unsigned int>(aP0.size()));
  }
}