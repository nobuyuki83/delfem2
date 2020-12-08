/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file geodesic & clustering on mesh using dijkstra method
 */

#ifndef DFM2_EXPMAP_H
#define DFM2_EXPMAP_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/dijkstra.h"
#include <queue>
#include <climits>
#include <random>

namespace delfem2 {

class CExpMap_DijkstraPoint {
public:
  const std::vector<double> &aXYZ;
  const std::vector<double> &aNorm;
  const std::vector<unsigned int> &aTri;
  std::vector<double> &aTex;
  const std::vector<unsigned int> &psup_ind;
  const std::vector<unsigned int> &psup;
  std::vector<double> aW;
  std::vector<CVec3d> aAxisX;
public:
  CExpMap_DijkstraPoint(
      unsigned int ip_ker,
      const std::vector<double> &aXYZ_,
      const std::vector<double> &aNorm_,
      const std::vector<unsigned int> &aTri_,
      std::vector<double> &aTex_,
      std::vector<unsigned int> &psup_ind_,
      std::vector<unsigned int> &psup_) :
      aXYZ(aXYZ_), aNorm(aNorm_), aTri(aTri_), aTex(aTex_),
      psup_ind(psup_ind_), psup(psup_) {
    const unsigned int np = aXYZ.size() / 3;
    aAxisX.resize(np, CVec3d(0, 0, 0));
    aTex.resize(np * 2);
    aW.assign(np, 0.0);

    { // set kernel point information
      aTex[ip_ker * 2 + 0] = 0.0;
      aTex[ip_ker * 2 + 1] = 0.0;
      CVec3d y0;
      GetVertical2Vector(
          CVec3d(aNorm.data() + ip_ker * 3),
          aAxisX[ip_ker], y0);
      aW[ip_ker] = 1.0;
    }
  }

  /**
   *
   * @param[in] ip0 point newly added
   * @param[in] io0 the order of newly added point. ( 0 if this is the first point)
   * @param[in] aOrder map from point index 2 the order of points added
   */
  void AddPoint(
      unsigned int ip0,
      std::vector<unsigned int> &aOrder) {
    assert(aOrder.size() == aXYZ.size() / 3);
    const CVec3d n0 = CVec3d(aNorm.data() + ip0 * 3).Normalize();
    aAxisX[ip0].SetNormalizedVector();
    const CVec3d x0 = aAxisX[ip0];
    const CVec3d y0 = n0 ^x0;
    aTex[ip0 * 2 + 0] /= aW[ip0];
    aTex[ip0 * 2 + 1] /= aW[ip0];
    for (unsigned int ipsup = psup_ind[ip0]; ipsup < psup_ind[ip0 + 1]; ++ipsup) {
      const unsigned int ip1 = psup[ipsup];
      if (aOrder[ip1] != UINT_MAX) { continue; } // effect propagate from fixed to unfixed
      const CVec3d n1 = CVec3d(aNorm.data() + ip1 * 3).Normalize();
      const CVec3d d01 = CVec3d(aXYZ.data() + ip1 * 3) - CVec3d(aXYZ.data() + ip0 * 3);
      const double len01 = d01.Length();
      const CVec3d e01 = (d01 - (d01 * n0) * n0).Normalize() * len01; // projected edge and same length
      const double w01 = 1.0 / len01;
      const CVec3d x1 = Mat3_MinimumRotation(n0, n1) * x0;
      aAxisX[ip1] += w01 * x1;
      aW[ip1] += w01;
      aTex[ip1 * 2 + 0] += w01 * (aTex[ip0 * 2 + 0] + (e01 * x0));
      aTex[ip1 * 2 + 1] += w01 * (aTex[ip0 * 2 + 1] + (e01 * y0));
    }
  }
};


class CExpMap_DijkstraElem{
public:
  std::vector<double>& aTex;
  const std::vector<double>& aXYZ;
  const std::vector<unsigned int>& aTri;
  const std::vector<unsigned int>& aTriSuTri;
  std::vector<double> aW;
  std::vector<CVec3d> aAxisX;
public:
  CExpMap_DijkstraElem(
      std::vector<double>& aTex_,
      unsigned int it_ker,
      const std::vector<double>& aXYZ_,
      const std::vector<unsigned int>& aTri_,
      std::vector<unsigned int>& aTriSuTri_) :
      aTex(aTex_), aXYZ(aXYZ_), aTri(aTri_), aTriSuTri(aTriSuTri_)
  {
    const unsigned int ntri = aTri.size() / 3;
    aAxisX.resize(ntri, CVec3d(0, 0, 0));
    aTex.resize(ntri * 2);
    aW.assign(ntri, 0.0);
    { // set kernel point information
      aTex[it_ker * 2 + 0] = 0.0;
      aTex[it_ker * 2 + 1] = 0.0;
      CVec3d n0 = Normal_Tri3(it_ker,aTri,aXYZ).Normalize();
      CVec3d y0;
      GetVertical2Vector(
          n0,
          aAxisX[it_ker], y0);
      aW[it_ker] = 1.0;
    }
  }

  /**
   *
   * @param[in] ip0 point newly added
   * @param[in] io0 the order of newly added point. ( 0 if this is the first point)
   * @param[in] aOrder map from point index 2 the order of points added
   */
  void AddElem(
      unsigned int it0,
      std::vector<unsigned int>& aOrder)
  {
    assert( aOrder.size() == aTri.size()/3 );
    assert( aOrder[it0] != UINT_MAX );
    const CVec3d n0 = Normal_Tri3(it0,aTri,aXYZ).Normalize();
    const CVec3d p0 = CG_Tri3(it0,aTri,aXYZ);
    aAxisX[it0].SetNormalizedVector();
    const CVec3d x0 = aAxisX[it0];
    const CVec3d y0 = n0^x0;
    aTex[it0 * 2 + 0] /= aW[it0];
    aTex[it0 * 2 + 1] /= aW[it0];
    for (unsigned int iedge = 0; iedge < 3; ++iedge) {
      const unsigned int it1 = aTriSuTri[it0*3+iedge];
      if( it1 == UINT_MAX ){ continue; }
      if ( aOrder[it1] != UINT_MAX ) { continue; } // effect propagate from fixed to unfixed
      const CVec3d n1 = Normal_Tri3(it1,aTri,aXYZ).Normalize();
      const CVec3d d01 = CG_Tri3(it1,aTri,aXYZ)-p0;
      const double len01 = d01.Length();
      const CVec3d e01 = (d01 - (d01 * n0) * n0).Normalize() * len01; // projected edge and same length
      const double w01 = 1.0 / len01;
      const CVec3d x1 = Mat3_MinimumRotation(n0, n1) * x0;
      aAxisX[it1] += w01 * x1;
      aW[it1] += w01;
      aTex[it1 * 2 + 0] += w01 * (aTex[it0 * 2 + 0] + (e01 * x0));
      aTex[it1 * 2 + 1] += w01 * (aTex[it0 * 2 + 1] + (e01 * y0));
    }
  }
};

}

#endif