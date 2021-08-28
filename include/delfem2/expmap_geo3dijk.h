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

#include <queue>
#include <climits>
#include <random>

#include "delfem2/dfm2_inline.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/dijkstra.h"

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
    const size_t np = aXYZ.size() / 3;
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
    const CVec3d n0 = CVec3d(aNorm.data() + ip0 * 3).normalized();
    aAxisX[ip0].normalize();
    const CVec3d x0 = aAxisX[ip0];
    const CVec3d y0 = n0 ^x0;
    aTex[ip0 * 2 + 0] /= aW[ip0];
    aTex[ip0 * 2 + 1] /= aW[ip0];
    for (unsigned int ipsup = psup_ind[ip0]; ipsup < psup_ind[ip0 + 1]; ++ipsup) {
      const unsigned int ip1 = psup[ipsup];
      if (aOrder[ip1] != UINT_MAX) { continue; } // effect propagate from fixed to unfixed
      const CVec3d n1 = CVec3d(aNorm.data() + ip1 * 3).normalized();
      const CVec3d d01 = CVec3d(aXYZ.data() + ip1 * 3) - CVec3d(aXYZ.data() + ip0 * 3);
      const double len01 = d01.norm();
      const CVec3d e01 = (d01 - (d01.dot(n0)) * n0).normalized() * len01; // projected edge and same length
      const double w01 = 1.0 / len01;
      const CVec3d x1 = Mat3_MinimumRotation(n0, n1) * x0;
      aAxisX[ip1] += w01 * x1;
      aW[ip1] += w01;
      aTex[ip1 * 2 + 0] += w01 * (aTex[ip0 * 2 + 0] + (e01.dot(x0)));
      aTex[ip1 * 2 + 1] += w01 * (aTex[ip0 * 2 + 1] + (e01.dot(y0)));
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
      const std::vector<unsigned int>& aTriSuTri_) :
      aTex(aTex_), aXYZ(aXYZ_), aTri(aTri_), aTriSuTri(aTriSuTri_)
  {
    const size_t ntri = aTri.size() / 3;
    assert(it_ker<ntri);
    aAxisX.resize(ntri, CVec3d(0, 0, 0));
    aTex.resize(ntri * 2);
    aW.assign(ntri, 0.0);
    { // set kernel point information
      aTex[it_ker * 2 + 0] = 0.0;
      aTex[it_ker * 2 + 1] = 0.0;
      CVec3d n0 = Normal_Tri3(it_ker,aTri,aXYZ).normalized();
      CVec3d y0;
      GetVertical2Vector(
          n0,
          aAxisX[it_ker], y0);
      aW[it_ker] = 1.0;
    }
  }

  virtual bool IsIncludeElem(
	  [[maybe_unused]] unsigned int ie){ return true; }

  /**
   *
   * @param[in] ip0 point newly added
   * @param[in] io0 the order of newly added point. ( 0 if this is the first point)
   * @param[in] aOrder map from point index 2 the order of points added
   */
  virtual void AddElem(
      unsigned int it0,
      std::vector<unsigned int>& aOrder)
  {
    assert( aOrder.size() == aTri.size()/3 );
    assert( aOrder[it0] != UINT_MAX );
    const CVec3d n0 = Normal_Tri3(it0,aTri,aXYZ).normalized();
    const CVec3d p0 = CG_Tri3(it0,aTri,aXYZ);
    aAxisX[it0].normalize();
    const CVec3d x0 = aAxisX[it0];
    const CVec3d y0 = n0^x0;
    aTex[it0 * 2 + 0] /= aW[it0];
    aTex[it0 * 2 + 1] /= aW[it0];
    for (unsigned int iedge = 0; iedge < 3; ++iedge) {
      const unsigned int it1 = aTriSuTri[it0*3+iedge];
      if( it1 == UINT_MAX ){ continue; }
      if ( aOrder[it1] != UINT_MAX ) { continue; } // effect propagate from fixed to unfixed
      const CVec3d n1 = Normal_Tri3(it1,aTri,aXYZ).normalized();
      const CVec3d d01 = CG_Tri3(it1,aTri,aXYZ)-p0;
      const double len01 = d01.norm();
      const CVec3d e01 = (d01 - (d01.dot(n0)) * n0).normalized() * len01; // projected edge and same length
      const double w01 = 1.0 / len01;
      const CVec3d x1 = Mat3_MinimumRotation(n0, n1) * x0;
      aAxisX[it1] += w01 * x1;
      aW[it1] += w01;
      aTex[it1 * 2 + 0] += w01 * (aTex[it0 * 2 + 0] + (e01.dot(x0)));
      aTex[it1 * 2 + 1] += w01 * (aTex[it0 * 2 + 1] + (e01.dot(y0)));
    }
  }
};


class CExpMap_DijkstraElemFlag
    : public CExpMap_DijkstraElem
{
public:
  const std::vector<unsigned int>& aFlgElem;
  unsigned int iflag0;
public:
  CExpMap_DijkstraElemFlag(
      std::vector<double>& aTex_,
      unsigned int it_ker,
      const std::vector<double>& aXYZ_,
      const std::vector<unsigned int>& aTri_,
      const std::vector<unsigned int>& aTriSuTri_,
      const std::vector<unsigned int>& aFlgElem_) :
      CExpMap_DijkstraElem(aTex_,it_ker,aXYZ_,aTri_,aTriSuTri_), aFlgElem(aFlgElem_)
  {
    iflag0 = aFlgElem[it_ker];
  }

  void AddElem(
      unsigned int it0,
      std::vector<unsigned int>& aOrder) override{
    assert( aFlgElem[it0] == iflag0 );
    CExpMap_DijkstraElem::AddElem(it0,aOrder);
  }

  bool IsIncludeElem(unsigned int ie) override{
    return aFlgElem[ie] == iflag0;
  }
};


void TexPoint_TexElemFlag(
    std::vector<double> &aTexP,
    const std::vector<double> &aXYZ,
    const std::vector<double> &aTexE,
    const std::vector<unsigned int> &aTri,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    unsigned int iflg,
    std::vector<unsigned int> aFlgTri)
{
  aTexP.resize(aXYZ.size() / 3 * 2);
  for (unsigned int ip0 = 0; ip0 < aXYZ.size() / 3; ++ip0) {
    const CVec3d p0(aXYZ.data() + ip0 * 3);
    double tex[2] = {0, 0};
    double w0 = 0.0;
    for (unsigned int ielsup = elsup_ind[ip0]; ielsup < elsup_ind[ip0 + 1]; ++ielsup) {
      const unsigned int it1 = elsup[ielsup];
      if( aFlgTri[it1] != iflg ){ continue; }
      CVec3d p1 = CG_Tri3(it1, aTri, aXYZ);
      double w1 = 1.0 / (p0 - p1).norm();
      tex[0] += w1 * aTexE[it1 * 2 + 0];
      tex[1] += w1 * aTexE[it1 * 2 + 1];
      w0 += w1;
    }
    tex[0] /= w0;
    tex[1] /= w0;
    aTexP[ip0 * 2 + 0] = tex[0];
    aTexP[ip0 * 2 + 1] = tex[1];
  }
}

void FlatteringPattern(
    unsigned int &ielm_ker,
    CVec3d coordLocal[4],
    std::vector<double> &aTexP,
    //
    const std::vector<double> &aXYZ,
    const std::vector<unsigned int> &aTri,
    const std::vector<unsigned int> &aTriSuTri,
    unsigned int iflg0,
    const std::vector<unsigned int> &aFlgTri,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup) {
  {
    class CProc {
    public:
      unsigned int ie0 = UINT_MAX;
    public:
      void AddElem(unsigned int iel, 
		  [[maybe_unused]] std::vector<unsigned> &aOrder) { ie0 = iel; }
    } proc;
    std::vector<unsigned int> aOrder, aDist;
    Dijkstra_FillFromBoundary(
        aOrder, aDist,
        proc, iflg0, aFlgTri, aTriSuTri);
    ielm_ker = proc.ie0;
    assert(ielm_ker<aTri.size()/3);
  }
  {
    std::vector<double> aTexE; // element-wise texture coordinate
    CExpMap_DijkstraElemFlag expmap(
        aTexE,
        ielm_ker, aXYZ, aTri, aTriSuTri,
        aFlgTri);
    std::vector<double> aDist;
    std::vector<unsigned int> aOrder;
    DijkstraElem_MeshElemGeo3(
        aDist, aOrder, expmap,
        ielm_ker,
        aTri, 
		static_cast<unsigned int>(aTri.size() / 3),
        aXYZ,
        aTriSuTri);
    coordLocal[0] = expmap.aAxisX[ielm_ker];
    coordLocal[2] = Normal_Tri3(ielm_ker, aTri, aXYZ).normalized();
    assert(fabs(coordLocal[0].dot(coordLocal[2])) < 1.0e-10);
    assert(fabs(coordLocal[0].norm() - 1.0) < 1.0e-10);
    assert(fabs(coordLocal[2].norm() - 1.0) < 1.0e-10);
    coordLocal[1] = coordLocal[2] ^ coordLocal[0];
    coordLocal[3] = CG_Tri3(ielm_ker, aTri, aXYZ);
    //
    TexPoint_TexElemFlag(aTexP,
        aXYZ, aTexE, aTri,
        elsup_ind, elsup,
        iflg0, aFlgTri);
  }
}

}

#endif
