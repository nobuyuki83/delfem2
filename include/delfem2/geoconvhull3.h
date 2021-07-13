/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @file convex hull computation in 3D space
 */

#ifndef DFM2_GEOCONVHULL3_H
#define DFM2_GEOCONVHULL3_H

#include "delfem2/dfm2_inline.h"
#include <cassert>
#include <memory>
#include <cmath>
#include <iostream>
#include <vector>
#include <climits>
#include <stack>

namespace delfem2 {

/**
 *
 * @tparam REAL float,double ...etc
 * @tparam VEC Eigen::Vector3f, DelFEM2::CVector3f ...etc
 * @tparam ALLOCATOR (typycally guessed by the compiler)
 * @param[out] aTri triangle index
 * @param[in] aXYZ std::vector array of vector type
 */
template<typename REAL, typename VEC, typename ALLOCATOR>
void ConvexHull(
    std::vector<unsigned int> &aTri,
    const std::vector<VEC, ALLOCATOR> &aXYZ);

} // end namespace delfem2

// ----------------------------------
// local functions

namespace delfem2 {
namespace convhull3 {

template<class VEC>
inline void Normal(
    VEC &vnorm,
    const VEC &v1,
    const VEC &v2,
    const VEC &v3)
{
  vnorm(0) = (v2(1) - v1(1)) * (v3(2) - v1(2)) - (v2(2) - v1(2)) * (v3(1) - v1(1));
  vnorm(1) = (v2(2) - v1(2)) * (v3(0) - v1(0)) - (v2(0) - v1(0)) * (v3(2) - v1(2));
  vnorm(2) = (v2(0) - v1(0)) * (v3(1) - v1(1)) - (v2(1) - v1(1)) * (v3(0) - v1(0));
}

template<typename REAL, class VEC>
inline REAL Dot3(
    const VEC &v1,
    const VEC &v2)
{
  return v1(0)*v2(0) + v1(1)*v2(1) + v1(2)*v2(2);
}

template<typename REAL, typename VEC, typename ALLOCATOR>
bool IsOut(
    unsigned int itri,
    const VEC &v,
    const std::vector<VEC, ALLOCATOR> &aXYZ,
    const std::vector<unsigned int> &aTri)
{
  const unsigned int i0 = aTri[itri * 3 + 0];
  const unsigned int i1 = aTri[itri * 3 + 1];
  const unsigned int i2 = aTri[itri * 3 + 2];
  const VEC &v0 = aXYZ[i0];
  const VEC &v1 = aXYZ[i1];
  const VEC &v2 = aXYZ[i2];
  VEC n; n.setZero();
  Normal(n, v0, v1, v2);
//  REAL dot = Dot3<REAL,VEC>(v-v0,n);
  REAL dot = n.dot(v-v0);
  return dot > 0;
}

//! Volume of a tetrahedra
template <typename REAL, typename VEC>
REAL Volume_Tet(
    const VEC &v0,
    const VEC &v1,
    const VEC &v2,
    const VEC &v3)
{
  REAL v
      = (v1(0) - v0(0)) * ((v2(1) - v0(1)) * (v3(2) - v0(2)) - (v3(1) - v0(1)) * (v2(2) - v0(2)))
      + (v1(1) - v0(1)) * ((v2(2) - v0(2)) * (v3(0) - v0(0)) - (v3(2) - v0(2)) * (v2(0) - v0(0)))
      + (v1(2) - v0(2)) * ((v2(0) - v0(0)) * (v3(1) - v0(1)) - (v3(0) - v0(0)) * (v2(1) - v0(1)));
  return v * 0.16666666666666666666666666666667;
}

}
}


template<typename REAL, typename VEC, typename ALLOCATOR>
void delfem2::ConvexHull(
    std::vector<unsigned int> &aTri,
    const std::vector<VEC, ALLOCATOR> &aXYZ)
{
  namespace lcl = ::delfem2::convhull3;
  std::vector<int> aBflg(aXYZ.size(), -1);
  aTri.reserve(aXYZ.size() * 6);
  aTri = {
      1,2,3,
      0,3,2,
      0,1,3,
      0,2,1
  };
  std::vector<std::pair<unsigned int, int> > aTriSur;
  aTriSur = {
      {1,0},{2,0},{3,0},{0,0},{3,2}, {2,1},
      {0,1},{1,2},{3,1},{0,2},{2,2},{1,1}
  };
  {
    REAL vol = lcl::Volume_Tet<REAL>(aXYZ[0], aXYZ[1], aXYZ[2], aXYZ[3]);
    if (vol < 0) {
      aTri = {
          3,2,1,
          2,3,0,
          3,1,0,
          1,2,0 };
      aTriSur = {
          {3,2,},{2,2},{1,2},{2,1},{3,0},{0,2},
          {3,1},{1,0},{0,1},{1,1},{2,0},{0,0} };
    }
  }
  const int triEd[3][2] = {
      {1, 2},
      {2, 0},
      {0, 1}};
  for (unsigned iv = 4; iv < static_cast<unsigned int>(aXYZ.size()); iv++) {
    VEC v = aXYZ[iv];
    int itri_ker = -1;
	const unsigned int ntri = static_cast<unsigned int>(aTri.size() / 3);
    for (unsigned int itri = 0; itri < ntri; itri++) {
      if (lcl::IsOut<REAL>(itri, v, aXYZ, aTri)) {
        itri_ker = itri;
        break;
      }
    }
#ifndef NDEBUG
    {
      for (std::size_t itri = 0; itri < aTri.size() / 3; itri++) {
        for (int ied = 0; ied < 3; ied++) {
          int ied1 = triEd[ied][0];
          int ied2 = triEd[ied][1];
          int itri_s = aTriSur[itri * 3 + ied].first;
          int ied_s0 = aTriSur[itri * 3 + ied].second;
          assert(aTriSur[itri_s * 3 + ied_s0].first == itri);
          assert(aTriSur[itri_s * 3 + ied_s0].second == ied);
          int ied_s1 = triEd[ied_s0][0];
          int ied_s2 = triEd[ied_s0][1];
          assert(aTri[itri * 3 + ied1] == aTri[itri_s * 3 + ied_s2]);
          assert(aTri[itri * 3 + ied2] == aTri[itri_s * 3 + ied_s1]);
        }
      }
    }
#endif
    if (itri_ker == -1) continue; // inside
    std::vector<std::pair<int, int> > aB;
    std::vector<int> isDelTri(aTri.size() / 3, -1);
    {
      std::vector<int> isLookedEdge(aTri.size(), -1);
      std::stack<std::pair<unsigned int, int> > sBound; // itri,iedge
      { // initialize
        sBound.push(aTriSur[itri_ker * 3 + 0]);
        sBound.push(aTriSur[itri_ker * 3 + 1]);
        sBound.push(aTriSur[itri_ker * 3 + 2]);
        isDelTri[itri_ker] = 1;
      }
      for (;;) {
        if (sBound.empty()) break;
        const unsigned int itri0 = sBound.top().first;
        const int ied0 = sBound.top().second;
        sBound.pop();
        if (isLookedEdge[itri0 * 3 + ied0] == 1) continue;
        isLookedEdge[itri0 * 3 + ied0] = 1;
        {
          const std::pair<int, int> &s0 = aTriSur[itri0 * 3 + ied0];
          isLookedEdge[s0.first * 3 + s0.second] = 1;
        }
        isDelTri[itri0] = (::delfem2::convhull3::IsOut<REAL>(itri0, v, aXYZ, aTri)) ? 1 : 0;
        if (isDelTri[itri0] == 1) { // expand this boundary
          int ied1 = triEd[ied0][0];
          int ied2 = triEd[ied0][1];
          const std::pair<int, int> &s1 = aTriSur[itri0 * 3 + ied1];
          const std::pair<int, int> &s2 = aTriSur[itri0 * 3 + ied2];
          assert(aTriSur[s1.first * 3 + s1.second].first == itri0);
          assert(aTriSur[s2.first * 3 + s2.second].first == itri0);
          sBound.push(s1);
          sBound.push(s2);
        } else { // this is a actuall boundary
          aB.emplace_back(itri0, ied0);
        }
      }
    }
    std::vector<unsigned int> aBSur(aB.size() * 2, UINT_MAX);
    {
      for (auto &ib : aB) {
        int itri0 = ib.first;
        int itn0 = ib.second;
        int itn1 = triEd[itn0][0];
        int itn2 = triEd[itn0][1];
        int iv1 = aTri[itri0 * 3 + itn1];
        int iv2 = aTri[itri0 * 3 + itn2];
        aBflg[iv1] = -1;
        aBflg[iv2] = -1;
      }
      for (std::size_t ib = 0; ib < aB.size(); ib++) {
        int itri0 = aB[ib].first;
        int itn0 = aB[ib].second;
        int itn1 = triEd[itn0][0];
        int itn2 = triEd[itn0][1];
        int iv1 = aTri[itri0 * 3 + itn1];
        int iv2 = aTri[itri0 * 3 + itn2];
        if (aBflg[iv1] == -2) {}
        else if (aBflg[iv1] == -1) { aBflg[iv1] = static_cast<int>(ib); }
        else {
          assert(aBflg[iv1] >= 0);
          int ib0 = aBflg[iv1];
          aBSur[ib * 2 + 1] = ib0;
          aBSur[ib0 * 2 + 0] = static_cast<int>(ib);
          aBflg[iv1] = -2;
        }
        if (aBflg[iv2] == -2) {}
        else if (aBflg[iv2] == -1) { aBflg[iv2] = static_cast<int>(ib); }
        else {
          assert(aBflg[iv2] >= 0);
          int ib0 = aBflg[iv2];
          aBSur[ib * 2 + 0] = ib0;
          aBSur[ib0 * 2 + 1] = static_cast<int>(ib);
          aBflg[iv2] = -2;
        }
      }
    }
#ifndef NDEBUG
    for (std::size_t ib = 0; ib < aB.size(); ib++) {
      for (int inb = 0; inb < 2; inb++) {
        int itri0 = aB[ib].first;
        int itn0 = aB[ib].second;
        int iv1 = aTri[itri0 * 3 + triEd[itn0][0]];
        int iv2 = aTri[itri0 * 3 + triEd[itn0][1]];
        int ib_s0 = aBSur[ib * 2 + inb];
        int itri_s0 = aB[ib_s0].first;
        int itn_s0 = aB[ib_s0].second;
        int iv_s1 = aTri[itri_s0 * 3 + triEd[itn_s0][0]];
        int iv_s2 = aTri[itri_s0 * 3 + triEd[itn_s0][1]];
        if (inb == 0) { assert(iv2 == iv_s1); }
        else { assert(iv1 == iv_s2); }
      }
    }
#endif
    std::vector<unsigned int> mapOld2New(aTri.size() / 3, UINT_MAX);
    std::vector<unsigned int> aTri1;
    std::vector<std::pair<unsigned int, int> > aTriSur1;
    aTri1.reserve(aTri.size() + 60);
    aTriSur1.reserve(aTriSur.size() + 60);
    for (unsigned int itri = 0; itri < ntri; itri++) {
      if (isDelTri[itri] == 1) continue;
      assert(!::delfem2::convhull3::IsOut<REAL>(itri, v, aXYZ, aTri));
      // itri is inside
      mapOld2New[itri] = static_cast<unsigned int>(aTri1.size() / 3);
      aTri1.push_back(aTri[itri * 3 + 0]);
      aTri1.push_back(aTri[itri * 3 + 1]);
      aTri1.push_back(aTri[itri * 3 + 2]);
      aTriSur1.emplace_back(UINT_MAX, 0);
      aTriSur1.emplace_back(UINT_MAX, 0);
      aTriSur1.emplace_back(UINT_MAX, 0);
    }
    for (unsigned int itri = 0; itri < aTri.size() / 3; itri++) { // set old relation
      if (isDelTri[itri] == 1) continue;
      unsigned int jtri0 = mapOld2New[itri];
      assert(jtri0 < aTri1.size() / 3);
      for (int iet = 0; iet < 3; iet++) {
        int itri_s = aTriSur[itri * 3 + iet].first;
        if (mapOld2New[itri_s] == UINT_MAX) continue;
        aTriSur1[jtri0 * 3 + iet].first = mapOld2New[itri_s];
        aTriSur1[jtri0 * 3 + iet].second = aTriSur[itri * 3 + iet].second;
      }
    }
    const unsigned int ntri_old = static_cast<unsigned int>(aTri1.size() / 3);
    for (std::size_t ib = 0; ib < aB.size(); ib++) {
      const unsigned int itri0 = aB[ib].first;
      int itn0 = aB[ib].second;
      int itn1 = triEd[itn0][0];
      int itn2 = triEd[itn0][1];
      assert(!lcl::IsOut<REAL>(itri0, v, aXYZ, aTri));
#ifndef NDEBUG
      {
        const unsigned int itri_s = aTriSur[itri0 * 3 + itn0].first;
        assert(lcl::IsOut<REAL>(itri_s, v, aXYZ, aTri));
        const int ied_s0 = aTriSur[itri0 * 3 + itn0].second;
        assert(aTriSur[itri_s * 3 + ied_s0].first == itri0);
        assert(aTriSur[itri_s * 3 + ied_s0].second == itn0);
        int ied_s1 = triEd[ied_s0][0];
        int ied_s2 = triEd[ied_s0][1];
        assert(aTri[itri0 * 3 + itn1] == aTri[itri_s * 3 + ied_s2]);
        assert(aTri[itri0 * 3 + itn2] == aTri[itri_s * 3 + ied_s1]);
      }
#endif
      assert(isDelTri[itri0] == 0);
      const unsigned int jtri0 = mapOld2New[itri0];
      assert(jtri0 != UINT_MAX);
      const unsigned int jtri1 = static_cast<unsigned int>(aTri1.size() / 3);
      assert(jtri1 == ntri_old + ib);
      aTri1.push_back(iv);
      aTri1.push_back(aTri[itri0 * 3 + itn2]);
      aTri1.push_back(aTri[itri0 * 3 + itn1]);
      aTriSur1[jtri0 * 3 + itn0] = std::make_pair(jtri1, 0);
      //
      const unsigned int jtri2 = aBSur[ib * 2 + 0] + ntri_old;
      const unsigned int jtri3 = aBSur[ib * 2 + 1] + ntri_old;
      aTriSur1.emplace_back(jtri0, itn0);
      aTriSur1.emplace_back(jtri3, 2);
      aTriSur1.emplace_back(jtri2, 1);
    }
    aTri = aTri1;
    aTriSur = aTriSur1;
  }
}


#endif // VEC3_H
