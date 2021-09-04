/*
 * Copyright (c) 2019-2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/dtri2_v2dtri.h"

#include <set>
#include <algorithm>
#include <climits>

#include "delfem2/geo_polygon2.h"

// ===========================================
// unexposed functions

namespace delfem2 {
namespace dtri2 {

DFM2_INLINE bool LaplacianArroundPoint
    (std::vector<CVec2d> &aVec2,
     const unsigned int ipoin,
     const std::vector<CDynPntSur> &aPo,
     const std::vector<CDynTri> &aTri) {
  assert(aVec2.size() == aPo.size());
  unsigned int itri0 = aPo[ipoin].e;
  unsigned int ino0 = aPo[ipoin].d;
  CVec2d vec_delta = aVec2[ipoin];
  unsigned int ntri_around = 1;
  for (;;) { // counter-clock wise
    assert(itri0 < aTri.size() && ino0 < 3 && aTri[itri0].v[ino0] == ipoin);
    vec_delta += aVec2[aTri[itri0].v[(ino0 + 1) % 3]];
    ntri_around++;
    if (!MoveCCW(itri0, ino0, UINT_MAX, aTri)) { return false; }
    if (itri0 == aPo[ipoin].e) { break; }
  }
  aVec2[ipoin] = vec_delta / ntri_around;
  return true;
}

DFM2_INLINE bool FindEdgePoint_AcrossEdge(
    unsigned int &itri0,
    unsigned int &inotri0,
    unsigned &inotri1,
    double &ratio,
    unsigned int ipo0,
    unsigned int ipo1,
    const std::vector<CDynPntSur> &aDPo,
    const std::vector<CDynTri> &aDTri,
    const std::vector<CVec2d> &aVec2) {
  const unsigned int itri_ini = aDPo[ipo0].e;
  const unsigned int inotri_ini = aDPo[ipo0].d;
  unsigned int inotri_cur = inotri_ini;
  unsigned int itri_cur = itri_ini;
  for (;;) {
    assert(aDTri[itri_cur].v[inotri_cur] == ipo0);
    {
      const unsigned int inotri2 = (inotri_cur + 1) % 3;
      const unsigned int inotri3 = (inotri_cur + 2) % 3;
      double area0 = Area_Tri(aVec2[ipo0],
                              aVec2[aDTri[itri_cur].v[inotri2]],
                              aVec2[ipo1]);
      if (area0 > -1.0e-20) {
        double area1 = Area_Tri(aVec2[ipo0],
                                aVec2[ipo1],
                                aVec2[aDTri[itri_cur].v[inotri3]]);
        if (area1 > -1.0e-20) {
          assert(area0 + area1 > 1.0e-20);
          ratio = area0 / (area0 + area1);
          itri0 = itri_cur;
          inotri0 = inotri2;
          inotri1 = inotri3;
          return true;
        }
      }
    }
    {
      const unsigned int inotri2 = (inotri_cur + 1) % 3;
      const unsigned int itri_nex = aDTri[itri_cur].s2[inotri2];
      if (itri_nex == UINT_MAX) { break; }
      const unsigned int jnob = FindAdjEdgeIndex(aDTri[itri_nex], inotri2, aDTri);
      const unsigned int inotri3 = (jnob + 1) % 3;
      assert(itri_nex < aDTri.size());
      assert(aDTri[itri_nex].v[inotri3] == ipo0);
      if (itri_nex == itri_ini) {
        itri0 = 0;
        inotri0 = 0;
        inotri1 = 0;
        ratio = 0.0;
        return false;
      }
      itri_cur = itri_nex;
      inotri_cur = inotri3;
    }
  }

  inotri_cur = inotri_ini;
  itri_cur = itri_ini;
  for (;;) {
    assert(aDTri[itri_cur].v[inotri_cur] == ipo0);
    {
      const unsigned int inotri2 = (inotri_cur + 1) % 3; // indexRot3[1][inotri_cur];
      const unsigned int inotri3 = (inotri_cur + 2) % 3; // indexRot3[2][inotri_cur];
      double area0 = Area_Tri(aVec2[ipo0],
                              aVec2[aDTri[itri_cur].v[inotri2]],
                              aVec2[ipo1]);
      if (area0 > -1.0e-20) {
        double area1 = Area_Tri(aVec2[ipo0],
                                aVec2[ipo1],
                                aVec2[aDTri[itri_cur].v[inotri3]]);
        if (area1 > -1.0e-20) {
          assert(area0 + area1 > 1.0e-20);
          ratio = area0 / (area0 + area1);
          itri0 = itri_cur;
          inotri0 = inotri2;
          inotri1 = inotri3;
          return true;
        }
      }
    }
    {
      const unsigned int inotri2 = (inotri_cur + 2) % 3;
      unsigned int itri_nex = aDTri[itri_cur].s2[inotri2];
      unsigned int jnob = FindAdjEdgeIndex(aDTri[itri_cur], inotri2, aDTri);
      const unsigned int inotri3 = (jnob + 1) % 3;
      assert(aDTri[itri_nex].v[inotri3] == ipo0);
      if (itri_nex == itri_ini) {
        assert(0);
      }
      itri_cur = itri_nex;
      inotri_cur = inotri3;
    }
  }

  itri0 = 0;
  inotri0 = 0;
  inotri1 = 0;
  ratio = 0.0;

  return false;
}

} // cad2
} // delfem2

// unexposed functions
// --------------------------------------------------------
// exposed functions

DFM2_INLINE void delfem2::CheckTri(
    [[maybe_unused]] const std::vector<CDynPntSur> &aPo3D,
    [[maybe_unused]] const std::vector<CDynTri> &aDTri,
    [[maybe_unused]] const std::vector<CVec2d> &aXYZ) {
#if !defined(NDEBUG)
  for (const auto &tri : aDTri) {
    const unsigned int i0 = tri.v[0];
    assert(i0 < aPo3D.size());
    if (i0 == UINT_MAX) continue;
    const unsigned int i1 = tri.v[1];
    assert(i1 < aPo3D.size());
    const unsigned int i2 = tri.v[2];
    assert(i2 < aPo3D.size());
    const double area = Area_Tri(aXYZ[i0], aXYZ[i1], aXYZ[i2]);
    if (area < 1.0e-10) {  // very small volume
      assert(0);
      abort();
    }
  }
#endif
}

DFM2_INLINE void delfem2::DelaunayAroundPoint(
    const unsigned int ipo0,
    std::vector<CDynPntSur> &aPo,
    std::vector<CDynTri> &aTri,
    const std::vector<CVec2d> &aVec2) {
  assert(aVec2.size() == aPo.size());
  assert(ipo0 < aPo.size());
  if (aPo[ipo0].e == UINT_MAX) { return; }

  unsigned int itri0 = aPo[ipo0].e;
  unsigned int ino0 = aPo[ipo0].d;

  // ---------------------------
  // go counter-clock-wise
  bool flag_is_wall = false;
  for (;;) {
    assert(itri0 < aTri.size() && ino0 < 3 && aTri[itri0].v[ino0] == ipo0);
    if (aTri[itri0].s2[ino0] < aTri.size()) {
      const unsigned int jtri0 = aTri[itri0].s2[ino0];
      const unsigned int jno0 = FindAdjEdgeIndex(aTri[itri0], ino0, aTri);
      assert(aTri[jtri0].s2[jno0] == itri0);
      const unsigned int jpo0 = aTri[jtri0].v[jno0];
      const int ires = DetDelaunay(
          aVec2[aTri[itri0].v[0]],
          aVec2[aTri[itri0].v[1]],
          aVec2[aTri[itri0].v[2]],
          aVec2[jpo0]);
      if (ires == 0) {
        FlipEdge(itri0, ino0, aPo, aTri); // this edge is not on the edge and should be successfull
        ino0 = 2;
        assert(aTri[itri0].v[ino0] == ipo0); // this is the rule from FlipEdge function
        continue; // need to check the fliped element
      }
    }
    if (!MoveCCW(itri0, ino0, UINT_MAX, aTri)) {
      flag_is_wall = true;
      break;
    }
    if (itri0 == aPo[ipo0].e) break;
  }
  if (!flag_is_wall) { return; }

  // ----------------------------
  // go clock-wise

  for (;;) {
    assert(itri0 < aTri.size() && ino0 < 3 && aTri[itri0].v[ino0] == ipo0);
    if (aTri[itri0].s2[ino0] < aTri.size()) {
      const unsigned int jtri0 = aTri[itri0].s2[ino0];
      const unsigned int jno0 = FindAdjEdgeIndex(aTri[itri0], ino0, aTri);
      assert(aTri[jtri0].s2[jno0] == itri0);
      const unsigned int ipo_dia = aTri[jtri0].v[jno0];
      const int ires = DetDelaunay(
          aVec2[aTri[itri0].v[0]],
          aVec2[aTri[itri0].v[1]],
          aVec2[aTri[itri0].v[2]],
          aVec2[ipo_dia]);
      if (ires == 0) { // Delaunay condition is not satisfiled
        FlipEdge(itri0, ino0, aPo, aTri);
        itri0 = jtri0;
        ino0 = 1;
        assert(aTri[itri0].v[ino0] == ipo0);
        continue;
      }
    }
    if (!MoveCW(itri0, ino0, UINT_MAX, aTri)) { return; }
  }
}

DFM2_INLINE void delfem2::MeshingInside(
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CDynTri> &aTri,
    std::vector<CVec2d> &aVec2,
    std::vector<int> &aFlagPnt,
    std::vector<unsigned int> &aFlagTri,
    const size_t nPointFix,
    const unsigned int nflgpnt_offset,
    const double len,
    const CInputTriangulation &mesh_density) {
  assert(aVec2.size() == aPo2D.size());
  assert(aFlagPnt.size() == aPo2D.size());
  assert(aFlagTri.size() == aTri.size());

  double ratio = 3.0;
  for (;;) {
    int nadd = 0;
    for (unsigned int itri = 0; itri < aTri.size(); itri++) {
      const double area = Area_Tri(aVec2[aTri[itri].v[0]],
                                   aVec2[aTri[itri].v[1]],
                                   aVec2[aTri[itri].v[2]]);
      const double pcnt[2] = {
          (aVec2[aTri[itri].v[0]].x + aVec2[aTri[itri].v[1]].x + aVec2[aTri[itri].v[2]].x) / 3.0,
          (aVec2[aTri[itri].v[0]].y + aVec2[aTri[itri].v[1]].y + aVec2[aTri[itri].v[2]].y) / 3.0
      };
      double len2 = len * mesh_density.edgeLengthRatio(pcnt[0], pcnt[1]);
      if (area < len2 * len2 * ratio) { continue; }
      const auto ipo0 = static_cast<unsigned int>(aPo2D.size());
      aPo2D.resize(aPo2D.size() + 1);
      aVec2.resize(aVec2.size() + 1);
      aVec2[ipo0].p[0] = (aVec2[aTri[itri].v[0]].x + aVec2[aTri[itri].v[1]].x + aVec2[aTri[itri].v[2]].x) / 3.0;
      aVec2[ipo0].p[1] = (aVec2[aTri[itri].v[0]].y + aVec2[aTri[itri].v[1]].y + aVec2[aTri[itri].v[2]].y) / 3.0;
      InsertPoint_Elem(ipo0, itri, aPo2D, aTri);
      const int iflgtri = aFlagTri[itri];
      aFlagTri.push_back(iflgtri);
      aFlagTri.push_back(iflgtri);
      aFlagPnt.push_back(iflgtri + nflgpnt_offset);
      DelaunayAroundPoint(ipo0, aPo2D, aTri, aVec2);
      nadd++;
    }
    for (size_t ip = nPointFix; ip < aVec2.size(); ++ip) {
      dtri2::LaplacianArroundPoint(
          aVec2,
          static_cast<unsigned int>(ip),
          aPo2D, aTri);
    }
    if (nadd != 0) { ratio *= 0.8; }
    else { ratio *= 0.5; }
    if (ratio < 0.65) break;
  }

  for (size_t ip = nPointFix; ip < aVec2.size(); ++ip) {
    dtri2::LaplacianArroundPoint(
        aVec2,
        static_cast<unsigned int>(ip),
        aPo2D, aTri);
    DelaunayAroundPoint(
        static_cast<unsigned int>(ip),
        aPo2D, aTri, aVec2);
  }
}

DFM2_INLINE void delfem2::MakeSuperTriangle
    (std::vector<CVec2d> &aVec2,
     std::vector<CDynPntSur> &aPo2D,
     std::vector<CDynTri> &aTri,
     const double bound_2d[4]) { // super triangle
  assert(aVec2.size() == aPo2D.size());
  double max_len;
  double center[2];
  {
    max_len =
        (bound_2d[1] - bound_2d[0] > bound_2d[3] - bound_2d[2]) ? bound_2d[1] - bound_2d[0] : bound_2d[3] - bound_2d[2];
    center[0] = (bound_2d[1] + bound_2d[0]) * 0.5;
    center[1] = (bound_2d[3] + bound_2d[2]) * 0.5;
  }

  const double tri_len = max_len * 4.0;
  const double tmp_len = tri_len * sqrt(3.0) / 6.0;

  const size_t npo = aPo2D.size();
  aPo2D.resize(npo + 3);
  aVec2.resize(npo + 3);
  aVec2[npo + 0] = CVec2d(center[0], center[1] + 2.0 * tmp_len);
  aPo2D[npo + 0].e = 0;
  aPo2D[npo + 0].d = 0;
  aVec2[npo + 1] = CVec2d(center[0] - 0.5 * tri_len, center[1] - tmp_len);
  aPo2D[npo + 1].e = 0;
  aPo2D[npo + 1].d = 1;
  aVec2[npo + 2] = CVec2d(center[0] + 0.5 * tri_len, center[1] - tmp_len);
  aPo2D[npo + 2].e = 0;
  aPo2D[npo + 2].d = 2;

  aTri.resize(1);
  CDynTri &tri = aTri[0];
  tri.v[0] = static_cast<unsigned int>(npo + 0);
  tri.v[1] = static_cast<unsigned int>(npo + 1);
  tri.v[2] = static_cast<unsigned int>(npo + 2);
  tri.s2[0] = UINT_MAX;
  tri.s2[1] = UINT_MAX;
  tri.s2[2] = UINT_MAX;
}

DFM2_INLINE void delfem2::AddPointsMesh
    (const std::vector<CVec2d> &aVec2,
     std::vector<CDynPntSur> &aPo2D,
     std::vector<CDynTri> &aTri,
     int ipoin,
     double MIN_TRI_AREA) {
  assert(aPo2D.size() == aVec2.size());
  if (aPo2D[ipoin].e != UINT_MAX) { return; } // already added
  const CVec2d &po_add = aVec2[ipoin];
  int itri_in = -1;
  int iedge = -1;
  int iflg1 = 0, iflg2 = 0;
  for (unsigned int itri = 0; itri < aTri.size(); itri++) {
    iflg1 = 0;
    iflg2 = 0;
    if (Area_Tri(po_add, aVec2[aTri[itri].v[1]], aVec2[aTri[itri].v[2]]) > MIN_TRI_AREA) {
      iflg1++;
      iflg2 += 0;
    }
    if (Area_Tri(po_add, aVec2[aTri[itri].v[2]], aVec2[aTri[itri].v[0]]) > MIN_TRI_AREA) {
      iflg1++;
      iflg2 += 1;
    }
    if (Area_Tri(po_add, aVec2[aTri[itri].v[0]], aVec2[aTri[itri].v[1]]) > MIN_TRI_AREA) {
      iflg1++;
      iflg2 += 2;
    }
    if (iflg1 == 3) { // add in triangle
      itri_in = itri;
      break;
    } else if (iflg1 == 2) { // add in edge
      const int ied0 = 3 - iflg2;
      const unsigned int ipo_e0 = aTri[itri].v[(ied0 + 1) % 3];
      const unsigned int ipo_e1 = aTri[itri].v[(ied0 + 2) % 3];
      const unsigned int itri_s = aTri[itri].s2[ied0];
      if (itri_s == UINT_MAX) { return; }
      const unsigned int jno0 = FindAdjEdgeIndex(aTri[itri], ied0, aTri);
      assert(aTri[itri_s].v[(jno0 + 2) % 3] == ipo_e0);
      assert(aTri[itri_s].v[(jno0 + 1) % 3] == ipo_e1);
      const unsigned int inoel_d = jno0;
      assert(aTri[itri_s].s2[inoel_d] == itri);
      const unsigned int ipo_d = aTri[itri_s].v[inoel_d];
      assert(Area_Tri(po_add, aVec2[ipo_e1], aVec2[aTri[itri].v[ied0]]) > MIN_TRI_AREA);
      assert(Area_Tri(po_add, aVec2[aTri[itri].v[ied0]], aVec2[ipo_e0]) > MIN_TRI_AREA);
      if (Area_Tri(po_add, aVec2[ipo_e0], aVec2[ipo_d]) < MIN_TRI_AREA) { continue; }
      if (Area_Tri(po_add, aVec2[ipo_d], aVec2[ipo_e1]) < MIN_TRI_AREA) { continue; }
      const int det_d = DetDelaunay(po_add, aVec2[ipo_e0], aVec2[ipo_e1], aVec2[ipo_d]);
      if (det_d == 2 || det_d == 1) continue;
      itri_in = itri;
      iedge = ied0;
      break;
    }
  }
  if (itri_in == -1) {
    std::cout << "super triangle failure " << iflg1 << " " << iflg2 << std::endl;
    assert(0);
    abort();
  }
  if (iedge == -1) {
    InsertPoint_Elem(ipoin, itri_in, aPo2D, aTri);
  } else {
    InsertPoint_ElemEdge(ipoin, itri_in, iedge, aPo2D, aTri);
  }
}

DFM2_INLINE void delfem2::EnforceEdge
    (std::vector<CDynPntSur> &aPo2D,
     std::vector<CDynTri> &aTri,
     unsigned int ip0,
     unsigned int ip1,
     const std::vector<CVec2d> &aVec2) { // enforce edge
  assert(ip0 < aPo2D.size());
  assert(ip1 < aPo2D.size());
  for (;;) {
    unsigned int itri0, inotri0, inotri1;
    if (FindEdge_LookAroundPoint(itri0,
                                 inotri0,
                                 inotri1,
                                 ip0,
                                 ip1,
                                 aPo2D,
                                 aTri)) { // this edge divide outside and inside
      assert(inotri0 != inotri1);
      assert(inotri0 < 3);
      assert(inotri1 < 3);
      assert(aTri[itri0].v[inotri0] == ip0);
      assert(aTri[itri0].v[inotri1] == ip1);
      const int ied0 = 3 - inotri0 - inotri1;
      {
        const unsigned int itri1 = aTri[itri0].s2[ied0];
        const unsigned int ied1 = FindAdjEdgeIndex(aTri[itri0], ied0, aTri);
        assert(aTri[itri1].s2[ied1] == itri0);
        aTri[itri1].s2[ied1] = UINT_MAX;
        aTri[itri0].s2[ied0] = UINT_MAX;
      }
      break;
    } else { // this edge is devided from connection outer triangle
      double ratio;
      if (!dtri2::FindEdgePoint_AcrossEdge(itri0, inotri0, inotri1, ratio,
                                           ip0, ip1,
                                           aPo2D, aTri, aVec2)) { assert(0); }
      assert(ratio > -1.0e-20 && ratio < 1.0 + 1.0e-20);
      assert(Area_Tri(aVec2[ip0], aVec2[aTri[itri0].v[inotri0]], aVec2[ip1]) > 1.0e-20);
      assert(Area_Tri(aVec2[ip0], aVec2[ip1], aVec2[aTri[itri0].v[inotri1]]) > 1.0e-20);
      //            std::cout << ratio << std::endl;
      if (ratio < 1.0e-20) {
        assert(0);
        abort();
      } else if (ratio > 1.0 - 1.0e-10) {
        assert(0);
        abort();
      } else {
        const int ied0 = 3 - inotri0 - inotri1;
        assert(aTri[itri0].s2[ied0] < aTri.size());
#if !defined(NDEBUG)
        const unsigned int itri1 = aTri[itri0].s2[ied0];
        const unsigned int ied1 = FindAdjEdgeIndex(aTri[itri0], ied0, aTri);
        assert(aTri[itri1].s2[ied1] == itri0);
#endif
        bool res = FlipEdge(itri0, ied0, aPo2D, aTri);
//        std::cout << itri0 << " " << ied0 << " " << ratio << " " << res << std::endl;
//        continue;
        if (!res) {
          break;
        }
      }
    }
  }
}

DFM2_INLINE void delfem2::FlagConnected
    (std::vector<int> &inout_flg,
     const std::vector<CDynTri> &aTri_in,
     unsigned int itri0_ker,
     int iflag) {
#ifndef NDEBUG
  const size_t ntri = aTri_in.size();
  assert(inout_flg.size() == ntri);
  assert(itri0_ker < inout_flg.size());
#endif
  inout_flg[itri0_ker] = iflag;
  std::stack<unsigned int> ind_stack;
  ind_stack.push(itri0_ker);
  for (;;) {
    if (ind_stack.empty()) break;
    const int itri_cur = ind_stack.top();
    ind_stack.pop();
    for (unsigned int jtri0 : aTri_in[itri_cur].s2) {
      if (jtri0 == UINT_MAX) continue;
      if (inout_flg[jtri0] != iflag) {
        inout_flg[jtri0] = iflag;
        ind_stack.push(jtri0);
      }
    }
  }
}

DFM2_INLINE void delfem2::DeleteTriFlag
    (std::vector<CDynTri> &aTri1,
     std::vector<int> &aFlg1,
     int iflag) {
  assert(aFlg1.size() == aTri1.size());
  const size_t ntri0 = aTri1.size();
  std::vector<unsigned int> map01(ntri0, UINT_MAX);
  unsigned int ntri1 = 0;
  for (unsigned int itri = 0; itri < ntri0; ++itri) {
    if (aFlg1[itri] != iflag) {
      map01[itri] = ntri1;
      ntri1++;
    }
  }
  const std::vector<CDynTri> aTri0 = aTri1;
  const std::vector<int> aFlg0 = aFlg1;
  aTri1.clear();
  aTri1.resize(ntri1);
  aFlg1.resize(ntri1);
  for (unsigned int itri0 = 0; itri0 < aTri0.size(); itri0++) {
    if (map01[itri0] != UINT_MAX) {
      const unsigned int itri1 = map01[itri0];
      assert(itri1 < ntri1);
      aTri1[itri1] = aTri0[itri0];
      aFlg1[itri1] = aFlg0[itri0];
      assert(aFlg1[itri1] != iflag);
    }
  }
  for (unsigned int itri1 = 0; itri1 < ntri1; itri1++) {
    for (unsigned int ifatri = 0; ifatri < 3; ++ifatri) {
      if (aTri1[itri1].s2[ifatri] == UINT_MAX) continue;
      const unsigned int itri_s0 = aTri1[itri1].s2[ifatri];
      assert(itri_s0 < aTri0.size());
      unsigned int jtri0 = map01[itri_s0];
      assert(jtri0 < aTri1.size());
      aTri1[itri1].s2[ifatri] = jtri0;
    }
  }
}

DFM2_INLINE void delfem2::DeleteUnrefPoints(
    std::vector<CVec2d> &aVec2,
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CDynTri> &aTri_in,
    const std::vector<unsigned int> &aPoDel) {
  assert(aPo2D.size() == aVec2.size());
  std::vector<unsigned int> map_po_del;
  unsigned int npo_pos;
  {
    map_po_del.resize(aPo2D.size(), UINT_MAX - 1);
    for (unsigned int ipo : aPoDel) {
      map_po_del[ipo] = UINT_MAX;
    }
    npo_pos = 0;
    for (unsigned int ipo = 0; ipo < aPo2D.size(); ipo++) {
      if (map_po_del[ipo] == UINT_MAX) continue;
      map_po_del[ipo] = npo_pos;
      npo_pos++;
    }
  }
  {
    std::vector<CDynPntSur> aPo_tmp = aPo2D;
    std::vector<CVec2d> aVec2_tmp = aVec2;
    aPo2D.resize(npo_pos);
    aVec2.resize(npo_pos);
    for (unsigned int ipo = 0; ipo < map_po_del.size(); ipo++) {
      if (map_po_del[ipo] == UINT_MAX) continue;
      unsigned int ipo1 = map_po_del[ipo];
      aPo2D[ipo1] = aPo_tmp[ipo];
      aVec2[ipo1] = aVec2_tmp[ipo];
    }
  }
  for (unsigned int itri = 0; itri < aTri_in.size(); itri++) {
    for (int ifatri = 0; ifatri < 3; ifatri++) {
      const unsigned int ipo = aTri_in[itri].v[ifatri];
      assert(map_po_del[ipo] != UINT_MAX);
      aTri_in[itri].v[ifatri] = map_po_del[ipo];
      aPo2D[ipo].e = itri;
      aPo2D[ipo].d = ifatri;
    }
  }
}

DFM2_INLINE void delfem2::DeletePointsFlag(
    std::vector<CVec2d> &aVec1,
    std::vector<CDynPntSur> &aPo1,
    std::vector<CDynTri> &aTri,
    std::vector<int> &aFlgPnt1,
    int iflg) {
  const size_t np0 = aVec1.size();
  assert(aPo1.size() == np0);
  assert(aFlgPnt1.size() == np0);
  std::vector<int> map01;
  int npo1;
  {
    map01.resize(np0, -1);
    npo1 = 0;
    for (unsigned int ipo = 0; ipo < np0; ipo++) {
      if (aFlgPnt1[ipo] == iflg) continue;
      map01[ipo] = npo1;
      npo1++;
    }
  }
  {
    std::vector<CDynPntSur> aPo0 = aPo1;
    std::vector<CVec2d> aVec0 = aVec1;
    std::vector<int> aFlgPnt0 = aFlgPnt1;
    aPo1.resize(npo1);
    aVec1.resize(npo1);
    aFlgPnt1.resize(npo1);
    for (unsigned int ipo0 = 0; ipo0 < np0; ipo0++) {
      if (map01[ipo0] == -1) continue;
      int ipo1 = map01[ipo0];
      assert(ipo1 >= 0 && ipo1 < npo1);
      aPo1[ipo1] = aPo0[ipo0];
      aVec1[ipo1] = aVec0[ipo0];
      aFlgPnt1[ipo1] = aFlgPnt0[ipo0];
    }
  }
  for (unsigned int itri = 0; itri < aTri.size(); itri++) {
    for (int ifatri = 0; ifatri < 3; ifatri++) {
      const unsigned int ipo = aTri[itri].v[ifatri];
      aTri[itri].v[ifatri] = map01[ipo];
      aPo1[ipo].e = itri;
      aPo1[ipo].d = ifatri;
    }
  }
}

DFM2_INLINE void delfem2::Meshing_Initialize(
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CDynTri> &aTri,
    std::vector<CVec2d> &aVec2) {
  aPo2D.resize(aVec2.size());
  for (size_t ixys = 0; ixys < aVec2.size(); ixys++) {
    aPo2D[ixys].e = UINT_MAX;
    aPo2D[ixys].d = 0;
  }
  {
    double bound_2d[4];
    bound_2d[0] = aVec2[0].x;
    bound_2d[1] = aVec2[0].x;
    bound_2d[2] = aVec2[0].y;
    bound_2d[3] = aVec2[0].y;
    for (int ipoin = 1; ipoin < (int) aPo2D.size(); ipoin++) {
      if (aVec2[ipoin].x < bound_2d[0]) { bound_2d[0] = aVec2[ipoin].x; }
      if (aVec2[ipoin].x > bound_2d[1]) { bound_2d[1] = aVec2[ipoin].x; }
      if (aVec2[ipoin].y < bound_2d[2]) { bound_2d[2] = aVec2[ipoin].y; }
      if (aVec2[ipoin].y > bound_2d[3]) { bound_2d[3] = aVec2[ipoin].y; }
    }
    MakeSuperTriangle(aVec2, aPo2D, aTri,
                      bound_2d);
  }
  {
    const double MIN_TRI_AREA = 1.0e-10;
    for (size_t ip = 0; ip < aPo2D.size() - 3; ++ip) {
      AddPointsMesh(
          aVec2, aPo2D, aTri,
          static_cast<unsigned int>(ip),
          MIN_TRI_AREA);
      DelaunayAroundPoint(
          static_cast<unsigned int>(ip),
          aPo2D, aTri, aVec2);
    }
  }
}



/*
 bool Triangulation
 (std::vector<int>& aTri_out,    // out
 std::vector<double>& aVec, // out
 const std::vector<int>& aPtrVtxInd, // out
 const std::vector<int>& aVtxInd, // out
 ////
 const double max_edge_length, // ind
 const CInputTriangulation& mesh_density)
 {
 
 const int nxy = aXY_in.size()/2;
 std::vector<CEPo2> aPo2D;
 std::vector<CVector2> aVec2;
 aPo2D.resize(nxy);
 aVec2.resize(nxy);
 for(int ixys=0;ixys<nxy;ixys++){
 aVec2[ixys] = CVector2(aXY_in[ixys*2+0], aXY_in[ixys*2+1]);
 aPo2D[ixys].e = -1;
 aPo2D[ixys].d = -1;
 }
 
 ////////////////////////////////
 std::vector<ETri> aTri_in;
 if( !MeshingOuterLoop(aPo2D,aTri_in,aVec2,    aPtrVtxInd, aVtxInd) ){
 return true;
 }
 if( max_edge_length > 0 ){
 MeshingInside2(aPo2D,aTri_in,aVec2, aVtxInd,max_edge_length,mesh_density);
 }
 
 ////////////////////////////////
 // pushing back to STL vector
 
 return true;
 }
 */

DFM2_INLINE void delfem2::MeshTri2D_Export
    (std::vector<double> &aXY_out,
     std::vector<unsigned int> &aTri_out,
     const std::vector<CVec2d> &aVec2,
     const std::vector<CDynTri> &aTri_in) {
  aTri_out.clear();
  aXY_out.clear();
  const size_t ntri = aTri_in.size();
  aTri_out.resize(ntri * 3);
  for (unsigned int itri = 0; itri < ntri; itri++) {
    aTri_out[itri * 3 + 0] = aTri_in[itri].v[0];
    aTri_out[itri * 3 + 1] = aTri_in[itri].v[1];
    aTri_out[itri * 3 + 2] = aTri_in[itri].v[2];
  }
  const size_t nxy_out = aVec2.size();
  aXY_out.resize(nxy_out * 2);
  for (unsigned int ixy = 0; ixy < nxy_out; ixy++) {
    aXY_out[ixy * 2 + 0] = aVec2[ixy].x;
    aXY_out[ixy * 2 + 1] = aVec2[ixy].y;
  }
}

// -------------------------------------------------------

/*
 ////////////////////////////////
 // resampling
 // no resampling edge if(max_edge_length < 0)
 std::vector< std::vector<int> > aPoInEd;
 aPoInEd.resize(nxys_presum);
 if( max_edge_length > 0 ){
 for(int iloop=0;iloop<nloop;++iloop){
 int nadd = 0;
 const int nbar = loop_ind[iloop+1]-loop_ind[iloop];
 for(int ibar=0;ibar<nbar;ibar++){
 int ipo0 = loop_ind[iloop]+ibar;
 int ipo1 = loop_ind[iloop]+ibar+1;
 if( ibar == nbar-1 ){ ipo1 = loop_ind[iloop]; }
 const double len = Distance( aVec2[ipo0], aVec2[ipo1] );
 nadd = (int)(len / max_edge_length);
 if( nadd == 0 || !is_add_point_boundary ) continue;
 const int ndiv = nadd+1;
 const double delx = (aVec2[ipo1].x - aVec2[ipo0].x)/ndiv;
 const double dely = (aVec2[ipo1].y - aVec2[ipo0].y)/ndiv;
 for(int iadd=0;iadd<nadd;++iadd){
 const unsigned int ipo = (int)aPo2D.size();
 CVector2 v2;
 v2.x = aVec2[ipo0].x + delx*(iadd+1);
 v2.y = aVec2[ipo0].y + dely*(iadd+1);
 CEPo2 po;
 po.e = -1;
 po.d = -1;
 aPo2D.push_back(po);
 aVec2.push_back(v2);
 aPoInEd[ loop_ind[iloop]+ibar ].push_back(ipo);
 }
 }
 }
 }
 
 ////////////////////////////////
 // invert the index if negative area
 {
 aPtrVtxInd.resize(nloop+1);
 aPtrVtxInd[0] = 0;
 for(int iloop=0;iloop<nloop;++iloop){
 const int nbar0 = loop_ind[iloop+1]-loop_ind[iloop];
 int nbar1 = nbar0;
 for(int ibar=0;ibar<nbar0;ibar++){
 nbar1 += aPoInEd[ loop_ind[iloop]+ibar].size();
 }
 aPtrVtxInd[iloop+1] = aPtrVtxInd[iloop] + nbar1;
 }
 // adding new vertices on the outline
 aVtxInd.resize(aPtrVtxInd[nloop]);
 {
 int ivtx0 = 0;
 for(int iloop=0;iloop<nloop;iloop++){
 double area_loop = 0;
 { // area of this loop
 CVector2 vtmp(0,0);
 const int nbar = aPtrVtxInd[iloop+1]-aPtrVtxInd[iloop];
 for(int ibar=0;ibar<nbar;ibar++){
 int ipo0 = aPtrVtxInd[iloop]+ibar;
 int ipo1 = aPtrVtxInd[iloop]+ibar+1;
 if( ibar == nbar-1 ){ ipo1 = aPtrVtxInd[iloop]; }
 area_loop += TriArea( vtmp, aVec2[ipo0], aVec2[ipo1] );
 }
 }
 const int nbar0 = loop_ind[iloop+1]-loop_ind[iloop];
 if( (area_loop > 0) == (iloop == 0) ){ // outer loop
 for(int ibar=0;ibar<nbar0;ibar++){
 int ie = loop_ind[iloop] + ibar;
 const std::vector<int>& add = aPoInEd[ie];
 aVtxInd[ivtx0] = ie;  ivtx0++;
 for(unsigned int iadd=0;iadd<add.size();iadd++){
 aVtxInd[ivtx0] = add[iadd];  ivtx0++;
 }
 }
 }
 else{
 for(int ibar=0;ibar<nbar0;ibar++){ // inner loop
 int ie = loop_ind[iloop+1] - 1 - ibar;
 const std::vector<int>& add = aPoInEd[ie];
 const int nadd = (int)add.size();
 for(int iadd=0;iadd<nadd;iadd++){
 aVtxInd[ivtx0] = add[nadd-1-iadd];  ivtx0++;
 }
 aVtxInd[ivtx0] = ie;  ivtx0++;
 }
 }
 }
 }
 }
 */

/*
void PrepareInput
(std::vector<int>& loop1_ind, // out
 std::vector<int>& loop1,
 const std::vector<int>& loop0_ind,
 const std::vector<double>& aXY)
{
  const int nloop = loop0_ind.size()-1;
  const int nxy = aXY.size()/2;
  std::vector< std::vector<int> > aPoInEd;
  aPoInEd.resize(nxy);
  loop1_ind.resize(nloop+1);
  loop1_ind[0] = 0;
  for(int iloop=0;iloop<nloop;++iloop){
    const int nbar0 = loop0_ind[iloop+1]-loop0_ind[iloop];
    int nbar1 = nbar0;
    for(int ibar=0;ibar<nbar0;ibar++){
      nbar1 += aPoInEd[ loop0_ind[iloop]+ibar].size();
    }
    loop1_ind[iloop+1] = loop1_ind[iloop] + nbar1;
  }
  // adding new vertices on the outline
  loop1.resize(loop1_ind[nloop]);
  {
    int ivtx0 = 0;
    for(int iloop=0;iloop<nloop;iloop++){
      double area_loop = 0;
      { // area of this loop
        CVector2 vtmp(0,0);
        const int nbar = loop1_ind[iloop+1]-loop1_ind[iloop];
        for(int ibar=0;ibar<nbar;ibar++){
          int ipo0 = loop1_ind[iloop]+ibar;
          int ipo1 = loop1_ind[iloop]+ibar+1;
          if( ibar == nbar-1 ){ ipo1 = loop1_ind[iloop]; }
          area_loop += TriArea(vtmp,
                               CVector2(aXY[ipo0*2+0],aXY[ipo0*2+1]),
                               CVector2(aXY[ipo1*2+0],aXY[ipo0*2+1]) );
        }
      }
      const int nbar0 = loop0_ind[iloop+1]-loop0_ind[iloop];
      if( (area_loop > 0) == (iloop == 0) ){ // outer loop
        for(int ibar=0;ibar<nbar0;ibar++){
          int ie = loop0_ind[iloop] + ibar;
          const std::vector<int>& add = aPoInEd[ie];
          loop1[ivtx0] = ie;  ivtx0++;
          for(unsigned int iadd=0;iadd<add.size();iadd++){
            loop1[ivtx0] = add[iadd];  ivtx0++;
          }
        }
      }
      else{
        for(int ibar=0;ibar<nbar0;ibar++){ // inner loop
          int ie = loop0_ind[iloop+1] - 1 - ibar;
          const std::vector<int>& add = aPoInEd[ie];
          const int nadd = (int)add.size();
          for(int iadd=0;iadd<nadd;iadd++){
            loop1[ivtx0] = add[nadd-1-iadd];  ivtx0++;
          }
          loop1[ivtx0] = ie;  ivtx0++;
        }
      }
    }
  }
}
 */


DFM2_INLINE void delfem2::Meshing_SingleConnectedShape2D(
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CVec2d> &aVec2,
    std::vector<CDynTri> &aDTri,
    const std::vector<int> &loopIP_ind,
    const std::vector<int> &loopIP) {
  std::vector<unsigned int> aPoDel;
  {
    const size_t npo = aVec2.size();
    aPoDel.push_back(static_cast<unsigned int>(npo + 0));
    aPoDel.push_back(static_cast<unsigned int>(npo + 1));
    aPoDel.push_back(static_cast<unsigned int>(npo + 2));
  }
  Meshing_Initialize(aPo2D, aDTri, aVec2);
  for (size_t iloop = 0; iloop < loopIP_ind.size() - 1; ++iloop) {
    const int np0 = loopIP_ind[iloop + 1] - loopIP_ind[iloop];
    for (int iip = loopIP_ind[iloop]; iip < loopIP_ind[iloop + 1]; ++iip) {
      const int ip0 = loopIP[loopIP_ind[iloop] + (iip + 0) % np0];
      const int ip1 = loopIP[loopIP_ind[iloop] + (iip + 1) % np0];
      EnforceEdge(aPo2D, aDTri,
                  ip0, ip1, aVec2);
    }
  }
  { // delete triangles
    std::vector<int> aflg(aDTri.size(), 0);
    unsigned int itri0_ker = UINT_MAX;
    {
      unsigned int iedtri;
      FindEdge_LookAllTriangles(itri0_ker, iedtri,
                                loopIP[0], loopIP[1], aDTri);
      assert(itri0_ker < aDTri.size());
      FlagConnected(aflg,
                    aDTri, itri0_ker, 1);
    }
    DeleteTriFlag(aDTri,
                  aflg, 0);
  }
  // delete points
  DeleteUnrefPoints(aVec2, aPo2D, aDTri,
                    aPoDel);
}


// -----------------------------------------

DFM2_INLINE void delfem2::CMeshTri2D(
    std::vector<double> &aXY,
    std::vector<unsigned int> &aTri,
    const std::vector<CVec2d> &aVec2,
    const std::vector<CDynTri> &aETri) {
  aXY.resize(aVec2.size() * 2);
  for (size_t ip = 0; ip < aVec2.size(); ++ip) {
    aXY[ip * 2 + 0] = aVec2[ip].x;
    aXY[ip * 2 + 1] = aVec2[ip].y;
  }
  aTri.resize(aETri.size() * 3);
  for (size_t it = 0; it < aETri.size(); ++it) {
    aTri[it * 3 + 0] = aETri[it].v[0];
    aTri[it * 3 + 1] = aETri[it].v[1];
    aTri[it * 3 + 2] = aETri[it].v[2];
  }
}

DFM2_INLINE void delfem2::RefinementPlan_EdgeLongerThan_InsideCircle(
    CCmdRefineMesh &aCmd,
    double elen,
    double px,
    double py,
    double rad,
    [[maybe_unused]] const std::vector<CDynPntSur> &aPo2D,
    const std::vector<CVec2d> &aVec2,
    const std::vector<CDynTri> &aETri) {
  std::set<CCmdRefineMesh::CCmdEdge> setCmd;
  for (const auto &itri : aETri) {
    const unsigned int i0 = itri.v[0];
    const unsigned int i1 = itri.v[1];
    const unsigned int i2 = itri.v[2];
    const CVec2d &p0 = aVec2[i0];
    const CVec2d &p1 = aVec2[i1];
    const CVec2d &p2 = aVec2[i2];
    CVec2d pc = (p0 + p1 + p2) * 0.333333;
    CVec2d ps = CVec2d(px, py);
    if (Distance(pc, ps) < rad) {
      const double d01 = Distance(p0, p1);
      const double d12 = Distance(p1, p2);
      const double d20 = Distance(p2, p0);
      if (d01 > elen) { setCmd.insert(CCmdRefineMesh::CCmdEdge(i0, i1, 0.5)); }
      if (d12 > elen) { setCmd.insert(CCmdRefineMesh::CCmdEdge(i1, i2, 0.5)); }
      if (d20 > elen) { setCmd.insert(CCmdRefineMesh::CCmdEdge(i2, i0, 0.5)); }
    }
  }
  aCmd.aCmdEdge.assign(setCmd.begin(), setCmd.end());
}


// TODO: implement this function
DFM2_INLINE void delfem2::RefineMesh(
    std::vector<CDynPntSur> &aEPo2,
    std::vector<CDynTri> &aSTri,
    std::vector<CVec2d> &aVec2,
    CCmdRefineMesh &aCmd) {
  assert(aVec2.size() == aEPo2.size());
  std::stack<int> aIV_free;
  for (unsigned int ip = 0; ip < aEPo2.size(); ++ip) {
    if (aEPo2[ip].e != UINT_MAX) { continue; }
    aIV_free.push(ip);
  }
  for (auto &cmd : aCmd.aCmdEdge) {
    int i0 = cmd.ipo0;
    int i1 = cmd.ipo1;
    double r0 = cmd.r0;
    CVec2d v01 = r0 * aVec2[i0] + (1.0 - r0) * aVec2[i1];
    if (aIV_free.empty()) {
      auto ipo = static_cast<unsigned int>(aVec2.size());
      aVec2.push_back(v01);
      aEPo2.emplace_back();
      cmd.ipo_new = ipo;
    } else {
      int ipo = aIV_free.top();
      aIV_free.pop();
      aVec2[ipo] = v01;
      aEPo2[ipo] = CDynPntSur();
      cmd.ipo_new = ipo;
    }
  }
  for (auto &icmd : aCmd.aCmdEdge) {
    const int ip0 = icmd.ipo_new;
    AddPointsMesh(aVec2, aEPo2, aSTri, ip0, 1.0e-10);
    DelaunayAroundPoint(ip0, aEPo2, aSTri, aVec2);
  }
}

DFM2_INLINE void delfem2::MakeInvMassLumped_Tri(
    std::vector<double> &aInvMassLumped,
    double rho,
    const std::vector<CVec2d> &aVec2,
    const std::vector<CDynTri> &aETri) {
  aInvMassLumped.assign(aVec2.size(), 0.0);
  for (const auto &it : aETri) {
    const unsigned int aIP[3] = {it.v[0], it.v[1], it.v[2]};
    double P[3][2] = {
        {aVec2[aIP[0]].x, aVec2[aIP[0]].y},
        {aVec2[aIP[1]].x, aVec2[aIP[1]].y},
        {aVec2[aIP[2]].x, aVec2[aIP[2]].y}};
    const double Area = Area_Tri2(P[0], P[1], P[2]);
    for (auto ip : aIP) {
      aInvMassLumped[ip] += Area * rho / 3.0;
    }
  }
  for (size_t ip = 0; ip < aVec2.size(); ++ip) {
    double m0 = aInvMassLumped[ip];
    if (m0 < 1.0e-10) { continue; }
    aInvMassLumped[ip] = 1.0 / m0;
  }
}

DFM2_INLINE void delfem2::MinMaxTriArea(
    double &min_area,
    double &max_area,
    const std::vector<CVec2d> &aVec2,
    const std::vector<CDynTri> &aETri) {
  for (size_t it = 0; it < aETri.size(); ++it) {
    const unsigned int i0 = aETri[it].v[0];
    const unsigned int i1 = aETri[it].v[1];
    const unsigned int i2 = aETri[it].v[2];
    double P[3][2] = {
        {aVec2[i0].x, aVec2[i0].y},
        {aVec2[i1].x, aVec2[i1].y},
        {aVec2[i2].x, aVec2[i2].y}};
    const double Area = Area_Tri2(P[0], P[1], P[2]);
    if (it == 0) {
      min_area = max_area = Area;
      continue;
    }
    if (Area < min_area) { min_area = Area; }
    if (Area > max_area) { max_area = Area; }
  }
}

DFM2_INLINE void delfem2::GenMesh(
    std::vector<CDynPntSur> &aPo2D,
    std::vector<CDynTri> &aETri,
    std::vector<CVec2d> &aVec2,
    const std::vector<std::vector<double> > &aaXY,
    double resolution_edge,
    double resolution_face) {
  std::vector<int> loopIP_ind, loopIP;
  {
    JArray_FromVecVec_XY(loopIP_ind, loopIP, aVec2,
                         aaXY);
    if (!CheckInputBoundaryForTriangulation(loopIP_ind, aVec2)) {
      return;
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind, aVec2);
    if (resolution_edge > 10e-10) {
      ResamplingLoop(loopIP_ind, loopIP, aVec2,
                     resolution_edge);
    }
  }
  // -----------------------
  Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                 loopIP_ind, loopIP);
  if (resolution_face > 1.0e-10) {
    CInputTriangulation_Uniform param(1.0);
    std::vector<int> flg_pnt(aVec2.size());
    std::vector<unsigned int> flg_tri(aETri.size(), 0);
    MeshingInside(aPo2D, aETri, aVec2, flg_pnt, flg_tri,
                  aVec2.size(), 0, resolution_face, param);
  }
}
