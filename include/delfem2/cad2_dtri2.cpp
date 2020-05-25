/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/cad2_dtri2.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/funcs.h"
#include <cstdio>
#include <deque>
#include <climits>


// ------------------------------------------

namespace delfem2 {
namespace cad2 {

DFM2_INLINE delfem2::CBoundingBox2D BB_LoopEdgeCad2D
    (const std::vector<CCad2D_EdgeGeo> &aEdge) {
  CBoundingBox2D bb;
  for (const auto &ie : aEdge) {
    CBoundingBox2D bb0 = ie.BB();
    bb += bb0;
  }
  return bb;
}

DFM2_INLINE void GetBound
    (double bound_2d[4],
     unsigned int ifc0,
     const CCadTopo &topo,
     const std::vector<CCad2D_EdgeGeo> &aEdgeGeo) {
  assert(ifc0 < topo.aFace.size());
  const int il0 = topo.aFace[ifc0].aIL[0];
  const std::vector<std::pair<int, bool> > &aIE = topo.aLoop[il0].aIE;
  {
    int ie0 = aIE[0].first;
    const CVec2d &p0 = aEdgeGeo[ie0].p0;
    bound_2d[0] = p0.x();
    bound_2d[1] = p0.x();
    bound_2d[2] = p0.y();
    bound_2d[3] = p0.y();
  }
  for (const auto &iie : aIE) {
    const unsigned int ie0 = iie.first;
    {
      const CVec2d &p0 = aEdgeGeo[ie0].p0;
      if (p0.x() < bound_2d[0]) { bound_2d[0] = p0.x(); }
      if (p0.x() > bound_2d[1]) { bound_2d[1] = p0.x(); }
      if (p0.y() < bound_2d[2]) { bound_2d[2] = p0.y(); }
      if (p0.y() > bound_2d[3]) { bound_2d[3] = p0.y(); }
    }
    for (const auto &p0 : aEdgeGeo[ie0].aP) {
      if (p0.x() < bound_2d[0]) { bound_2d[0] = p0.x(); }
      if (p0.x() > bound_2d[1]) { bound_2d[1] = p0.x(); }
      if (p0.y() < bound_2d[2]) { bound_2d[2] = p0.y(); }
      if (p0.y() > bound_2d[3]) { bound_2d[3] = p0.y(); }
    }
  }
}

DFM2_INLINE void GenMeshCadFace
    (std::vector<CVec2d> &aVec2,
     std::vector<CDynTri> &aETri,
     const CCad2D_FaceGeo &facegeo, unsigned int iface0,
     const CCadTopo &topo,
     const std::vector<CCad2D_VtxGeo> &aVtxGeo,
     const std::vector<CCad2D_EdgeGeo> &aEdgeGeo) {
  assert(iface0 < topo.aFace.size());
  std::vector<CDynPntSur> aPo2D;
  {
    aPo2D.resize(aVec2.size());
    for (size_t ixys = 0; ixys < aVec2.size(); ixys++) {
      aPo2D[ixys].e = UINT_MAX;
      aPo2D[ixys].d = 0;
    }
  }
  {
    double bound_2d[4];
    GetBound(bound_2d,
             iface0, topo, aEdgeGeo);
    MakeSuperTriangle(aVec2, aPo2D, aETri,
                      bound_2d);
  }
  std::vector<int> aIP;
  { // make list of index of point involved this face
    for (size_t iil = 0; iil < topo.aFace[iface0].aIL.size(); ++iil) {
      const int il0 = topo.aFace[iface0].aIL[iil];
      if (topo.aLoop[il0].iv != -1) {
        aIP.push_back(topo.aLoop[il0].iv);
        continue;
      }
      const std::vector<std::pair<int, bool> > &aIE = topo.aLoop[il0].aIE;
      for (const auto &iie : aIE) {
        const int ie = iie.first;
        aIP.push_back(topo.aEdge[ie].iv0);
        for (size_t iip = 0; iip < aEdgeGeo[ie].aP.size(); ++iip) {
          aIP.push_back(aEdgeGeo[ie].ip0 + iip);
        }
      }
    }
  }
  {
    const double MIN_TRI_AREA = 1.0e-10;
    for (int ip : aIP) {
      AddPointsMesh(aVec2, aPo2D, aETri,
                    ip, MIN_TRI_AREA);
      DelaunayAroundPoint(ip, aPo2D, aETri, aVec2);
    }
  }
  {
    for (size_t iil = 0; iil < topo.aFace[iface0].aIL.size(); ++iil) {
      const int il0 = topo.aFace[iface0].aIL[iil];
      if (topo.aLoop[il0].iv != -1) { continue; }
      const std::vector<std::pair<int, bool> > &aIE = topo.aLoop[il0].aIE;
      for (const auto &iie : aIE) {
        const int ie0 = iie.first;
        const unsigned int np = aEdgeGeo[ie0].aP.size();
        const unsigned int nseg = np + 1;
        for (unsigned int iseg = 0; iseg < nseg; ++iseg) {
          const int ip0 = (iseg == 0) ? topo.aEdge[ie0].iv0 : aEdgeGeo[ie0].ip0 + iseg - 1;
          const int ip1 = (iseg == nseg - 1) ? topo.aEdge[ie0].iv1 : aEdgeGeo[ie0].ip0 + iseg;
          EnforceEdge(aPo2D, aETri,
                      ip0, ip1, aVec2);
        }
      }
    }
  }
  std::vector<int> aFlgTri(aETri.size(), -1);
  {
    int ip0 = -1, ip1 = -1;
    {
      const int il0 = topo.aFace[iface0].aIL[0];
      const std::vector<std::pair<int, bool> > &aIE = topo.aLoop[il0].aIE;
      unsigned int ie0 = aIE[0].first;
      const unsigned int np = aEdgeGeo[ie0].aP.size();
      ip0 = topo.aEdge[ie0].iv0;
      ip1 = (np == 0) ? topo.aEdge[ie0].iv1 : aEdgeGeo[ie0].ip0;
    }
    assert(ip0 != -1 && ip1 != -1);
    unsigned int itri0_ker=UINT_MAX, iedtri;
    FindEdge_LookAllTriangles(itri0_ker, iedtri,
                              ip0, ip1, aETri);
    assert( itri0_ker < aETri.size() );
    FlagConnected(aFlgTri,
                  aETri, itri0_ker, (int) iface0);
  }
  { // Delete Outer Triangles & Points
    assert(aFlgTri.size() == aETri.size());
    DeleteTriFlag(aETri, aFlgTri,
                  -1);
    aPo2D.resize(aPo2D.size() - 3);
    aVec2.resize(aVec2.size() - 3);
  }
}

DFM2_INLINE std::vector<std::string> SVG_Split_Path_d
    (std::string &s0) {
  unsigned int imark = 0;
  std::vector<std::string> aS;
  for (size_t i = 0; i < s0.size(); ++i) {
    if (isNumeric(s0[i])) continue;
    if (s0[i] == ',') {
      std::string s1(s0.begin() + imark, s0.begin() + i);
      aS.push_back(s1);
      imark = i + 1;  // mark shoud be the begining position of the string so move next
    }
    if (s0[i] == ' ') { // sometimes the space act as delimiter in the SVG (inkscape version)
      if (i > imark) {
        std::string s1(s0.begin() + imark, s0.begin() + i);
        aS.push_back(s1);
      }
      imark = i+1; // mark shoud be the begining position of the string so move next
    }
    if (s0[i] == '-') {
      if (i > imark) {
        std::string s1(s0.begin() + imark, s0.begin() + i);
        aS.push_back(s1);
      }
      imark = i;
    }
    if (isAlphabet(s0[i])) {
      if (i > imark) {
        std::string s1(s0.begin() + imark, s0.begin() + i);
        aS.push_back(s1);
      }
      const char s2[2] = {s0[i], '\0'};
      aS.emplace_back(s2);
      imark = i + 1;
    }
  }
  return aS;
}

DFM2_INLINE void LoopEdgeCad2D_SVGPathD
    (std::vector<CCad2D_EdgeGeo> &aEdge,
     std::vector<std::string> &aStr1) {
  assert(aStr1[0] == "M");
  assert(aStr1[aStr1.size() - 1] == "Z" || aStr1[aStr1.size() - 1] == "z");
  CVec2d pos_cur;
  for (int is = 0;;) {
    if (aStr1[is] == "M") {
      pos_cur.p[0] = myStod(aStr1[is + 1]);
      pos_cur.p[1] = myStod(aStr1[is + 2]);
      is += 3;
      for (;;) {
        if (isAlphabet(aStr1[is][0])) { break; }
        CCad2D_EdgeGeo e;
        e.p0 = pos_cur;
        e.p1.p[0] = myStod(aStr1[is + 0]);
        e.p1.p[1] = myStod(aStr1[is + 1]);
        pos_cur = e.p1;
        aEdge.push_back(e);
        is += 2;
      }
    } else if (aStr1[is] == "C") { // cubic Bezier absolute coordinates
      ++is;
      for (;;) { // loop for poly Bezier
        CCad2D_EdgeGeo e;
        e.p0 = pos_cur;
        e.p1 = CVec2d(myStod(aStr1[is + 4]), myStod(aStr1[is + 5]));
        double len01 = (e.p1 - e.p0).Length();
        const CVec2d lx = (e.p1 - e.p0) / (len01 * len01);
        const CVec2d ly = CVec2d(lx.y(), -lx.x());
        e.type_edge = CCad2D_EdgeGeo::BEZIER_CUBIC;
        e.param.resize(4, 0.0);
        CVec2d p2(myStod(aStr1[is + 0]), myStod(aStr1[is + 1]));
        CVec2d p3(myStod(aStr1[is + 2]), myStod(aStr1[is + 3]));
        e.param[0] = (p2 - e.p0) * lx;
        e.param[1] = (p2 - e.p0) * ly;
        e.param[2] = (p3 - e.p0) * lx;
        e.param[3] = (p3 - e.p0) * ly;
        aEdge.push_back(e);
        pos_cur = e.p1;
        is += 6;
        if (isAlphabet(aStr1[is][0])) { break; }
      }
    } else if (aStr1[is] == "c") {
      ++is; // 'c'
      for (;;) { // loop for poly-Bezeir curve
        CCad2D_EdgeGeo e;
        e.p0 = pos_cur;
        e.p1 = CVec2d(pos_cur.x() + myStod(aStr1[is + 4]),
                            pos_cur.y() + myStod(aStr1[is + 5]));
        const CVec2d p2(pos_cur.x() + myStod(aStr1[is + 0]),
                              pos_cur.y() + myStod(aStr1[is + 1]));
        const CVec2d p3(pos_cur.x() + myStod(aStr1[is + 2]),
                              pos_cur.y() + myStod(aStr1[is + 3]));
        const double len01 = (e.p1 - e.p0).Length();
        const CVec2d lx = (e.p1 - e.p0) / (len01 * len01);
        const CVec2d ly = CVec2d(lx.y(), -lx.x());
        e.type_edge = CCad2D_EdgeGeo::BEZIER_CUBIC;
        e.param.resize(4, 0.0);
        e.param[0] = (p2 - e.p0) * lx;
        e.param[1] = (p2 - e.p0) * ly;
        e.param[2] = (p3 - e.p0) * lx;
        e.param[3] = (p3 - e.p0) * ly;
        aEdge.push_back(e);
        pos_cur = e.p1;
        is += 6;
        if (isAlphabet(aStr1[is][0])) { break; }
      }
    } else if (aStr1[is] == "l") {
      CCad2D_EdgeGeo e;
      e.p0 = pos_cur;
      e.p1.p[0] = pos_cur.x() + myStod(aStr1[is + 1]);
      e.p1.p[1] = pos_cur.y() + myStod(aStr1[is + 2]);
      pos_cur = e.p1;
      aEdge.push_back(e);
      is += 3;
      for (;;) {
        if (isAlphabet(aStr1[is][0])) { break; }
        CCad2D_EdgeGeo e0;
        e0.p0 = pos_cur;
        e0.p1.p[0] = pos_cur.x() + myStod(aStr1[is + 0]);
        e0.p1.p[1] = pos_cur.y() + myStod(aStr1[is + 1]);
        pos_cur = e0.p1;
        aEdge.push_back(e0);
        is += 2;
      }
    } else if (aStr1[is] == "L") {
      CCad2D_EdgeGeo e;
      e.p0 = pos_cur;
      e.p1.p[0] = myStod(aStr1[is + 1]);
      e.p1.p[1] = myStod(aStr1[is + 2]);
      pos_cur = e.p1;
      aEdge.push_back(e);
      is += 3;
      for (;;) {
        if (isAlphabet(aStr1[is][0])) { break; }
        CCad2D_EdgeGeo e0;
        e0.p0 = pos_cur;
        e0.p1.p[0] = myStod(aStr1[is + 0]);
        e0.p1.p[1] = myStod(aStr1[is + 1]);
        pos_cur = e0.p1;
        aEdge.push_back(e0);
        is += 2;
      }
    } else if (aStr1[is] == "S") {
      CCad2D_EdgeGeo e;
      e.p0 = pos_cur;
      e.p1 = CVec2d(myStod(aStr1[is + 3]), myStod(aStr1[is + 4]));
      double len01 = (e.p1 - e.p0).Length();
      const CVec2d lx = (e.p1 - e.p0) / (len01 * len01);
      const CVec2d ly = CVec2d(lx.y(), -lx.x());
      e.type_edge = CCad2D_EdgeGeo::BEZIER_CUBIC;
      CVec2d p2(myStod(aStr1[is + 1]), myStod(aStr1[is + 2]));
      CVec2d p3 = e.p1;
      e.param.resize(4, 0.0);
      e.param[0] = (p2 - e.p0) * lx;
      e.param[1] = (p2 - e.p0) * ly;
      e.param[2] = (p3 - e.p0) * lx;
      e.param[3] = (p3 - e.p0) * ly;
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 5;
    } else if (aStr1[is] == "v") {
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = CVec2d(pos_cur.x(), pos_cur.y() + myStod(aStr1[is + 1]));
        e.type_edge = CCad2D_EdgeGeo::LINE;
      }
      pos_cur = e.p1;
      aEdge.push_back(e);
      is += 2;
    } else if (aStr1[is] == "V") {
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = CVec2d(pos_cur.x(), myStod(aStr1[is + 1]));
        e.type_edge = CCad2D_EdgeGeo::LINE;
      }
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 2;
    } else if (aStr1[is] == "H") {
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = CVec2d(myStod(aStr1[is + 1]), pos_cur.y());
        e.type_edge = CCad2D_EdgeGeo::LINE;
      }
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 2;
    } else if (aStr1[is] == "h") {
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = CVec2d(e.p0.x() + myStod(aStr1[is + 1]), pos_cur.y());
        e.type_edge = CCad2D_EdgeGeo::LINE;
      }
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 2;
    } else if (aStr1[is] == "q") {
      const CVec2d p0 = pos_cur;
      const CVec2d p1(p0.x() + myStod(aStr1[is + 1]), p0.y() + myStod(aStr1[is + 2]));
      const CVec2d p2(p0.x() + myStod(aStr1[is + 3]), p0.y() + myStod(aStr1[is + 4]));
      const CVec2d lx = (p2 - p0) / (p2 - p0).SqLength();
      const CVec2d ly = -CVec2d(lx.y(), -lx.x());
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = p2;
        e.type_edge = CCad2D_EdgeGeo::BEZIER_QUADRATIC;
        e.param.resize(2);
        e.param[0] = (p1 - p0) * lx;
        e.param[1] = (p1 - p0) * ly;
      }
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 5;
    } else if (aStr1[is] == "Q") {
      const CVec2d p0 = pos_cur;
      const CVec2d p1(myStod(aStr1[is + 1]), myStod(aStr1[is + 2]));
      const CVec2d p2(myStod(aStr1[is + 3]), myStod(aStr1[is + 4]));
      const CVec2d lx = (p2 - p0) / (p2 - p0).SqLength();
      const CVec2d ly = -CVec2d(lx.y(), -lx.x());
      CCad2D_EdgeGeo e;
      {
        e.p0 = pos_cur;
        e.p1 = p2;
        e.type_edge = CCad2D_EdgeGeo::BEZIER_QUADRATIC;
        e.param.resize(2);
        e.param[0] = (p1 - p0) * lx;
        e.param[1] = (p1 - p0) * ly;
      }
      aEdge.push_back(e);
      pos_cur = e.p1;
      is += 5;
    } else if (aStr1[is] == "z" || aStr1[is] == "Z") {
      const CVec2d p1 = aEdge[0].p0;
      const CVec2d p0 = aEdge[aEdge.size() - 1].p1;
      double dist0 = (p0 - p1).Length();
      if (dist0 > 1.0e-9) {
        CCad2D_EdgeGeo e;
        e.p0 = p0;
        e.p1 = p1;
        e.type_edge = CCad2D_EdgeGeo::LINE;
        aEdge.push_back(e);
      }
      break;
    } else {
      std::cout << "error!--> " << aStr1[is] << std::endl;
      break;
    }
  }
}

DFM2_INLINE void LoopEdgeCad2D_SVGPolygonPoints
    (std::vector<CCad2D_EdgeGeo> &aEdge,
     std::vector<std::string> &aS) {
  const size_t np = aS.size() / 2;
  std::vector<CVec2d> aP;
  for (size_t ip = 0; ip < np; ++ip) {
    aP.emplace_back(myStod(aS[ip * 2 + 0]), myStod(aS[ip * 2 + 1]));
  }
  for (size_t ie = 0; ie < np; ++ie) {
    CCad2D_EdgeGeo e;
    e.p0 = aP[(ie + 0) % np];
    e.p1 = aP[(ie + 1) % np];
    aEdge.push_back(e);
  }
}

} // cad2
} // delfem2

// =========================================================

DFM2_INLINE void delfem2::CCad2D::Pick
 (double x0, double y0,
  double view_height)
{
  CVec2d pin(x0,y0);
  if( this->iedge_picked != -1 ){
    const CCad2D_EdgeGeo& edge = aEdge[iedge_picked];
    if( edge.type_edge == 1 ){
      assert( edge.param.size() == 4 );
      const CVec2d lx = (edge.p1 - edge.p0).Normalize();
      const CVec2d ly = CVec2d(lx.y(), -lx.x());
      const CVec2d q0 = edge.p0 + edge.param[0]*lx + edge.param[1]*ly;
      const CVec2d q1 = edge.p1 + edge.param[2]*lx + edge.param[3]*ly;
      if( Distance(pin, q0) < view_height*0.05 ){ this->ipicked_elem = 1; return; }
      if( Distance(pin, q1) < view_height*0.05 ){ this->ipicked_elem = 2; return; }
    }
  }
  this->ipicked_elem = 0;
  this->ivtx_picked = -1;
  this->iedge_picked = -1;
  this->iface_picked = -1;
  for(unsigned int ivtx=0;ivtx<aVtx.size();++ivtx){
    double x1 = aVtx[ivtx].pos.x();
    double y1 = aVtx[ivtx].pos.y();
    double dist = sqrt( (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1) );
    if( dist < view_height*0.05 ){
      this->ivtx_picked = ivtx;
      return;
    }
  }
  ////
  for(unsigned int iedge=0;iedge<aEdge.size();++iedge){
    double dist = aEdge[iedge].Distance(x0, y0);
    if( dist < view_height*0.05 ){
      this->iedge_picked = iedge;
      return;
    }
  }
  ////
  for(unsigned int iface=0;iface<aFace.size();++iface){
    bool is_inside = aFace[iface].IsInside(x0, y0, aVec2_Tessellation);
    if( is_inside ){
      this->iface_picked = iface;
      return;
    }
  }
}



DFM2_INLINE void delfem2::CCad2D::DragPicked
 (double p1x, double p1y, double p0x, double p0y)
{
  if( ivtx_picked >= 0 && ivtx_picked < (int)aVtx.size() ){
    aVtx[ivtx_picked].pos.p[0] = p1x;
    aVtx[ivtx_picked].pos.p[1] = p1y;
    Tessellation(); // tesselation for visualization
    return;
  }
  if( iedge_picked >= 0 && iedge_picked < (int)aEdge.size() ){
    if( ipicked_elem == 0 ){
      int iv0 = topo.aEdge[iedge_picked].iv0;
      int iv1 = topo.aEdge[iedge_picked].iv1;
      aVtx[iv0].pos.p[0] += p1x-p0x;
      aVtx[iv0].pos.p[1] += p1y-p0y;
      aVtx[iv1].pos.p[0] += p1x-p0x;
      aVtx[iv1].pos.p[1] += p1y-p0y;
    }
    else{
      CCad2D_EdgeGeo& edge = aEdge[iedge_picked];
      if( edge.type_edge == 1 ){
        assert( edge.param.size() == 4 );
        const CVec2d lx = (edge.p1 - edge.p0).Normalize();
        const CVec2d ly = CVec2d(lx.y(),-lx.x());
        if( ipicked_elem == 1 ){
          edge.param[0] = (CVec2d(p1x,p1y)- edge.p0 )*lx;
          edge.param[1] = (CVec2d(p1x,p1y)- edge.p0 )*ly;
        }
        else if( ipicked_elem == 2 ){
          edge.param[2] = (CVec2d(p1x,p1y)- edge.p1 )*lx;
          edge.param[3] = (CVec2d(p1x,p1y)- edge.p1 )*ly;
        }
      }
    }
    Tessellation(); // tesselation for visualization
    return;
  }
  if( iface_picked >= 0 && iface_picked < (int)aFace.size() ){
    std::vector<int> aIdV = topo.aLoop[iface_picked].GetArray_IdVertex(topo.aEdge);
    for(int iv1 : aIdV){
      aVtx[iv1].pos.p[0] += p1x-p0x;
      aVtx[iv1].pos.p[1] += p1y-p0y;
    }
    Tessellation(); // tesselation for visualization
    return;
  }
}

DFM2_INLINE bool delfem2::CCad2D::Check() const
{
  if( !this->topo.Check() ){ assert(0); return false; }
  if( (int)aVtx.size() != topo.nVertex ){ assert(0); return false; }
  if( aEdge.size() != topo.aEdge.size() ){ assert(0); return false; }
  if( aFace.size() != topo.aFace.size() ){ assert(0); return false; }
  return true;
}

DFM2_INLINE void delfem2::CCad2D::AddPolygon(const std::vector<double>& aXY)
{
  const size_t np = aXY.size()/2;
  topo.AddPolygon(np);
  for(unsigned int ip=0;ip<np;++ip){
    aVtx.emplace_back(CVec2d(aXY[ip*2+0], aXY[ip*2+1]));
  }
//  const unsigned int iedge0 = aEdge.size();
//  const unsigned int iface0 = aFace.size();
  for(unsigned int ie=0;ie<np;++ie){
    aEdge.emplace_back();
  }
  aFace.emplace_back();
  ////
  assert( this->Check() );
  Tessellation();
}

DFM2_INLINE void delfem2::CCad2D::AddFace(const std::vector<CCad2D_EdgeGeo>& aEdgeIn)
{
  if( aEdgeIn.empty() ){ return; }
  const int np = aEdgeIn.size();
  topo.AddPolygon(np);
  for(int ip=0;ip<np;++ip){
    aVtx.emplace_back(aEdgeIn[ip].p0);
  }
  for(int ie=0;ie<np;++ie){
    aEdge.push_back(aEdgeIn[ie]);
  }
  aFace.emplace_back();
  assert( this->Check() );
  Tessellation();
}

DFM2_INLINE void delfem2::CCad2D::AddVtxFace(double x0, double y0, unsigned int ifc_add)
{
  if( ifc_add >= topo.aFace.size() ){ return; }
  topo.AddVtx_Face(ifc_add);
  assert( topo.Check() );
  aVtx.emplace_back(CVec2d(x0,y0) );
  assert( this->Check() );
  Tessellation();
}

DFM2_INLINE void delfem2::CCad2D::AddVtxEdge(double x, double y, unsigned int ie_add)
{
  if( ie_add >= topo.aEdge.size() ){ return; }
  topo.AddVtx_Edge(ie_add);
  assert( topo.Check() );
  aVtx.emplace_back(CVec2d(x,y) );
  aEdge.emplace_back( );
  Tessellation();
}

DFM2_INLINE std::vector<double> delfem2::CCad2D::MinMaxXYZ() const
{
  CBoundingBox2D bb = this->BB();
  return bb.MinMaxXYZ();
}

DFM2_INLINE delfem2::CBoundingBox2D delfem2::CCad2D::BB() const
{
  CBoundingBox2D bb;
  for(const auto & ie : aEdge){
    bb += ie.BB();
  }
  return bb;
}

DFM2_INLINE void delfem2::CCad2D::GetPointsEdge
(std::vector<int>& aIdP,
 const double* pXY, int np,
 const std::vector<int>& aIE,
 double tolerance ) const
{
  aIdP.clear();
  for(int ip=0;ip<np;++ip){
    const double x = pXY[ip*2+0];
    const double y = pXY[ip*2+1];
    for(int ie0 : aIE){
      const CCad2D_EdgeGeo& eg = this->aEdge[ie0];
      const double dist = eg.Distance(x,y);
      if( dist > tolerance ){ continue; }
      aIdP.push_back(ip);
    }
  }
}

DFM2_INLINE std::vector<double> delfem2::CCad2D::XY_VtxCtrl_Face
(int iface) const
{
  std::vector<double> aXY;
  for(const auto & iie : topo.aLoop[iface].aIE){
    int ie0 = iie.first;
    bool dir = iie.second;
    if( dir ){
      const int iv = topo.aEdge[ie0].iv0;
      aXY.push_back( aVtx[iv].pos.x() );
      aXY.push_back( aVtx[iv].pos.y() );
      if( aEdge[ie0].type_edge == 1 ){
        const CCad2D_EdgeGeo& edge = aEdge[ie0];
        const CVec2d lx = (edge.p1 - edge.p0).Normalize();
        const CVec2d ly = CVec2d(lx.y(),-lx.x());
        const CVec2d q0 = edge.p0 + edge.param[0]*lx + edge.param[1]*ly;
        const CVec2d q1 = edge.p1 + edge.param[2]*lx + edge.param[3]*ly;
        aXY.push_back( q0.x() );
        aXY.push_back( q0.y() );
        aXY.push_back( q1.x() );
        aXY.push_back( q1.y() );
      }
    }
    else{
      const int iv = topo.aEdge[ie0].iv1;
      aXY.push_back( aVtx[iv].pos.x() );
      aXY.push_back( aVtx[iv].pos.y() );
    }
  }
  /*
  std::vector<int> aIdV = topo.aFace[iface].GetArray_IdVertex(topo.aEdge);
  for(unsigned int iiv=0;iiv<aIdV.size();++iiv){
    const int iv = aIdV[iiv];
    aXY.push_back( aVtx[iv].pos.x );
    aXY.push_back( aVtx[iv].pos.y );
  }
   */
  return aXY;
}

DFM2_INLINE std::vector<int> delfem2::CCad2D::Ind_Vtx_Face
 (int iface) const
{
  return topo.aLoop[iface].GetArray_IdVertex(topo.aEdge);
}

DFM2_INLINE std::vector<double> delfem2::CCad2D::XY_Vtx(int ivtx) const
{
  std::vector<double> xy;
  xy.push_back(aVtx[ivtx].pos.x());
  xy.push_back(aVtx[ivtx].pos.y());
  return xy;
}

DFM2_INLINE std::vector<int> delfem2::CCad2D::Ind_Vtx_Edge(int iedge) const
{
  std::vector<int> aRes;
  if( iedge < 0 || iedge > (int)topo.aEdge.size() ){ return aRes; }
  aRes.push_back(topo.aEdge[iedge].iv0);
  aRes.push_back(topo.aEdge[iedge].iv1);
  return aRes;
}

DFM2_INLINE std::vector<std::pair<int,bool> > delfem2::CCad2D::Ind_Edge_Face(int iface) const
{
//  std::vector<int> aIdE;
  std::vector<std::pair<int,bool> > aIdE;
  for(const auto & iie : topo.aLoop[iface].aIE){
    const int ie0 = iie.first;
    const bool dir0 = iie.second;
    aIdE.emplace_back(ie0,dir0 );
  }
  return aIdE;
}

// ----------------------------------------------


DFM2_INLINE void delfem2::CCad2D_EdgeGeo::GenMeshNDiv
 (unsigned int ndiv)
{
  assert( ndiv > 0 );
  aP.clear();
  if( type_edge == LINE ){
    for(unsigned int ip=1;ip<ndiv;++ip){
      double r2 = (double)ip/ndiv;
      CVec2d v2 = (1-r2)*p0 + r2*p1;
      aP.push_back(v2);
    }
  }
  else if( type_edge == BEZIER_QUADRATIC ){
    const CVec2d lx = (this->p1 - this->p0);
    const CVec2d ly = CVec2d(lx.y(), -lx.x());
    const CVec2d q0 = p0 + param[0]*lx + param[1]*ly;
    for(unsigned int ip=1;ip<ndiv;++ip){
      double t = (double)ip/ndiv;
      CVec2d pos = pointCurve_BezierQuadratic(t, p0, q0, p1);
      aP.push_back(pos);
    }
  }
  else if( type_edge == BEZIER_CUBIC ){
    const CVec2d lx = (this->p1 - this->p0);
    const CVec2d ly = CVec2d(lx.y(), -lx.x());
    const CVec2d q0 = p0 + param[0]*lx + param[1]*ly;
    const CVec2d q1 = p0 + param[2]*lx + param[3]*ly;
    for(unsigned int ip=1;ip<ndiv;++ip){
      double t = (double)ip/ndiv;
      CVec2d pos = pointCurve_BezierCubic(t, p0, q0, q1, p1);
      aP.push_back(pos);
    }
  }
}


DFM2_INLINE double delfem2::CCad2D_EdgeGeo::LengthNDiv
 (unsigned int ndiv) const
{
  if( type_edge == LINE ){
    return (this->p0 - this->p1).Length();
  }
  else if( this->type_edge == BEZIER_QUADRATIC ){
    const CVec2d lx = (this->p1 - this->p0);
    const CVec2d ly = CVec2d(lx.y(), -lx.x());
    const CVec2d q0 = p0 + param[0]*lx + param[1]*ly;
    std::vector<CVec2d> aP0;
    Polyline_BezierQuadratic(aP0,
                             ndiv, p0,q0,p1);
    return Length_Polygon(aP0);
  }
  else if( this->type_edge == BEZIER_CUBIC ){
    const CVec2d lx = (this->p1 - this->p0);
    const CVec2d ly = CVec2d(lx.y(), -lx.x());
    const CVec2d q0 = p0 + param[0]*lx + param[1]*ly;
    const CVec2d q1 = p0 + param[2]*lx + param[3]*ly;
    std::vector<CVec2d> aP0;
    Polyline_BezierCubic(aP0,
                         ndiv, p0,q0,q1,p1);
    return Length_Polygon(aP0);
  }
  assert(0);
  return 0;
}

// for visualization
DFM2_INLINE void delfem2::CCad2D_EdgeGeo::GenMeshLength
(double elen)
{
  if( elen <= 0 ){
    this->GenMeshNDiv(1);
    return;
  }
  double len = LengthNDiv(20);
  double ndiv = len/elen+1;
  this->GenMeshNDiv(ndiv);
}

/*
// for meshing
void CCad2D_EdgeGeo::GetInternalPoints_ElemLen
(std::vector<CVector2>& aV, double elen) const
{
  aV.clear();
  const int nadd = (int)( this->Length() / elen);
  if( nadd == 0 ) return;
  /////
  if( type_edge == 0 ){
    for(int iadd=0;iadd<nadd;++iadd){
      double r2 = (double)(iadd+1)/(nadd+1);
      CVector2 v2 = (1-r2)*p0 + r2*p1;
      aV.push_back(v2);
    }
  }
  else if( type_edge == 1 ){
    const CVector2 lx = (this->p1 - this->p0).Normalize();
    const CVector2 ly = CVector2(lx.y,-lx.x);
    const CVector2 q0 = p0 + param[0]*lx + param[1]*ly;
    const CVector2 q1 = p1 + param[2]*lx + param[3]*ly;
    for(int iadd=0;iadd<nadd;++iadd){
      double t = (double)(iadd+1)/(nadd+1);
      CVector2 pos = pointCurve_BezierCubic(t, p0, q0, q1, p1);
      aV.push_back(pos);
    }
  }
}
 */

DFM2_INLINE double delfem2::CCad2D_EdgeGeo::Distance
 (double x, double y) const
{
  const CVec2d q(x,y);
  if( type_edge == 0 ){
    const CVec2d pn = GetNearest_LineSeg_Point(q,p0,p1);
    return delfem2::Distance(pn,q);
  }
  else if( type_edge == 1 ){
    assert( param.size() == 4 );
    double min_dist = -1;
    for(size_t ie=0;ie<aP.size()+1;++ie){
      CVec2d q0 = (ie==0) ? p0 : aP[ie];
      CVec2d q1 = (ie==aP.size()-1) ? p1 : aP[ie+1];
      double dist = delfem2::Distance(q, GetNearest_LineSeg_Point(q,q0,q1));
      if( min_dist < 0 || dist < min_dist ){ min_dist = dist; }
    }
    return min_dist;
  }
  assert(0);
  return 0;
}

DFM2_INLINE void delfem2::CCad2D::MeshForVisualization
(std::vector<float>& aXYf,
 std::vector<std::vector<unsigned int> >& aaLine,
 std::vector<unsigned int>& aTri) const
{
  const std::vector<CVec2d>& aVec2 = aVec2_Tessellation;
  aXYf.reserve(aVec2.size()*2);
  for(const auto & iv : aVec2){
    aXYf.push_back(iv.x());
    aXYf.push_back(iv.y());
  }
  aaLine.resize(aEdge.size());
  for(size_t ie=0;ie<aEdge.size();++ie){
    std::vector<unsigned int>& aLine = aaLine[ie];
    aLine.clear();
    aLine.push_back(topo.aEdge[ie].iv0);
    for(size_t ip=0;ip<aEdge[ie].aP.size();++ip){
      aLine.push_back(ip+aEdge[ie].ip0);
    }
    aLine.push_back(topo.aEdge[ie].iv1);
  }
  for(const auto & fc : aFace){
    aTri.insert(aTri.end(),fc.aTri.begin(),fc.aTri.end());
  }
}

DFM2_INLINE double delfem2::CCad2D_EdgeGeo::LengthMesh() const
{
  double len0 = 0.0;
  for(size_t ie=0;ie<aP.size()+1;++ie){
    const CVec2d q0 = (ie==0) ? p0 : aP[ie-1];
    const CVec2d q1 = (ie==aP.size()) ? p1 : aP[ie];
    double dist0 = delfem2::Distance(q0, q1);
    len0 += dist0;
  }
  return len0;
}


DFM2_INLINE double delfem2::AreaLoop
(const std::vector<CCad2D_EdgeGeo>& aEdge)
{
  double a0 = 0;
  CVec2d qo(0,0);
  for(const auto & ie : aEdge){
    const std::vector<CVec2d>& aP = ie.aP;
    const unsigned int nseg = aP.size()+1;
    for(unsigned int iseg=0;iseg<nseg;++iseg){
      const CVec2d q0 = (iseg==0) ? ie.p0 : aP[iseg-1];
      const CVec2d q1 = (iseg==nseg-1) ?ie.p1 : aP[iseg];
      a0 += Area_Tri(qo, q0, q1);
    }
  }
  return a0;
}


DFM2_INLINE std::vector<delfem2::CCad2D_EdgeGeo> delfem2::InvertLoop
(const std::vector<CCad2D_EdgeGeo>& aEdge)
{
  const unsigned int ne = aEdge.size();
  std::vector<CCad2D_EdgeGeo> aEdgeOut(ne);
  for(unsigned int ie=0;ie<ne;++ie){
    const CCad2D_EdgeGeo& ei = aEdge[ie];
    CCad2D_EdgeGeo& eo = aEdgeOut[ne-ie-1];
    eo.p1 = ei.p0;
    eo.p0 = ei.p1;
    eo.type_edge = ei.type_edge;
    eo.param.resize(ei.param.size());
    for(unsigned int ip=0;ip<ei.param.size()/2;++ip){
      int jp = ei.param.size()/2-1-ip;
      eo.param[ip*2+0] = 1-ei.param[jp*2+0];
      eo.param[ip*2+1] = -ei.param[jp*2+1];
    }
  }
  return aEdgeOut;
}


DFM2_INLINE std::vector<delfem2::CCad2D_EdgeGeo> delfem2::RemoveEdgeWithZeroLength
(const std::vector<CCad2D_EdgeGeo>& aEdge)
{
  const unsigned int ne = aEdge.size();
  std::vector<CCad2D_EdgeGeo> aEdgeOut;
  aEdgeOut.reserve(ne);
  for(unsigned int ie=0;ie<ne;++ie){
    if( aEdge[ie].LengthMesh() < 1.0e-10 ) continue;
    aEdgeOut.push_back(aEdge[ie]);
  }
  return aEdgeOut;
}


// ---------------------------------------------------------------

DFM2_INLINE void delfem2::CCad2D::Tessellation()
{
  for(size_t ie=0;ie<topo.aEdge.size();++ie){
    const int iv0 = topo.aEdge[ie].iv0;
    const int iv1 = topo.aEdge[ie].iv1;
    aEdge[ie].p0 = aVtx[iv0].pos;
    aEdge[ie].p1 = aVtx[iv1].pos;
    aEdge[ie].GenMeshNDiv(20);
  }
  std::vector<CVec2d>& aVec2 = this->aVec2_Tessellation;
  aVec2.clear();
  for(auto & iv : aVtx){
    aVec2.push_back( iv.pos );
  }
  for(auto & ie : aEdge){
    ie.ip0 = aVec2.size();
    for(const auto & ip : ie.aP){
      aVec2.push_back(ip);
    }
  }
  for(size_t ifc=0;ifc<topo.aFace.size();++ifc){
    std::vector<CDynTri> aETri;
    cad2::GenMeshCadFace(aVec2, aETri,
        aFace[ifc],ifc,
        topo,
        aVtx,aEdge);
    std::vector<unsigned int>& aTri = aFace[ifc].aTri;
    aTri.clear();
    const int ntri = (int)aETri.size();
    aTri.resize(ntri*3);
    for(int itri=0;itri<ntri;itri++){
      aTri[itri*3+0] = aETri[itri].v[0];
      aTri[itri*3+1] = aETri[itri].v[1];
      aTri[itri*3+2] = aETri[itri].v[2];
    }
  }
}


// above: CCad2
// ======================================================================
// below: CMesher_Cad2

DFM2_INLINE void delfem2::CMesher_Cad2D::Meshing
(CMeshDynTri2D& dmsh,
 const CCad2D& cad)
{
  std::vector<CCad2D_EdgeGeo> aEdgeGeo = cad.aEdge;
  for(size_t ie=0;ie<aEdgeGeo.size();++ie){
    const int iv0 = cad.topo.aEdge[ie].iv0;
    const int iv1 = cad.topo.aEdge[ie].iv1;
    aEdgeGeo[ie].p0 = cad.aVtx[iv0].pos;
    aEdgeGeo[ie].p1 = cad.aVtx[iv1].pos;
    const auto itr = this->mapIdEd_NDiv.find(ie);
    if( itr == this->mapIdEd_NDiv.end() ){
      aEdgeGeo[ie].GenMeshLength(edge_length);
    }
    else{
      aEdgeGeo[ie].GenMeshNDiv(itr->second);
    }
  }
  //
  aFlgPnt.clear();
  dmsh.Clear();
  for(size_t iv=0;iv<cad.aVtx.size();++iv){
    dmsh.aVec2.push_back( cad.aVtx[iv].pos );
    aFlgPnt.push_back(iv);
  }
  for(size_t ie=0;ie<aEdgeGeo.size();++ie){
    aEdgeGeo[ie].ip0 = dmsh.aVec2.size();
    for(size_t ip=0;ip<aEdgeGeo[ie].aP.size();++ip){
      dmsh.aVec2.push_back(aEdgeGeo[ie].aP[ip]);
      aFlgPnt.push_back(cad.aVtx.size()+ie);
    }
  }
  //
  aFlgTri.clear();
  for(size_t ifc=0;ifc<cad.aFace.size();++ifc){ // add face to dmsh
    std::vector<CDynTri> aDTri;
    cad2::GenMeshCadFace(dmsh.aVec2, aDTri,
        cad.aFace[ifc], ifc,
        cad.topo,cad.aVtx,aEdgeGeo);
    const unsigned int ntri0 = dmsh.aETri.size();
    for(auto & tri : aDTri){
      if( tri.s2[0]!=UINT_MAX ){ tri.s2[0] += ntri0; }
      if( tri.s2[1]!=UINT_MAX ){ tri.s2[1] += ntri0; }
      if( tri.s2[2]!=UINT_MAX ){ tri.s2[2] += ntri0; }
      dmsh.aETri.push_back(tri);
    }
    aFlgTri.resize(dmsh.aETri.size(),ifc);
  }
  { // make EPo
    dmsh.aEPo.resize(dmsh.aVec2.size());
    for(size_t it=0;it<dmsh.aETri.size();++it){
      for(unsigned int inotri=0;inotri<3;++inotri){
        dmsh.aEPo[ dmsh.aETri[it].v[inotri] ].e = it;
        dmsh.aEPo[ dmsh.aETri[it].v[inotri] ].d = inotri;
      }
    }
  }
  if( edge_length > 1.0e-10 ){
    CInputTriangulation_Uniform param(1.0);
    MeshingInside(dmsh.aEPo,dmsh.aETri,dmsh.aVec2,
                  aFlgPnt,aFlgTri,
                  dmsh.aVec2.size(), cad.aVtx.size()+cad.aEdge.size(),
                  edge_length, param);
  }
  nvtx = cad.aVtx.size();
  nedge = cad.aEdge.size();
  nface = cad.aFace.size();
}


DFM2_INLINE std::vector<unsigned int> delfem2::CMesher_Cad2D::IndPoint_IndEdge
(const int iedge,
 bool is_end_point,
 const CCad2D& cad2d)
{
  std::vector<unsigned int> res;
  if( iedge >= cad2d.aEdge.size() ){ return res; }
  //    std::cout << nvtx << " " << nedge << " " << nface << std::endl;
  std::vector<int> aflg(nvtx+nedge+nface,0);
  {
    aflg[nvtx+iedge] = 1;
  }
  std::vector<int> aIP_E = cad2d.Ind_Vtx_Edge(iedge);
  if( is_end_point ){ res.push_back(aIP_E[0]); }
  for(size_t ip=0;ip<this->aFlgPnt.size();++ip){
    int iflg = aFlgPnt[ip]; assert(iflg<int(nvtx+nedge+nface));
    if( iflg >= (int)(nvtx+nedge) ){ break; }
    if( aflg[iflg] == 1 ){ res.push_back(ip); }
  }
  if( is_end_point ){ res.push_back(aIP_E[1]); }
  return res;
}

DFM2_INLINE std::vector<unsigned int> delfem2::CMesher_Cad2D::IndPoint_IndEdgeArray
(const std::vector<int>& aIndEd,
 const CCad2D& cad2d)
{
  //    std::cout << nvtx << " " << nedge << " " << nface << std::endl;
  std::vector<int> aflg(nvtx+nedge+nface,0);
  {
    for(unsigned int iedge : aIndEd){
      assert(iedge<nedge);
      aflg[iedge+nvtx] = 1;
      {
        const std::vector<int>& aIV = cad2d.Ind_Vtx_Edge(iedge);
        for(int iv0 : aIV){
          aflg[iv0] = 1;
        }
      }
    }
  }
  std::vector<unsigned int> res;
  for(size_t ip=0;ip<this->aFlgPnt.size();++ip){
    int iflg = aFlgPnt[ip]; assert(iflg<int(nvtx+nedge+nface));
    if( aflg[iflg] == 1 ){ res.push_back(ip); }
  }
  return res;
}


DFM2_INLINE std::vector<int> delfem2::CMesher_Cad2D::IndPoint_IndFaceArray
(const std::vector<int>& aIndFc,
 const CCad2D& cad2d)
{
  std::vector<int> aflg(nvtx+nedge+nface,0);
  {
    for(unsigned int iface : aIndFc){
      assert(iface<nface);
      aflg[nvtx+nedge+iface] = 1;
      {
        const std::vector<std::pair<int,bool> >& aIE = cad2d.Ind_Edge_Face(iface);
        for(const auto & iie : aIE){
          const int ie0 = iie.first;
          aflg[nvtx+ie0] = 1;
        }
      }
      {
        const std::vector<int>& aIV = cad2d.Ind_Vtx_Face(iface);
        for(int iv0 : aIV){
          aflg[iv0] = 1;
        }
      }
    }
  }
  std::vector<int> res;
  for(size_t ip=0;ip<this->aFlgPnt.size();++ip){
    int iflg = aFlgPnt[ip]; assert(iflg<int(nvtx+nedge+nface));
    if( aflg[iflg] == 1 ){ res.push_back(ip); }
  }
  return res;
}


DFM2_INLINE bool delfem2::WriteCAD_DXF
(const std::string& file_name,
 const CCad2D& cad,
 double scale)
{
  FILE *fp;
  if( (fp = ::fopen(file_name.c_str(),"w"))== nullptr ){
    fclose(fp);
    assert(0);
    return false;
  }
  CBoundingBox2D bb;
  {  // Get Bounding Box of this Object
    for(const auto & edge : cad.aEdge){
      bb += edge.BB();
    }
  }
  // header section
  fprintf(fp, "  0\nSECTION\n");
  fprintf(fp, "  2\nHEADER\n");
  fprintf(fp, "  9\n$ACADVER\n  1\nAC1009\n");
  fprintf(fp, "  9\n$EXTMIN\n  10\n%lf\n  20\n%lf\n",bb.x_min*scale,bb.y_min*scale);
  fprintf(fp, "  9\n$EXTMAX\n  10\n%lf\n  20\n%lf\n",bb.x_max*scale,bb.y_max*scale);
  fprintf(fp, "  0\nENDSEC\n");
  // table section
  fprintf(fp, "  0\nSECTION\n");
  fprintf(fp, "  2\nTABLES\n");
  fprintf(fp, "  0\nENDSEC\n");
  // block section
  fprintf(fp, "  0\nSECTION\n");
  fprintf(fp, "  2\nBLOCKS\n");
  fprintf(fp, "  0\nENDSEC\n");
  // entity section
  fprintf(fp,"  0\nSECTION\n");
  fprintf(fp,"  2\nENTITIES\n");
  for(size_t ifc=0;ifc<cad.aFace.size();++ifc){
    const std::vector<int>& aIL = cad.topo.aFace[ifc].aIL;
    for(int il0 : aIL){
      const std::vector< std::pair<int,bool> >& aIE = cad.topo.aLoop[il0].aIE;
      for(const auto & iie : aIE){
        unsigned int ie0 = iie.first;
//        bool dir0 = aIE[iie].second;
        unsigned int id_vs = cad.topo.aEdge[ie0].iv0;
        unsigned int id_ve = cad.topo.aEdge[ie0].iv1;
        const CVec2d& ps = cad.aVtx[id_vs].pos;
        const CVec2d& pe = cad.aVtx[id_ve].pos;
        if( cad.aEdge[ie0].type_edge == 0 ){
          fprintf(fp,"  0\nLINE\n  8\n%d\n  6\nCONTINUOUS\n  62\n7\n",il0);
          fprintf(fp,"  10\n%lf\n",ps.x()*scale);
          fprintf(fp,"  20\n%lf\n",ps.y()*scale);
          fprintf(fp,"  11\n%lf\n",pe.x()*scale);
          fprintf(fp,"  21\n%lf\n",pe.y()*scale);
        }
        /*
        else if( this->GetEdgeCurveType(id_e) == 1 ){ // Arc
          const CEdge2D& edge = this->GetEdge(id_e);
          CVector2D pc;  double r;
          edge.GetCenterRadius(pc,r);
          double d1, d2;
          {
            CVector2D vs = ps - pc;
            CVector2D ve = pe - pc;
            double ds = atan2(vs.y,vs.x); ds = ds * 180.0 / 3.14159265; if( ds < 0.0 ) ds += 360;
            double de = atan2(ve.y,ve.x); de = de * 180.0 / 3.14159265; if( de < 0.0 ) de += 360;
            if( edge.is_left_side ){ d1 = de; d2 = ds; }
            else{                    d1 = ds; d2 = de; }
          }
          fprintf(fp,"  0\nARC\n  8\n%d\n  6\nCONTINUOUS\n  62\n7\n  100\nAcDbCircle\n",id_l);
          fprintf(fp,"  10\n%lf\n",pc.x*scale);  // x coord
          fprintf(fp,"  20\n%lf\n",pc.y*scale);  // y coord
          fprintf(fp,"  40\n%lf\n",r*scale);  // radius
          fprintf(fp,"  100\nAcDbArc\n");
          fprintf(fp,"  50\n%lf\n",d1);
          fprintf(fp,"  51\n%lf\n",d2);
        }
        else if( this->GetEdgeCurveType(id_e) == 2 ){ // polyline
          const CEdge2D& edge = this->GetEdge(id_e);
          fprintf(fp,"  0\nPOLYLINE\n  8\n%d\n  6\nCONTINUOUS\n",id_l);
          fprintf(fp,"  10\n0.0\n");
          fprintf(fp,"  20\n0.0\n");
          fprintf(fp,"  30\n0.0\n");
          fprintf(fp,"  70\n8\n");
          fprintf(fp,"  66\n1\n");
          ////
          const std::vector<double>& axys = edge.aRelCoMesh;
          assert( axys.size() % 2 == 0 );
          const unsigned int nno = axys.size()/2;
          const Com::CVector2D& po_s = this->GetVertexCoord( edge.id_v_s );
          const Com::CVector2D& po_e = this->GetVertexCoord( edge.id_v_e );
          Com::CVector2D v0 = po_e-po_s;
          Com::CVector2D v1(-v0.y,v0.x);
          fprintf(fp,"  0\nVERTEX\n 8\n0\n 10\n%lf\n 20\n%lf\n 30\n%lf\n", po_s.x*scale, po_s.y*scale, 0.0);
          for(unsigned int ino=0;ino<nno;ino++){
            const Com::CVector2D& p = po_s + v0*axys[ino*2+0] + v1*axys[ino*2+1];
            fprintf(fp,"  0\nVERTEX\n 8\n0\n 10\n%lf\n 20\n%lf\n 30\n%lf\n", p.x*scale, p.y*scale, 0.0);
          }
          fprintf(fp,"  0\nVERTEX\n 8\n0\n 10\n%lf\n 20\n%lf\n 30\n%lf\n", po_e.x*scale, po_e.y*scale, 0.0);
          fprintf(fp,"  0\nSEQEND\n");
        }
         */
      }
    }
  }
  fprintf(fp, "  0\nENDSEC\n  0\nEOF\n");
  fclose(fp);
  return true;
}


DFM2_INLINE void delfem2::ReadSVG_LoopEdgeCCad2D
(std::vector< std::vector<CCad2D_EdgeGeo>> & aaEdge,
 const std::string& fname)
{
  aaEdge.clear();
  std::vector<char> aC;
  if( !GetFileContents(aC, fname) ){ return; }
  
  // ----
  /*
  std::cout << "svg file content: ";
  for(unsigned int ic=0;ic<aC.size();++ic){ std::cout << aC[ic]; }
  std::cout << std::endl;
  */
  // ----
  
  std::vector< std::string > aStrTagContent;
  XML_SeparateTagContent(aStrTagContent,
                         aC);
  
  { // get path
    for(auto & sTagContent : aStrTagContent){
//      std::cout << "tagcontent: " << sTagContent << std::endl;
      std::string str_path;
      if( sTagContent.compare(0,4,"path") == 0 || // adobe illustrator
          sTagContent.compare(0,5,"path\r") == 0){ // inkscape
        str_path = std::string(sTagContent.begin()+5,sTagContent.end());
      }
      if( str_path == "" ){ continue; }
      // remove new line codes as inkscape svg has new line in side path tag
      // don't remove spaces here as inkscape svg uses space for the delimiter
      str_path = Remove(str_path, "\n\r");
//      std::cout << "str_path: " << str_path << std::endl;
      std::map< std::string, std::string > mapAttr;
      ParseAttributes(mapAttr,
                      str_path);
      std::string str_path_d = mapAttr["d"];
//      std::cout << "str_path_d: " << str_path_d << std::endl;
      if( str_path_d.empty() ){ continue; }
//      std::cout << "str_path_d: " << str_path_d << std::endl;
      std::vector<std::string> aStr1 = cad2::SVG_Split_Path_d(str_path_d);
      /*
      for(unsigned int is=0;is<aStr1.size();++is){
        std::cout << is << " " << aStr1[is] << std::endl;
      }
       */
      std::vector<CCad2D_EdgeGeo> aEdge;
      cad2::LoopEdgeCad2D_SVGPathD(aEdge,
          aStr1);
      /*
      for(int ie=0;ie<aEdge.size();++ie){
        std::cout << ie << " " << aEdge.size() << " " << aEdge[ie].param.size() << std::endl;
      }
       */
      aaEdge.push_back(aEdge);
    }
  }
  
  { // get polygon
    for(auto & sTagContent : aStrTagContent){
      std::string str_polygon;
      if( sTagContent.compare(0,8,"polygon ") == 0 ){
        str_polygon = std::string(sTagContent.begin()+8,sTagContent.end());
      }
//    std::cout << "str_polygon: " << str_polygon << std::endl;
      if( str_polygon == "" ){ continue; }
      std::map< std::string, std::string > mapAttr;
      ParseAttributes(mapAttr,
                      str_polygon);
      std::string str_polygon_points = mapAttr["points"];
      std::vector<std::string> aS = Split(str_polygon_points, "  ,");
      /*
      for(unsigned int is=0;is<aS.size();++is){
        std::cout << is << " " << aS[is] << std::endl;
      }
       */
      std::vector<CCad2D_EdgeGeo> aEdge;
      cad2::LoopEdgeCad2D_SVGPolygonPoints(aEdge,
          aS);
      aaEdge.push_back(aEdge);
    }
  }
}


DFM2_INLINE void delfem2::ReadSVG_Cad2D
(delfem2::CCad2D& cad,
const std::string& fpath,
double scale)
{
  std::vector< std::vector<delfem2::CCad2D_EdgeGeo> > aaEdge;
  ReadSVG_LoopEdgeCCad2D(aaEdge,
                         fpath);
  cad.Clear();
  for(unsigned int iae=0;iae<aaEdge.size();++iae){
    std::vector<delfem2::CCad2D_EdgeGeo> aEdge = aaEdge[iae];
    Transform_LoopEdgeCad2D(aEdge,false,true,scale,scale);
    if( AreaLoop(aEdge) < 0 ){ aEdge = InvertLoop(aEdge); }
    aEdge = RemoveEdgeWithZeroLength(aEdge);
    for(auto & ie : aEdge){ ie.GenMeshLength(-1); }
    cad.AddFace(aEdge);
  }
}


DFM2_INLINE void delfem2::Transform_LoopEdgeCad2D
(std::vector<CCad2D_EdgeGeo>& aEdge,
 bool is_flip_holizontal,
 bool is_flip_vertical,
 double scale_x,
 double scale_y)
{
  double A[4] = {scale_x,0,0,scale_y};
  if( is_flip_holizontal ){ A[0] *= -1; }
  if( is_flip_vertical ){ A[3] *= -1; }
  bool is_det_inv = (is_flip_holizontal != is_flip_vertical );
  for(auto & ie : aEdge){
    ie.p0 = Mat2Vec(A, ie.p0);
    ie.p1 = Mat2Vec(A, ie.p1);
    if( ie.type_edge == 1 && is_det_inv ){
      assert( ie.param.size() == 4);
      ie.param[1] *= -1;
      ie.param[3] *= -1;
    }
  }
}
