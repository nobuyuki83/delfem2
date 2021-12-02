/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#include "delfem2/cad2.h"

#include <deque>
#include <climits>

#include "delfem2/geo_polyline2.h"
#include "delfem2/geo_edge.h"
#include "delfem2/geo_bezier_quadratic.h"

// ------------------------------------------

namespace delfem2::cad2 {

DFM2_INLINE delfem2::CBoundingBox2<double> BB_LoopEdgeCad2D(
  const std::vector<CCad2D_EdgeGeo> &aEdge) {
  CBoundingBox2<double> bb;
  for (const auto &ie: aEdge) {
    CBoundingBox2<double> bb0 = ie.BB();
    bb += bb0;
  }
  return bb;
}

DFM2_INLINE void GetBound(
  double bound_2d[4],
  unsigned int ifc0,
  const CadTopo &topo,
  const std::vector<CCad2D_EdgeGeo> &aEdgeGeo) {
  assert(ifc0 < topo.faces.size());
  const int il0 = topo.faces[ifc0].aIL[0];
  const std::vector<std::pair<unsigned int, bool> > &aIE = topo.loops[il0].aIE;
  {
    int ie0 = aIE[0].first;
    const CVec2d &p0 = aEdgeGeo[ie0].p0;
    bound_2d[0] = p0.x;
    bound_2d[1] = p0.x;
    bound_2d[2] = p0.y;
    bound_2d[3] = p0.y;
  }
  for (const auto &iie: aIE) {
    const unsigned int ie0 = iie.first;
    {
      const CVec2d &p0 = aEdgeGeo[ie0].p0;
      if (p0.x < bound_2d[0]) { bound_2d[0] = p0.x; }
      if (p0.x > bound_2d[1]) { bound_2d[1] = p0.x; }
      if (p0.y < bound_2d[2]) { bound_2d[2] = p0.y; }
      if (p0.y > bound_2d[3]) { bound_2d[3] = p0.y; }
    }
    /*
    for (const auto &p0 : aEdgeGeo[ie0].aP) {
      if (p0.x < bound_2d[0]) { bound_2d[0] = p0.x; }
      if (p0.x > bound_2d[1]) { bound_2d[1] = p0.x; }
      if (p0.y < bound_2d[2]) { bound_2d[2] = p0.y; }
      if (p0.y > bound_2d[3]) { bound_2d[3] = p0.y; }
    }
     */
  }
}

} // delfem2::cad2

// =========================================================

DFM2_INLINE bool delfem2::CCad2D::Check() const {
  this->topo.Assert();
  if (aVtx.size() != topo.num_vertex) {
    assert(0);
    return false;
  }
  if (aEdge.size() != topo.edges.size()) {
    assert(0);
    return false;
  }
  return true;
}

DFM2_INLINE void delfem2::CCad2D::AddPolygon(
  const std::vector<double> &aXY) {
  const size_t np = aXY.size() / 2;
  topo.AddPolygon(static_cast<unsigned int>(np));
  for (unsigned int ip = 0; ip < np; ++ip) {
    aVtx.emplace_back(CVec2d(aXY[ip * 2 + 0], aXY[ip * 2 + 1]));
  }
  for (unsigned int ie = 0; ie < np; ++ie) {
    aEdge.emplace_back();
  }
  assert(this->Check());
  this->CopyVertexPositionsToEdges();
}

DFM2_INLINE void delfem2::CCad2D::AddFace(
  const std::vector<CCad2D_EdgeGeo> &aEdgeIn) {
  if (aEdgeIn.empty()) { return; }
  const size_t np = aEdgeIn.size();
  topo.AddPolygon(static_cast<unsigned int>(np));
  for (unsigned int ip = 0; ip < np; ++ip) {
    aVtx.emplace_back(aEdgeIn[ip].p0);
  }
  for (unsigned int ie = 0; ie < np; ++ie) {
    aEdge.push_back(aEdgeIn[ie]);
  }
  assert(this->Check());
  this->CopyVertexPositionsToEdges();
}

DFM2_INLINE void delfem2::CCad2D::AddVtxFace(
  double x0, double y0, unsigned int ifc_add) {
  if (ifc_add >= topo.faces.size()) { return; }
  topo.AddVtx_Face(ifc_add);
  topo.Assert();
  aVtx.emplace_back(CVec2d(x0, y0));
  assert(this->Check());
  this->CopyVertexPositionsToEdges();
}

DFM2_INLINE void delfem2::CCad2D::AddVtxEdge(
  double x, double y, unsigned int ie_add) {
  if (ie_add >= topo.edges.size()) { return; }
  topo.AddVtx_Edge(ie_add);
  topo.Assert();
  aVtx.emplace_back(CVec2d(x, y));
  aEdge.emplace_back();
  this->CopyVertexPositionsToEdges();
}

DFM2_INLINE std::vector<double> delfem2::CCad2D::MinMaxXYZ() const {
  CBoundingBox2<double> bb = this->BB();
  return bb.MinMaxXYZ();
}

DFM2_INLINE delfem2::CBoundingBox2<double> delfem2::CCad2D::BB() const {
  CBoundingBox2<double> bb;
  for (const auto &ie: aEdge) {
    bb += ie.BB();
  }
  return bb;
}

// ----------------------------------------------

DFM2_INLINE std::vector<double> delfem2::CCad2D_EdgeGeo::GenMesh(
  unsigned int ndiv) const {
  std::vector<double> xyw;
  assert(ndiv > 0);
  if (type_edge == LINE) {
    for (unsigned int ip = 1; ip < ndiv; ++ip) {
      double r2 = static_cast<double>(ip) / ndiv;
      CVec2d v2 = (1 - r2) * p0 + r2 * p1;
      xyw.push_back(v2.x);
      xyw.push_back(v2.y);
      xyw.push_back(r2);
    }
  } else if (type_edge == BEZIER_QUADRATIC) {
    const CVec2d q0 = p0 + CVec2d(param[0], param[1]);
    for (unsigned int ip = 1; ip < ndiv; ++ip) {
      double t = (double) ip / ndiv;
      CVec2d pos = PointOnQuadraticBezierCurve(t, p0, q0, p1);
      xyw.push_back(pos.x);
      xyw.push_back(pos.y);
      xyw.push_back(t);
    }
  } else if (type_edge == BEZIER_CUBIC) {
    const CVec2d q0 = p0 + CVec2d(param[0], param[1]);
    const CVec2d q1 = p1 + CVec2d(param[2], param[3]);
    for (unsigned int ip = 1; ip < ndiv; ++ip) {
      double t = (double) ip / ndiv;
      CVec2d pos = PointOnCubicBezierCurve(t, p0, q0, q1, p1);
      xyw.push_back(pos.x);
      xyw.push_back(pos.y);
      xyw.push_back(t);
    }
  }
  return xyw;
}

DFM2_INLINE double delfem2::CCad2D_EdgeGeo::Distance(
  double x, double y) const {
  const CVec2d q(x, y);
  if (type_edge == LINE) {
    const CVec2d pn = Nearest_Edge_Point(q, p0, p1);
    return delfem2::Distance(pn, q);
  } else if (type_edge == BEZIER_QUADRATIC) {
    assert(param.size() == 2);
    CVec2d p2 = CVec2d(param[0], param[1]) + p0;
    double t = Nearest_QuadraticBezierCurve(q, p0, p2, p1, 5, 3);
    CVec2d q0 = PointOnQuadraticBezierCurve(t, p0, p2, p1);
    return (q0 - q).norm();
  } else if (type_edge == BEZIER_CUBIC) {
    assert(param.size() == 4);
    const CVec2d r0 = CVec2d(param[0], param[1]) + p0;
    const CVec2d r1 = CVec2d(param[2], param[3]) + p1;
    double t = Nearest_CubicBezierCurve(q, p0, r0, r1, p1, 5, 3);
    CVec2d q0 = PointOnCubicBezierCurve(t, p0, r0, r1, p1);
    // std::cout << t << " " << q0 << " " << (q0-q).norm() << std::endl;
    return (q0 - q).norm();
  }
  assert(0);
  return 0;
}

DFM2_INLINE double delfem2::CCad2D_EdgeGeo::ArcLength() const {
  if (this->type_edge == BEZIER_QUADRATIC) {
    return Length_QuadraticBezierCurve_Analytic<CVec2d>(
      p0,
      p0 + CVec2d(param[0], param[1]),
      p1);
  } else if (this->type_edge == BEZIER_CUBIC) {
    return Length_CubicBezierCurve_Quadrature<CVec2d>(
      p0,
      p0 + CVec2d(param[0], param[1]),
      p1 + CVec2d(param[2], param[3]),
      p1, 3);
  }
  assert(this->type_edge == LINE);
  return delfem2::Distance(p0, p1);
}

double delfem2::CCad2D_EdgeGeo::AreaEnclosingOrigin() const {
  if (this->type_edge == BEZIER_QUADRATIC) {
    return delfem2::Area_QuadraticBezierCurve2<CVec2d>(
      p0,
      p0 + CVec2d(param[0], param[1]),
      p1);
  }
  else if (this->type_edge == BEZIER_CUBIC ) {
    return delfem2::Area_CubicBezierCurve2<CVec2d>(
      p0,
      p0 + CVec2d(param[0], param[1]),
      p1 + CVec2d(param[2], param[3]),
      p1);
  }
  assert(this->type_edge == LINE);
  return delfem2::Area_Tri2(CVec2d(0,0), p0, p1);
}

DFM2_INLINE double delfem2::AreaLoop(
  const std::vector<CCad2D_EdgeGeo> &aEdge) {
  double a0 = 0;
  for (const auto &ie: aEdge) {
    a0 += ie.AreaEnclosingOrigin();
  }
  return a0;
}

DFM2_INLINE std::vector<delfem2::CCad2D_EdgeGeo> delfem2::InvertLoop(
  const std::vector<CCad2D_EdgeGeo> &aEdge) {
  const size_t ne = aEdge.size();
  std::vector<CCad2D_EdgeGeo> aEdgeOut(ne);
  for (unsigned int ie = 0; ie < ne; ++ie) {
    const CCad2D_EdgeGeo &ei = aEdge[ie];
    CCad2D_EdgeGeo &eo = aEdgeOut[ne - ie - 1];
    eo.p1 = ei.p0;
    eo.p0 = ei.p1;
    eo.type_edge = ei.type_edge;
    eo.param.resize(ei.param.size());
    if (ei.type_edge == CCad2D_EdgeGeo::BEZIER_CUBIC) {
      CVec2d pA = ei.p0 + CVec2d(ei.param[0], ei.param[1]);
      CVec2d pB = ei.p1 + CVec2d(ei.param[2], ei.param[3]);
      eo.param[0] = (pB - eo.p0).x;
      eo.param[1] = (pB - eo.p0).y;
      eo.param[2] = (pA - eo.p1).x;
      eo.param[3] = (pA - eo.p1).y;
    }
  }
  return aEdgeOut;
}

DFM2_INLINE std::vector<delfem2::CCad2D_EdgeGeo>
delfem2::RemoveEdgeWithZeroLength(
  const std::vector<CCad2D_EdgeGeo> &aEdge) {
  const size_t ne = aEdge.size();
  std::vector<CCad2D_EdgeGeo> aEdgeOut;
  aEdgeOut.reserve(ne);
  for (unsigned int ie = 0; ie < ne; ++ie) {
    if (aEdge[ie].ArcLength() < 1.0e-10) continue;
    aEdgeOut.push_back(aEdge[ie]);
  }
  return aEdgeOut;
}


// ---------------------------------------------------------------

/*
DFM2_INLINE void delfem2::CCad2D::Tessellation()
{
  std::vector< std::vector< std::pair<unsigned int, double> > > edge_point(aEdge.size());
  std::vector< std::vector<double> > edge_xyw(aEdge.size());
  for(size_t ie=0;ie<topo.edges.size();++ie){
    const int iv0 = topo.edges[ie].iv0;
    const int iv1 = topo.edges[ie].iv1;
    aEdge[ie].p0 = aVtx[iv0].pos;
    aEdge[ie].p1 = aVtx[iv1].pos;
    edge_xyw[ie] = aEdge[ie].GenMesh(20);
    for(size_t ip=0;ip<edge_xyw[ie].size()/3;++ip){
      edge_point[ie].push_back(std::make_pair(UINT_MAX,edge_xyw[ie][ip*3+2]));
    }
  }
  std::vector<CVec2d>& aVec2 = this->aVec2_Tessellation;
  aVec2.clear();
  for(auto & iv : aVtx){
    aVec2.push_back( iv.pos );
  }
  for(unsigned int ie=0;ie<aEdge.size();++ie){
    for(unsigned int ip=0;ip<edge_xyw[ie].size()/3;++ip){
      edge_point[ie][ip].first = aVec2.size();
      aVec2.push_back({edge_xyw[ie][ip*3+0],edge_xyw[ie][ip*3+1]});
    }
  }
  for(unsigned int ifc=0;ifc<topo.faces.size();++ifc){
    std::vector<CDynTri> aETri;
    cad2::GenMeshCadFace(aVec2, aETri,
        aFace[ifc],ifc,
        topo,
        aVtx, aEdge, edge_point);
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
 */


// ================================


DFM2_INLINE void delfem2::Cad2_Ui::Pick(
  double x0,
  double y0,
  double view_height,
  const CCad2D &cad) {
  CVec2d pin(x0, y0);
  if (this->iedge_picked != UINT_MAX) {  // try to pick handles
    const CCad2D_EdgeGeo &edge = cad.aEdge[iedge_picked];
    if (edge.type_edge == CCad2D_EdgeGeo::BEZIER_CUBIC) {
      const CVec2d q0 = edge.p0 + CVec2d(edge.param[0], edge.param[1]);
      const CVec2d q1 = edge.p1 + CVec2d(edge.param[2], edge.param[3]);
      if (Distance(pin, q0) < view_height * 0.05) {
        this->ipicked_elem = 1;
        return;
      }
      if (Distance(pin, q1) < view_height * 0.05) {
        this->ipicked_elem = 2;
        return;
      }
    }
  }
  this->ipicked_elem = 0;
  this->ivtx_picked = -1;
  this->iedge_picked = -1;
  this->iface_picked = -1;
  for (unsigned int ivtx = 0; ivtx < cad.aVtx.size(); ++ivtx) {
    double x1 = cad.aVtx[ivtx].pos.x;
    double y1 = cad.aVtx[ivtx].pos.y;
    double dist = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
    if (dist < view_height * 0.05) {
      this->ivtx_picked = ivtx;
      picked_pos = {
        static_cast<float>(cad.aVtx[ivtx_picked].pos.x),
        static_cast<float>(cad.aVtx[ivtx_picked].pos.y)};
      return;
    }
  }
  for (unsigned int iedge = 0; iedge < cad.aEdge.size(); ++iedge) {
    double dist = cad.aEdge[iedge].Distance(x0, y0);
    if (dist < view_height * 0.05) {
      this->iedge_picked = iedge;
      picked_pos = {
        static_cast<float>(x0),
        static_cast<float>(y0)};
      return;
    }
  }
  //
  /*
  for(unsigned int iface=0;iface<aFace.size();++iface){
    bool is_inside = aFace[iface].IsInside(x0, y0, aVec2_Tessellation);
    if( is_inside ){
      this->iface_picked = iface;
      return;
    }
  }
   */
}

DFM2_INLINE void delfem2::Cad2_Ui::DragPicked(
  CCad2D &cad,
  double p1x,
  double p1y,
  double p0x,
  double p0y) {
  if (ivtx_picked < cad.aVtx.size()) {
    cad.aVtx[ivtx_picked].pos.p[0] = p1x;
    cad.aVtx[ivtx_picked].pos.p[1] = p1y;
    cad.CopyVertexPositionsToEdges();
    return;
  }
  if (iedge_picked < cad.aEdge.size()) {
    if (ipicked_elem == 0) {
      const unsigned int iv0 = cad.topo.edges[iedge_picked].iv0;
      const unsigned int iv1 = cad.topo.edges[iedge_picked].iv1;
      cad.aVtx[iv0].pos.p[0] += p1x - p0x;
      cad.aVtx[iv0].pos.p[1] += p1y - p0y;
      cad.aVtx[iv1].pos.p[0] += p1x - p0x;
      cad.aVtx[iv1].pos.p[1] += p1y - p0y;
    } else {
      CCad2D_EdgeGeo &edge = cad.aEdge[iedge_picked];
      if (edge.type_edge == 1) {
        assert(edge.param.size() == 4);
        const CVec2d lx = (edge.p1 - edge.p0).normalized();
        const CVec2d ly = CVec2d(lx.y, -lx.x);
        if (ipicked_elem == 1) {
          edge.param[0] = (CVec2d(p1x, p1y) - edge.p0).dot(lx);
          edge.param[1] = (CVec2d(p1x, p1y) - edge.p0).dot(ly);
        } else if (ipicked_elem == 2) {
          edge.param[2] = (CVec2d(p1x, p1y) - edge.p1).dot(lx);
          edge.param[3] = (CVec2d(p1x, p1y) - edge.p1).dot(ly);
        }
      }
    }
    cad.CopyVertexPositionsToEdges();
    return;
  }
  if (iface_picked < cad.topo.faces.size()) {
    std::vector<unsigned int> aIdV = cad.topo.loops[iface_picked].GetArray_IdVertex(cad.topo.edges);
    for (unsigned int iv1: aIdV) {
      cad.aVtx[iv1].pos.p[0] += p1x - p0x;
      cad.aVtx[iv1].pos.p[1] += p1y - p0y;
    }
    cad.CopyVertexPositionsToEdges();
    return;
  }
}