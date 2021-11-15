//
// Created by Nobuyuki Umetani on 2021-11-07.
//

#ifndef CAD2_MESH_DEFORMATION_H_
#define CAD2_MESH_DEFORMATION_H_

#include "delfem2/cagedef.h"
#include "delfem2/cad2.h"
#include "delfem2/cad2_mesher.h"

void SetCadMeshDeformationWeightVtx(
    std::vector<double> &vtx_w,
    unsigned int vertex_idx,
    const delfem2::CCad2D &cad,
    const delfem2::CMesher_Cad2D &mesher,
    const std::vector<double> &vtx_xy) {
  namespace dfm2 = delfem2;
  vtx_w.assign(vtx_xy.size() / 2, 0.);
  if (vertex_idx >= cad.nVtx()) { return; }
  std::vector<unsigned int> face_ids = cad.topo.FindFaceIndexes_IncludeVeretx(vertex_idx);
  if (face_ids.size() != 1) { return; }
  unsigned int face_idx = face_ids[0];
  unsigned int iloop = cad.topo.faces[face_idx].aIL[0];
  std::vector<double> xyws;
  for (auto[ie, dir]: cad.topo.loops[iloop].aIE) {
    unsigned int iv0 = cad.topo.edges[ie].IndexVertex(dir);
    unsigned int iv1 = cad.topo.edges[ie].IndexVertex(!dir);
    {
      vtx_w[iv0] = (iv0 == vertex_idx) ? 1. : 0.;
      xyws.push_back(cad.aVtx[iv0].pos.x);
      xyws.push_back(cad.aVtx[iv0].pos.y);
      xyws.push_back(vtx_w[iv0]);
    }
    if (iv0 != vertex_idx && iv1 != vertex_idx) { continue; }
    for (auto &point: mesher.edge_point[ie]) {
      unsigned int ip = point.first;
      double w0 = (iv1 == vertex_idx) ? point.second : 1. - point.second;
      if (cad.aEdge[ie].type_edge == dfm2::CCad2D_EdgeGeo::LINE) {
        vtx_w[ip] = w0;
      } else if (cad.aEdge[ie].type_edge == dfm2::CCad2D_EdgeGeo::BEZIER_CUBIC) {
        vtx_w[ip] = 3 * w0 * w0 - 2 * w0 * w0 * w0;
      }
      xyws.push_back(vtx_xy[ip * 2 + 0]);
      xyws.push_back(vtx_xy[ip * 2 + 1]);
      xyws.push_back(vtx_w[ip]);
    }
  }
  for (unsigned int ip = 0; ip < vtx_xy.size() / 2; ++ip) {
    if (mesher.aFlgPnt[ip] != face_idx + cad.nVtx() + cad.nEdge()) { continue; }
    vtx_w[ip] = dfm2::MeanValueCoordinate_Polygon2<dfm2::CVec2d>(
        xyws, {vtx_xy[ip * 2 + 0], vtx_xy[ip * 2 + 1]});
  }
}

void SetCadMeshDeformationWeightEdge(
    std::vector<double> &vtx_w,
    unsigned int edge_idx,
    const delfem2::CCad2D &cad,
    const delfem2::CMesher_Cad2D &mesher,
    const std::vector<double> &vtx_xy) {
  namespace dfm2 = delfem2;
  vtx_w.assign(vtx_xy.size() / 2, 0.);
  if (edge_idx >= cad.nVtx()) { return; }
  std::vector<unsigned int> face_ids = cad.topo.FindFaceIndexes_IncludeEdge(edge_idx);
  std::cout << face_ids.size() << std::endl;
  assert(face_ids.size() == 1);
  unsigned int face_idx = face_ids[0];
  unsigned int iloop = cad.topo.faces[face_idx].aIL[0];
  std::vector<double> xyws;
  for (auto[ie, dir]: cad.topo.loops[iloop].aIE) {
    unsigned int iv0 = cad.topo.edges[ie].IndexVertex(dir);
    unsigned int iv1 = cad.topo.edges[ie].IndexVertex(!dir);
    if (iv0 == cad.topo.edges[edge_idx].IndexVertex(true) ||
        iv0 == cad.topo.edges[edge_idx].IndexVertex(false)) {
      vtx_w[iv0] = 1.;
      xyws.push_back(cad.aVtx[iv0].pos.x);
      xyws.push_back(cad.aVtx[iv0].pos.y);
      xyws.push_back(vtx_w[iv0]);
    }
    for (auto &point: mesher.edge_point[ie]) {
      unsigned int ip = point.first;
      double w0 = 0.0;
      if (ie == edge_idx) {
        w0 = 1.0;
      } else if (cad.topo.edges[edge_idx].IndexVertex(true) == iv1) {
        w0 = point.second;
      } else if (cad.topo.edges[edge_idx].IndexVertex(false) == iv0) {
        w0 = 1.0 - point.second;
      }
      if (cad.aEdge[ie].type_edge == dfm2::CCad2D_EdgeGeo::LINE) {
        vtx_w[ip] = w0;
      } else if (cad.aEdge[ie].type_edge == dfm2::CCad2D_EdgeGeo::BEZIER_CUBIC) {
        vtx_w[ip] = 3 * w0 * w0 - 2 * w0 * w0 * w0;
      }
      xyws.push_back(vtx_xy[ip * 2 + 0]);
      xyws.push_back(vtx_xy[ip * 2 + 1]);
      xyws.push_back(vtx_w[ip]);
    }
  }
  for (unsigned int ip = 0; ip < vtx_xy.size() / 2; ++ip) {
    if (mesher.aFlgPnt[ip] != face_idx + cad.nVtx() + cad.nEdge()) { continue; }
    vtx_w[ip] = dfm2::MeanValueCoordinate_Polygon2<dfm2::CVec2d>(
        xyws, {vtx_xy[ip * 2 + 0], vtx_xy[ip * 2 + 1]});
  }
}

#endif //CAD2_MESH_DEFORMATION_H_
