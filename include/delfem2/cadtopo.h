/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_CADTOPO_H
#define DFM2_CADTOPO_H

#include <vector>
#include <climits>

namespace delfem2 {

/**
 * edge is directed
 */
class CadTopo_Edge {
 public:
  unsigned int IndexVertex(bool is_root) const { return is_root ? iv0 : iv1; }
 public:
  unsigned int iv0, iv1;
};

// --------------------

class CadTopo_Loop {
 public:
  CadTopo_Loop() = default;

  std::vector<unsigned int> GetArray_IdVertex(
      const std::vector<CadTopo_Edge> &edges) const {
    std::vector<unsigned int> res;
    for ( auto [ie0,dir] : aIE ) {
      assert( ie0 < edges.size() );
      res.push_back( edges[ie0].IndexVertex(dir) );
    }
    return res;
  }

  void Assert(
      const std::vector<CadTopo_Edge> &edges) const {
    assert( iv == UINT_MAX || aIE.empty() );
    const size_t ne = aIE.size();
    for (unsigned int iie = 0; iie < ne; ++iie) {
      const unsigned int ie0 = aIE[(iie + 0) % ne].first;
      const unsigned int ie1 = aIE[(iie + 1) % ne].first;
      assert( ie0 < edges.size() && ie1 < edges.size());
      const bool flg0 = aIE[(iie + 0) % ne].second;
      const bool flg1 = aIE[(iie + 1) % ne].second;
      // (root of ie1) == (end of ie0)
      assert( edges[ie0].IndexVertex(!flg0) == edges[ie1].IndexVertex(flg1) );
    }
  }

  bool IsIncludeVertex(
      unsigned int ivtx,
      const std::vector<CadTopo_Edge> &edges) const {
    if( iv != UINT_MAX ){
      assert( aIE.empty() );
      return true;
    }
    for (auto [ie0,dir] : aIE) {
      assert( ie0 < edges.size() );
      if( ivtx == edges[ie0].IndexVertex(dir) ){ return true; }
    }
    return false;
  }

  bool IsIncludeEdge(
      unsigned int iedge) const {
    for (auto [ie0,dir] : aIE) {
      if( ie0 == iedge ){ return true; }
    }
    return false;
  }

 public:
  unsigned int iv = UINT_MAX;
  std::vector<std::pair<unsigned int, bool> > aIE; // index of edge, is this edge ccw or not
};

// ---------------------------

class CadTopo_Face {
 public:
  bool IsIncludeVertex(
      unsigned int ivtx,
      const std::vector<CadTopo_Edge> &edges,
      const std::vector<CadTopo_Loop> &loops) const {
    for(unsigned int il : aIL ){
      assert( il < loops.size() );
      if( loops[il].IsIncludeVertex(ivtx, edges) ){ return true; }
    }
    return false;
  }

  bool IsIncludeEdge(
    unsigned int iedge,
    const std::vector<CadTopo_Loop> &loops) const {
      for(unsigned int il : aIL ){
        assert( il < loops.size() );
        if( loops[il].IsIncludeEdge(iedge) ){ return true; }
      }
      return false;
  }
 public:
  std::vector<unsigned int> aIL;  // island loop indeces
};

// ---------------------------

class CadTopo {
public:
  CadTopo() {
    num_vertex = 0;
  }

  void Clear() {
    num_vertex = 0;
    edges.clear();
    loops.clear();
    faces.clear();
  }

  void AddPolygon(unsigned int np) {
    const unsigned int iv0 = num_vertex;
    num_vertex += np;
    const unsigned int ie0 = static_cast<unsigned int>(edges.size());
    for (unsigned int iie = 0; iie < np; ++iie) {
      edges.emplace_back( CadTopo_Edge{
          iv0 + (iie + 0) % np,
          iv0 + (iie + 1) % np } );
    }
    { // loop
      CadTopo_Loop loop0;
      for (unsigned int iie = 0; iie < np; ++iie) {
        loop0.aIE.push_back(std::make_pair(ie0 + iie, true));
      }
      loops.push_back(loop0);
    }
    faces.emplace_back( CadTopo_Face{
      {static_cast<unsigned int>(loops.size() - 1)} } );
  }

  bool AddVtx_Face(unsigned int ifc) {
    if (ifc >= faces.size()) { return false; }
    const unsigned int ivn = num_vertex;
    num_vertex += 1;
    CadTopo_Loop loop0;
    loop0.iv = ivn;
    loops.push_back(loop0);
    const unsigned int il0 = static_cast<unsigned int>(loops.size() - 1);
    faces[ifc].aIL.push_back(il0);
    return true;
  }

  bool AddVtx_Edge(unsigned int ieo) {
    if (ieo >= edges.size()) { return false; }
    const unsigned int ivn = num_vertex;
    num_vertex += 1;
    const unsigned int iv0 = edges[ieo].iv0;
    const unsigned int iv1 = edges[ieo].iv1;
    const unsigned int ien = static_cast<unsigned int>(edges.size());
    edges.resize(edges.size() + 1);
    edges[ieo].iv0 = iv0;
    edges[ieo].iv1 = ivn;
    edges[ien].iv0 = ivn;
    edges[ien].iv1 = iv1;
    for (unsigned int il = 0; il < loops.size(); ++il) {
      const size_t ne = loops[il].aIE.size();
      unsigned int iie = 0;
      for (; iie < ne; ++iie) {
        if (loops[il].aIE[iie].first == ieo) { break; }
      }
      if (iie == ne) { continue; }
      if (loops[il].aIE[iie].second) {
        loops[il].aIE.insert(
            loops[il].aIE.begin() + iie + 1,
            std::make_pair(ien, true));
      } else {
        std::cerr << "TODO: implement this" << std::endl;
      }
    }
    return true;
  }

  void Assert() const {
    for (const auto& l : loops) {
      l.Assert(edges);
    }
  }

  std::pair<unsigned int, unsigned int> VertexIndexs_EdgeEndPoints(unsigned int edge_idx) const {
    assert(edge_idx < edges.size());
    return { edges[edge_idx].iv0, edges[edge_idx].iv1 };
  }

  std::vector<unsigned int> VertexIndexs_Face(
      unsigned int iface) const
  {
    assert(iface < faces.size());
    return loops[iface].GetArray_IdVertex(this->edges);
  }

  std::vector<std::pair<unsigned int,bool> > EdgeIndexes_Face(unsigned int iface) const
  {
    std::vector<std::pair<unsigned int,bool> > aIdE;
    for(unsigned int il : faces[iface].aIL) {
      aIdE.insert( aIdE.end(), loops[il].aIE.begin(), loops[il].aIE.end() );
    }
    return aIdE;
  }

  std::vector<unsigned int> FindFaceIndexes_IncludeVeretx(unsigned int ivtx) const {
    assert(ivtx < num_vertex);
    std::vector<unsigned int> res;
    for (unsigned int ifc = 0; ifc < faces.size(); ++ifc) {
      if( !faces[ifc].IsIncludeVertex(ivtx,edges,loops) ){ continue; }
      res.push_back(ifc);
    }
    return res;
  }

  std::vector<unsigned int> FindFaceIndexes_IncludeEdge(unsigned int iedge) const {
    assert(iedge < edges.size());
    std::vector<unsigned int> res;
    for (unsigned int ifc = 0; ifc < faces.size(); ++ifc) {
      if( !faces[ifc].IsIncludeEdge(iedge,loops) ){ continue; }
      res.push_back(ifc);
    }
    return res;
  }

public:
  unsigned int num_vertex;
  std::vector<CadTopo_Edge> edges;
  std::vector<CadTopo_Loop> loops;
  std::vector<CadTopo_Face> faces;
};

} // namespace delfem2


#endif /* DFM2_CADTOPO_H */
