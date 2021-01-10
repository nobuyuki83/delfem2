/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_CADTOPO_H
#define DFM2_CADTOPO_H

namespace delfem2 {

class CCadTopo {
public:
  CCadTopo() {
    nVertex = 0;
  }

  void Clear() {
    nVertex = 0;
    aEdge.clear();
    aLoop.clear();
    aFace.clear();
  }

  void AddPolygon(unsigned int np) {
    const int iv0 = nVertex;
    nVertex += np;
    const unsigned int ie0 = static_cast<unsigned int>(aEdge.size());
    for (unsigned int iie = 0; iie < np; ++iie) {
      CEdge edge0;
      edge0.iv0 = iv0 + (iie + 0) % np;
      edge0.iv1 = iv0 + (iie + 1) % np;
      aEdge.push_back(edge0);
    }
    { // loop
      CLoop loop0;
      for (unsigned int iie = 0; iie < np; ++iie) {
        loop0.aIE.push_back(std::make_pair(ie0 + iie, true));
      }
      aLoop.push_back(loop0);
    }
    { // face
      const unsigned int il0 = static_cast<unsigned int>(aLoop.size() - 1);
      CFace face0;
      face0.aIL.push_back(il0);
      aFace.push_back(face0);
    }
  }

  bool AddVtx_Face(unsigned int ifc) {
    if (ifc >= aFace.size()) { return false; }
    const int ivn = nVertex;
    nVertex += 1;
    CLoop loop0;
    loop0.iv = ivn;
    aLoop.push_back(loop0);
    const unsigned int il0 = static_cast<unsigned int>(aLoop.size() - 1);
    aFace[ifc].aIL.push_back(il0);
    return true;
  }

  bool AddVtx_Edge(unsigned int ieo) {
    if (ieo >= aEdge.size()) { return false; }
    const int ivn = nVertex;
    nVertex += 1;
    const int iv0 = aEdge[ieo].iv0;
    const int iv1 = aEdge[ieo].iv1;
    const unsigned int ien = static_cast<unsigned int>(aEdge.size());
    aEdge.resize(aEdge.size() + 1);
    aEdge[ieo].iv0 = iv0;
    aEdge[ieo].iv1 = ivn;
    aEdge[ien].iv0 = ivn;
    aEdge[ien].iv1 = iv1;
    for (unsigned int il = 0; il < aLoop.size(); ++il) {
      const size_t ne = aLoop[il].aIE.size();
      unsigned int iie = 0;
      for (; iie < ne; ++iie) {
        if (aLoop[il].aIE[iie].first == (int) ieo) { break; }
      }
      if (iie == ne) { continue; }
      if (aLoop[il].aIE[iie].second) {
        aLoop[il].aIE.insert(aLoop[il].aIE.begin() + iie + 1, std::make_pair(ien, true));
      } else {
        std::cout << "TODO: implement this" << std::endl;
      }
    }
    return true;
  }

  bool Check() const {
    for (unsigned int il = 0; il < aLoop.size(); ++il) {
      if (!aLoop[il].Check(aEdge)) {
        assert(0);
        return false;
      }
    }
    return true;
  }

public:
  class CEdge {
  public:
    int iv0, iv1;
  };

  class CLoop {
  public:
    CLoop() {
      iv = -1;
    }

    std::vector<int> GetArray_IdVertex(const std::vector<CEdge> &aEdge) const {
      std::vector<int> res;
      for (unsigned int ie = 0; ie < aIE.size(); ++ie) {
        const int ie0 = aIE[ie].first;
        const bool dir = aIE[ie].second;
        const int iv0 = dir ? aEdge[ie0].iv0 : aEdge[ie0].iv1;
        res.push_back(iv0);
      }
      return res;
    }

    bool Check(const std::vector<CEdge> &aEdge) const {
      if (iv == -1 && aIE.empty()) {
        assert(0);
        return false;
      }
      if (iv != -1 && !aIE.empty()) {
        assert(0);
        return false;
      }
      if (iv == -1 && aIE.empty()) {
        assert(0);
        return false;
      }
      const size_t ne = aIE.size();
      for (unsigned int iie = 0; iie < ne; ++iie) {
        const int ie0 = aIE[(iie + 0) % ne].first;
        const int ie1 = aIE[(iie + 1) % ne].first;
        const bool flg0 = aIE[(iie + 0) % ne].second;
        const bool flg1 = aIE[(iie + 1) % ne].second;
        if (ie0 < 0 || ie0 >= (int) aEdge.size()) { return false; }
        if (ie1 < 0 || ie1 >= (int) aEdge.size()) { return false; }
        int iv0a = (flg0) ? aEdge[ie0].iv1 : aEdge[ie0].iv0;
        int iv0b = (flg1) ? aEdge[ie1].iv0 : aEdge[ie1].iv1;
        assert(iv0a == iv0b);
        if (iv0a != iv0b) return false;
      }
      return true;
    }

  public:
    int iv;
    std::vector<std::pair<int, bool> > aIE; // index of edge, is this edge ccw?
  };

  class CFace {
  public:
    std::vector<int> aIL;
  };

public:
  int nVertex;
  std::vector<CEdge> aEdge;
  std::vector<CLoop> aLoop;
  std::vector<CFace> aFace;
};

} // namespace delfem2


#endif /* DFM2_CADTOPO_H */
