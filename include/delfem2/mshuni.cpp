/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mshuni.h"

#include <vector>
#include <cassert>
#include <stack>
#include <set>
#include <climits>

// ---------------------------------------------

DFM2_INLINE void delfem2::JArray_ElSuP_MeshElem(
    std::vector<unsigned int> &elsup_ind,
    std::vector<unsigned int> &elsup,
    // ----------
    const unsigned int *elem_vtx_idx,
    size_t num_elem,
    unsigned int num_vtx_par_elem,
    size_t num_vtx) {
  elsup_ind.assign(num_vtx + 1, 0);
  for (unsigned int ielem = 0; ielem < num_elem; ielem++) {
    for (unsigned int inoel = 0; inoel < num_vtx_par_elem; inoel++) {
      const unsigned int ino1 = elem_vtx_idx[ielem * num_vtx_par_elem + inoel];
      elsup_ind[ino1 + 1] += 1;
    }
  }
  for (unsigned int ino = 0; ino < num_vtx; ++ino) {
    elsup_ind[ino + 1] += elsup_ind[ino];
  }
  unsigned int nelsup = elsup_ind[num_vtx];
  elsup.resize(nelsup);
  for (unsigned int ielem = 0; ielem < num_elem; ielem++) {
    for (unsigned int inoel = 0; inoel < num_vtx_par_elem; inoel++) {
      int unsigned ino1 = elem_vtx_idx[ielem * num_vtx_par_elem + inoel];
      unsigned int ind1 = elsup_ind[ino1];
      elsup[ind1] = ielem;
      elsup_ind[ino1] += 1;
    }
  }
  for (int ino = static_cast<int>(num_vtx); ino >= 1; --ino) {
    elsup_ind[ino] = elsup_ind[ino - 1];
  }
  elsup_ind[0] = 0;
}

DFM2_INLINE unsigned int delfem2::FindAdjEdgeIndex(
    unsigned int itri0,
    unsigned int ied0,
    unsigned int jtri0,
    const unsigned int *aTri) {
  const unsigned int iv0 = aTri[itri0 * 3 + (ied0 + 1) % 3];
  const unsigned int iv1 = aTri[itri0 * 3 + (ied0 + 2) % 3];
  assert(iv0 != iv1);
  assert(jtri0 != UINT_MAX);
  if (aTri[jtri0 * 3 + 1] == iv1 && aTri[jtri0 * 3 + 2] == iv0) { return 0; }
  if (aTri[jtri0 * 3 + 2] == iv1 && aTri[jtri0 * 3 + 0] == iv0) { return 1; }
  if (aTri[jtri0 * 3 + 0] == iv1 && aTri[jtri0 * 3 + 1] == iv0) { return 2; }
  return UINT_MAX;
}

DFM2_INLINE void delfem2::ElemQuad_DihedralTri(
    std::vector<unsigned int> &quad_vtx_idx,
    const unsigned int *tri_vtx_idx,
    const size_t num_tri,
    const size_t num_vtx) {
  std::vector<unsigned int> aElSuEl;
  ElSuEl_MeshElem(
      aElSuEl,
      tri_vtx_idx, num_tri, MESHELEM_TRI,
      num_vtx);
  assert(aElSuEl.size() == num_tri * 3);
  for (unsigned int itri = 0; itri < num_tri; ++itri) {
    for (int iedtri = 0; iedtri < 3; ++iedtri) {
      const unsigned int jtri = aElSuEl[itri * 3 + iedtri];
      if (jtri == UINT_MAX) continue; // on the boundary
      if (jtri < itri) continue;
      const unsigned int jedtri = FindAdjEdgeIndex(itri, iedtri, jtri, tri_vtx_idx);
      assert(jedtri != UINT_MAX);
      const unsigned int ipo0 = tri_vtx_idx[itri * 3 + iedtri];
      const unsigned int ipo1 = tri_vtx_idx[jtri * 3 + jedtri];
      const unsigned int ipo2 = tri_vtx_idx[itri * 3 + (iedtri + 1) % 3];
      const unsigned int ipo3 = tri_vtx_idx[itri * 3 + (iedtri + 2) % 3];
      assert(tri_vtx_idx[jtri * 3 + (jedtri + 2) % 3] == ipo2);
      assert(tri_vtx_idx[jtri * 3 + (jedtri + 1) % 3] == ipo3);
      quad_vtx_idx.push_back(ipo0);
      quad_vtx_idx.push_back(ipo1);
      quad_vtx_idx.push_back(ipo2);
      quad_vtx_idx.push_back(ipo3);
    }
  }
}

// ---------------------------------

DFM2_INLINE void delfem2::convert2Tri_Quad(
    std::vector<unsigned int> &aTri,
    const std::vector<unsigned int> &aQuad) {
  const size_t nq = aQuad.size() / 4;
  aTri.resize(nq * 6);
  for (unsigned int iq = 0; iq < nq; ++iq) {
    const unsigned int i0 = aQuad[iq * 4 + 0];
    const unsigned int i1 = aQuad[iq * 4 + 1];
    const unsigned int i2 = aQuad[iq * 4 + 2];
    const unsigned int i3 = aQuad[iq * 4 + 3];
    aTri[iq * 6 + 0] = i0;
    aTri[iq * 6 + 1] = i1;
    aTri[iq * 6 + 2] = i2;
    aTri[iq * 6 + 3] = i2;
    aTri[iq * 6 + 4] = i3;
    aTri[iq * 6 + 5] = i0;
  }
}

DFM2_INLINE void delfem2::FlipElement_Tri(
    std::vector<unsigned int> &tri_vtx_ind) {
  for (unsigned int itri = 0; itri < tri_vtx_ind.size() / 3; itri++) {
    //    int i0 = aTri[itri*3+0];
    unsigned int i1 = tri_vtx_ind[itri * 3 + 1];
    unsigned int i2 = tri_vtx_ind[itri * 3 + 2];
    tri_vtx_ind[itri * 3 + 1] = i2;
    tri_vtx_ind[itri * 3 + 2] = i1;
  }
}


// -------------------------------------

DFM2_INLINE void delfem2::JArray_ElSuP_MeshTri(
    std::vector<unsigned int> &elsup_ind,
    std::vector<unsigned int> &elsup,
    // --
    const std::vector<unsigned int> &aTri,
    int nXYZ) {
  JArray_ElSuP_MeshElem(
      elsup_ind, elsup,
      aTri.data(), aTri.size() / 3, 3,
      nXYZ);
}

// ----------------------------------------------------------------------------------------------------------

DFM2_INLINE void delfem2::ElSuEl_MeshElem(
    std::vector<unsigned int> &elsuel,
    const unsigned int *elem_vtx_idx,
    size_t num_elem,
    int nNoEl,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    const int num_face_par_elem,
    const int num_vtx_on_face,
    const int (*vtx_on_elem_face)[4]) {
  assert(!elsup_ind.empty());
  const std::size_t np = elsup_ind.size() - 1;

  elsuel.assign(
      num_elem * num_face_par_elem,
      UINT_MAX);

  std::vector<int> flg_point(np, 0);
  std::vector<unsigned int> inpofa(num_vtx_on_face);
  for (unsigned int iel = 0; iel < num_elem; iel++) {
    for (int ifael = 0; ifael < num_face_par_elem; ifael++) {
      for (int ipofa = 0; ipofa < num_vtx_on_face; ipofa++) {
        int int0 = vtx_on_elem_face[ifael][ipofa];
        const unsigned int ip = elem_vtx_idx[iel * nNoEl + int0];
        assert(ip < np);
        inpofa[ipofa] = ip;
        flg_point[ip] = 1;
      }
      const unsigned int ipoin0 = inpofa[0];
      bool iflg = false;
      for (unsigned int ielsup = elsup_ind[ipoin0]; ielsup < elsup_ind[ipoin0 + 1]; ielsup++) {
        const unsigned int jelem0 = elsup[ielsup];
        if (jelem0 == iel) continue;
        for (int jfael = 0; jfael < num_face_par_elem; jfael++) {
          iflg = true;
          for (int jpofa = 0; jpofa < num_vtx_on_face; jpofa++) {
            int jnt0 = vtx_on_elem_face[jfael][jpofa];
            const unsigned int jpoin0 = elem_vtx_idx[jelem0 * nNoEl + jnt0];
            if (flg_point[jpoin0] == 0) {
              iflg = false;
              break;
            }
          }
          if (iflg) {
            elsuel[iel * num_face_par_elem + ifael] = jelem0;
            break;
          }
        }
        if (iflg) break;
      }
      if (!iflg) {
        elsuel[iel * num_face_par_elem + ifael] = UINT_MAX;
      }
      for (int ipofa = 0; ipofa < num_vtx_on_face; ipofa++) {
        flg_point[inpofa[ipofa]] = 0;
      }
    }
  }
}

/*
void makeSurroundingRelationship
(std::vector<int>& aElSurRel,
 const int* aEl, int nEl, int nNoEl,
 MESHELEM_TYPE type,
 const std::vector<int>& elsup_ind,
 const std::vector<int>& elsup)
{
  const int nfael = nFaceElem(type);
  const int nnofa = nNodeElemFace(type, 0);
  assert( nNoEl == nNodeElem(type) );
  makeSurroundingRelationship(aElSurRel,
                              aEl, nEl, nNoEl,
                              elsup_ind, elsup,
                              nfael, nnofa, noelElemFace(type));
}
 */

DFM2_INLINE void delfem2::ElSuEl_MeshElem(
    std::vector<unsigned int> &aElSuEl,
    const unsigned int *elem_vtx_idx,
    size_t num_elem,
    MESHELEM_TYPE type,
    const size_t num_vtx) {
  const int nNoEl = nNodeElem(type);
  std::vector<unsigned int> elsup_ind, elsup;
  JArray_ElSuP_MeshElem(
      elsup_ind, elsup,
      elem_vtx_idx, num_elem, nNoEl, num_vtx);
  const int nfael = nFaceElem(type);
  const int nnofa = nNodeElemFace(type, 0);
  ElSuEl_MeshElem(
      aElSuEl,
      elem_vtx_idx, num_elem, nNoEl,
      elsup_ind, elsup,
      nfael, nnofa, noelElemFace(type));
  assert(aElSuEl.size() == num_elem * nfael);
}


// -------------------------------------------------------------------------

DFM2_INLINE void delfem2::JArrayPointSurPoint_MeshOneRingNeighborhood(
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    //
    const unsigned int *elem_vtx_idx,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    unsigned int num_vtx_par_elem,
    size_t num_vtx) {
  std::vector<unsigned int> aflg(num_vtx, UINT_MAX);
  psup_ind.assign(num_vtx + 1, 0);
  for (unsigned int ipoint = 0; ipoint < num_vtx; ipoint++) {
    aflg[ipoint] = ipoint;
    for (unsigned int ielsup = elsup_ind[ipoint]; ielsup < elsup_ind[ipoint + 1]; ielsup++) {
      unsigned int jelem = elsup[ielsup];
      for (unsigned int jnoel = 0; jnoel < num_vtx_par_elem; jnoel++) {
        unsigned int jnode = elem_vtx_idx[jelem * num_vtx_par_elem + jnoel];
        if (aflg[jnode] != ipoint) {
          aflg[jnode] = ipoint;
          psup_ind[ipoint + 1]++;
        }
      }
    }
  }
  for (unsigned int ipoint = 0; ipoint < num_vtx; ipoint++) {
    psup_ind[ipoint + 1] += psup_ind[ipoint];
  }
  const unsigned int npsup = psup_ind[num_vtx];
  psup.resize(npsup);
  for (unsigned int ipoint = 0; ipoint < num_vtx; ipoint++) { aflg[ipoint] = UINT_MAX; }
  for (unsigned int ipoint = 0; ipoint < num_vtx; ipoint++) {
    aflg[ipoint] = ipoint;
    for (unsigned int ielsup = elsup_ind[ipoint]; ielsup < elsup_ind[ipoint + 1]; ielsup++) {
      unsigned int jelem = elsup[ielsup];
      for (unsigned int jnoel = 0; jnoel < num_vtx_par_elem; jnoel++) {
        unsigned int jnode = elem_vtx_idx[jelem * num_vtx_par_elem + jnoel];
        if (aflg[jnode] != ipoint) {
          aflg[jnode] = ipoint;
          const unsigned int ind = psup_ind[ipoint];
          psup[ind] = jnode;
          psup_ind[ipoint]++;
        }
      }
    }
  }
  for (int ipoint = (int) num_vtx; ipoint > 0; --ipoint) {
    psup_ind[ipoint] = psup_ind[ipoint - 1];
  }
  psup_ind[0] = 0;
}

DFM2_INLINE void delfem2::JArray_PSuP_MeshElem(
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    //
    const unsigned int *elem_vtx_idx,
    size_t num_elm,
    unsigned int num_vtx_par_elem,
    size_t num_vtx) {
  std::vector<unsigned int> elsup_ind, elsup;
  JArray_ElSuP_MeshElem(
      elsup_ind, elsup,
      elem_vtx_idx, num_elm, num_vtx_par_elem, num_vtx);
  JArrayPointSurPoint_MeshOneRingNeighborhood(
      psup_ind, psup,
      elem_vtx_idx, elsup_ind, elsup, num_vtx_par_elem, num_vtx);
}

DFM2_INLINE void delfem2::makeOneRingNeighborhood_TriFan(
    std::vector<int> &psup_ind,
    std::vector<int> &psup,
    // ----------------------
    const std::vector<int> &aTri,
    const std::vector<int> &aTriSurRel,
    const std::vector<int> &elsup_ind,
    const std::vector<int> &elsup,
    int npoint) {
  psup_ind.resize(npoint + 1);
  psup_ind[0] = 0;
  psup.clear();
  for (int ipoint = 0; ipoint < npoint; ++ipoint) {
    int iel0 = -1;
    int inoel0 = -1;
    {
      int ielsup0 = elsup_ind[ipoint];
      iel0 = elsup[ielsup0];
      if (aTri[iel0 * 3 + 0] == ipoint) { inoel0 = 0; }
      if (aTri[iel0 * 3 + 1] == ipoint) { inoel0 = 1; }
      if (aTri[iel0 * 3 + 2] == ipoint) { inoel0 = 2; }
      assert(inoel0 != -1);
    }
    int iel_cur = iel0;
    int inoel_cur = inoel0;
    for (;;) {
      int jnoel_cur = (inoel_cur + 1) % 3;
      int jp0 = aTri[iel_cur * 3 + jnoel_cur];
      psup.push_back(jp0);
      int iel_next = aTriSurRel[iel_cur * 6 + 2 * jnoel_cur + 0];
      int inoel_next = -1;
      if (aTri[iel_next * 3 + 0] == ipoint) { inoel_next = 0; }
      if (aTri[iel_next * 3 + 1] == ipoint) { inoel_next = 1; }
      if (aTri[iel_next * 3 + 2] == ipoint) { inoel_next = 2; }
      assert(inoel_next != -1);
      if (iel_next == iel0) break;
      iel_cur = iel_next;
      inoel_cur = inoel_next;
    }
    psup_ind[ipoint + 1] = (int) psup.size();
  }
}

DFM2_INLINE void delfem2::JArrayEdge_MeshElem(
    std::vector<unsigned int> &edge_ind,
    std::vector<unsigned int> &edge,
    //
    const unsigned int *elem_vtx_idx,
    MESHELEM_TYPE elem_type,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    bool is_bidirectional) {
  const int neElm = mapMeshElemType2NEdgeElem[elem_type];
  const int nnoelElm = mapMeshElemType2NNodeElem[elem_type];
  const int (*aNoelEdge)[2] = noelElemEdge(elem_type);
  const std::size_t nPoint0 = elsup_ind.size() - 1;
  edge_ind.resize(nPoint0 + 1);
  edge_ind[0] = 0;
  for (unsigned int ip = 0; ip < nPoint0; ++ip) {
    std::set<unsigned int> setIP;
    for (unsigned int ielsup = elsup_ind[ip]; ielsup < elsup_ind[ip + 1]; ++ielsup) {
      const unsigned int iq0 = elsup[ielsup];
      for (int ie = 0; ie < neElm; ++ie) {
        int inoel0 = aNoelEdge[ie][0];
        int inoel1 = aNoelEdge[ie][1];
        unsigned int ip0 = elem_vtx_idx[iq0 * nnoelElm + inoel0];
        unsigned int ip1 = elem_vtx_idx[iq0 * nnoelElm + inoel1];
        if (ip0 != ip && ip1 != ip) continue;
        if (ip0 == ip) {
          if (is_bidirectional || ip1 > ip) { setIP.insert(ip1); }
        } else {
          if (is_bidirectional || ip0 > ip) { setIP.insert(ip0); }
        }
      }
    }
    for (unsigned int itr : setIP) {
      edge.push_back(itr);
    }
    edge_ind[ip + 1] = edge_ind[ip] + static_cast<unsigned int>(setIP.size());
  }
}

DFM2_INLINE void delfem2::MeshLine_JArrayEdge(
    std::vector<unsigned int> &line_vtx_idx,
    //
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup) {
  line_vtx_idx.reserve(psup.size() * 2);
  const std::size_t np = psup_ind.size() - 1;
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      unsigned int jp = psup[ipsup];
      line_vtx_idx.push_back(ip);
      line_vtx_idx.push_back(jp);
    }
  }
}

DFM2_INLINE void delfem2::MeshLine_MeshElem(
    std::vector<unsigned int> &line_vtx_idx,
    const unsigned int *elem_vtx_idx,
    size_t num_elem,
    MESHELEM_TYPE elem_type,
    size_t num_vtx) {
  std::vector<unsigned int> elsup_ind, elsup;
  const unsigned int nPoEl = mapMeshElemType2NNodeElem[elem_type];
  JArray_ElSuP_MeshElem(
      elsup_ind, elsup,
      elem_vtx_idx, num_elem, nPoEl, num_vtx);
  std::vector<unsigned int> edge_ind, edge;
  JArrayEdge_MeshElem(
      edge_ind, edge,
      elem_vtx_idx,
      elem_type,
      elsup_ind, elsup, false);
  MeshLine_JArrayEdge(
      line_vtx_idx,
      edge_ind, edge);
}

// ---------------------------------------

DFM2_INLINE void delfem2::MarkConnectedElements(
    std::vector<unsigned int> &aFlagElem,
    unsigned int itri_ker,
    unsigned int igroup,
    const std::vector<unsigned int> &aElSuEl) {
  const std::size_t nel = aFlagElem.size();
  const std::size_t nfael = aElSuEl.size() / nel;
  aFlagElem[itri_ker] = igroup;
  std::stack<unsigned int> next;
  next.push(itri_ker);
  while (!next.empty()) {
    unsigned int itri0 = next.top();
    next.pop();
    for (unsigned int ie = 0; ie < nfael; ++ie) {
      const unsigned int ita = aElSuEl[itri0 * nfael + ie];
      if (ita == UINT_MAX) continue;
      if (aFlagElem[ita] != igroup) {
        aFlagElem[ita] = igroup;
        next.push(ita);
      }
    }
  }
}

DFM2_INLINE void delfem2::MakeGroupElem(
    int &ngroup,
    std::vector<unsigned int> &aIndGroup,
    const std::vector<unsigned int> &aTri,
    const std::vector<unsigned int> &aTriSurRel,
    [[maybe_unused]] const int nfael,
    const int nnoel) {
  const std::size_t nelem = aTri.size() / nnoel;
  aIndGroup.assign(nelem, UINT_MAX);
  int igroup = -1;
  for (;;) {
    unsigned int itri_ker = 0;
    for (; itri_ker < nelem; ++itri_ker) {
      if (aIndGroup[itri_ker] == UINT_MAX) break;
    }
    if (itri_ker == nelem) break;
    igroup++;
    MarkConnectedElements(aIndGroup, itri_ker, igroup, aTriSurRel);
  }
  ngroup = igroup + 1;
}
