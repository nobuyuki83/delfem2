/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file functions to analyze mesh topology for static meshes
 * @details the functions only care about the topology. Geometry (coordinate) information is not handled in this file
 */

// DONE(2020/12/23): change name mshuni.h
// DODO(2020/12/23): separaete jarray.h
// DONE(2020/12/12): separated mshsubdiv
// DONE(2020/12/09): separate mixed elem

#ifndef DFM2_MSH_TOPOLOGY_UNIFORM_H
#define DFM2_MSH_TOPOLOGY_UNIFORM_H

#include <cstdio>
#include <vector>

#include "delfem2/mshelm.h"
#include "delfem2/dfm2_inline.h"

namespace delfem2 {

// ---------------------------------------------------

DFM2_INLINE unsigned FindAdjEdgeIndex(
    unsigned int itri,
    unsigned int ied,
    unsigned int jtri,
    const unsigned int *aTri);

DFM2_INLINE void convert2Tri_Quad(
    std::vector<unsigned int> &tri_vtx,
    const std::vector<unsigned int> &quad_vtx);


/**
 * @brief Make quad mesh from tri mesh by merging adjacent triangle elements
 * @param[out] quad_vtx_idx element index of quad mesh
 * @param[in] tri_vtx_idx element index of tri mesh
 * @param[in] num_tri number of triangle mesh
 * @param[in] num_vtx number of points
 */
DFM2_INLINE void ElemQuad_DihedralTri(
    std::vector<unsigned int> &quad_vtx_idx,
    const unsigned int *tri_vtx_idx,
    size_t num_tri,
    size_t num_vtx);

DFM2_INLINE void FlipElement_Tri(
    std::vector<unsigned int> &tri_vtx_ind);


/**
 * make elem surrounding point
 */
DFM2_INLINE void JArray_ElSuP_MeshElem(
    std::vector<unsigned int> &elsup_ind,
    std::vector<unsigned int> &elsup,
    //
    const unsigned int *elem_vtx_idx,
    size_t num_elem,
    unsigned int num_vtx_par_elem,
    size_t num_vtx);

/**
 * @brief make elem surrounding point for triangle mesh
 */
DFM2_INLINE void JArray_ElSuP_MeshTri(
    std::vector<unsigned int> &elsup_ind,
    std::vector<unsigned int> &elsup,
    const std::vector<unsigned int> &tri_vtx,
    size_t num_vtx);



// -----------------
// elem sur elem

/**
 * @brief compute adjacent element index for mesh element
 * @param[out] elsuel neighbouring element index (UINT_MAX for boundary)
 * @param[in] elem_vtx_idx array of connectivity
 * @param[in] num_elem number of elements
 * @param[in] nNoEl number of nodes in a element
 * @param[in] elsup_ind jagged array index of "elem surrounding point"
 * @param[in] elsup jagged array value of "elem surrounding point"
 * @param[in] num_face_par_elem number of neibouring elements
 * @param[in] num_vtx_on_face how many nodes are shared with a nighbouring element
 * @param vtx_on_elem_face
 */
DFM2_INLINE void ElSuEl_MeshElem(
    std::vector<unsigned int> &elsuel,
    const unsigned int *elem_vtx_idx,
    size_t num_elem,
    int nNoEl,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    const int num_face_par_elem,
    const int num_vtx_on_face,
    const int (*vtx_on_elem_face)[4]);

/**
 * @brief compute adjacent element index for mesh element
 * @param[out] aElSuEl adjacent element index for element edge/face (UINT_MAX if face/edge is on the boundary)
 * @param[in] elem_vtx_idx elemnet index
 * @param[in] num_elem number of elements
 * @param[in] type type of element
 * @param[in] num_vtx number of points
 */
DFM2_INLINE void ElSuEl_MeshElem(
    std::vector<unsigned int> &aElSuEl,
    const unsigned int *elem_vtx_idx,
    size_t num_elem,
    delfem2::MESHELEM_TYPE type,
    const size_t num_vtx);

/**
 * @brief make point surrounding point
 * @details psup -> edge bidirectional
 * edge unidir (ip0<ip1)
 * line (array of 2)
 */
DFM2_INLINE void JArrayPointSurPoint_MeshOneRingNeighborhood(
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    //
    const unsigned int *elem_vtx,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    unsigned int num_vtx_par_elem,
    size_t num_vtx);

/**
 * @brief compute indexes of points surrounding a point as a jagged array
 * @param num_vtx_par_elem number of nodes in an element
 */
DFM2_INLINE void JArray_PSuP_MeshElem(
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    //
    const unsigned int *elem_vtx,
    size_t num_elm,
    unsigned int num_vtx_par_elem,
    size_t num_vtx);

DFM2_INLINE void makeOneRingNeighborhood_TriFan(
    std::vector<int> &psup_ind,
    std::vector<int> &psup,
    //
    const std::vector<int> &aTri,
    const std::vector<int> &aTriSurRel,
    const std::vector<int> &elsup_ind,
    const std::vector<int> &elsup,
    int np);

DFM2_INLINE void JArrayEdge_MeshElem(
    std::vector<unsigned int> &edge_ind,
    std::vector<unsigned int> &edge,
    //
    const unsigned int *elem_vtx,
    delfem2::MESHELEM_TYPE elem_type,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup,
    bool is_bidirectional);

DFM2_INLINE void MeshLine_JArrayEdge(
    std::vector<unsigned int> &line_vtx,
    //
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup);

/**
 * Extracting line element from mesh (e.g., triangle mesh).
 * The edges of the element becomes the line
 * @param[out] line_vtx
 * @param[in] elem_vtx
 * @param[in] num_elem
 * @param[in] elem_type
 * @param[in] num_vtx the number of the vertices
 */
DFM2_INLINE void MeshLine_MeshElem(
    std::vector<unsigned int> &line_vtx,
    //
    const unsigned int *elem_vtx,
    size_t num_elem,
    delfem2::MESHELEM_TYPE elem_type,
    size_t num_vtx);

// ------------------------------------

DFM2_INLINE unsigned int FindIndexTri(
    const unsigned int *tri_vtx,
    unsigned int ixyz);

// ------------------------------------

/**
 *
 * @tparam REAL float and double is supported for static_library
 * @param uni_xyz0
 * @param uni_tex0
 * @param tri_uni
 * @param uni_xyz
 * @param uni_tex
 * @param vtx_xyz
 * @param vtx_tex
 * @param elem_vtx_xyz
 * @param elem_vtx_tex
 */
template <typename REAL>
void UnifySeparateIndexing_PosTex(
    std::vector<REAL> &uni_xyz0,
    std::vector<REAL> &uni_tex0,
    std::vector<unsigned int>& tri_uni,
    std::vector<unsigned int>& uni_xyz,
    std::vector<unsigned int>& uni_tex,
    const std::vector<REAL> &vtx_xyz,
    const std::vector<REAL> &vtx_tex,
    const std::vector<unsigned int> &elem_vtx_xyz,
    const std::vector<unsigned int> &elem_vtx_tex);

// ------------------------------------

DFM2_INLINE void JArray_PSuP_Hair(
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    const std::vector<unsigned int> &aIP_HairRoot);

// ------------------------------------

DFM2_INLINE void MarkConnectedElements(
    std::vector<unsigned int> &elem_flag,
    unsigned int ielem_kernel,
    unsigned int flag,
    const std::vector<unsigned int> &elem_adjelem);

DFM2_INLINE void MakeGroupElem(
    int &num_group,
    std::vector<unsigned int> &elem_flag,
    const std::vector<unsigned int> &elem_vtx,
    const std::vector<unsigned int> &elem_adjelem,
    int num_face_par_elem,
    int num_vtx_par_elem);

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/msh_topology_uniform.cpp"
#endif

#endif /* DFM2_MSHUNI_H */
