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

// DONE(2020/12/09): separate mixed elem
// DONE(2020/12/12): separeted from mshtopo.h


#ifndef DFM2_MSHSUBDIV_H
#define DFM2_MSHSUBDIV_H

#include <cstdio>
#include <vector>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

// -----------------------------------------------

/**
 * @brief making topology for subdivision of quad
 * @details new points is in the order of [old points], [edge points], [face points]
 * @param quad_vtxidx1 (out) new connectivity
 * @param aEdgeFace0 (out) two end points on a edge and two quads touching the edge
 */
DFM2_INLINE void SubdivTopo_MeshQuad(
    std::vector<unsigned int> &quad_vtxidx1,
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    std::vector<unsigned int> &aEdgeFace0,
    const unsigned int *aQuad0,
    size_t nQuad0,
    size_t nPoint0);

DFM2_INLINE void SubdivTopo_MeshHex(
    std::vector<unsigned int> &aHex1,
    std::vector<unsigned int> &psupIndHex0,
    std::vector<unsigned int> &psupHex0,
    std::vector<unsigned int> &aQuadHex0,
    const unsigned int *aHex0,
    size_t nHex0,
    size_t nhp0);

DFM2_INLINE void SubdivTopo_MeshTet(
    std::vector<unsigned int> &aTet1,
    std::vector<unsigned int> &psup_ind,
    std::vector<unsigned int> &psup,
    const unsigned int *aTet0,
    int nTet0,
    unsigned int nPoint0);

DFM2_INLINE unsigned int findEdge(
    unsigned int ip0,
    unsigned int ip1,
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup);

DFM2_INLINE int findFace(
    unsigned int ip0,
    unsigned int ip1,
    unsigned int ip2,
    unsigned int ip3,
    const std::vector<unsigned int> &quad_vtx_idx,
    const std::vector<unsigned int> &elsup_ind,
    const std::vector<unsigned int> &elsup);

DFM2_INLINE void SubdivisionPoints_QuadCatmullClark(
    std::vector<double> &vtx_xyz1,
    //
    const std::vector<unsigned int> &quad_vtxidx1,
    const std::vector<unsigned int> &aEdgeFace0,
    const std::vector<unsigned int> &psupIndQuad0,
    const std::vector<unsigned int> &psupQuad0,
    const unsigned int *quad_vtxidx0,
    size_t num_quad0,
    const double *vtx_xyz0,
    size_t num_vtx0);

DFM2_INLINE void SubdivPoints3_MeshQuad(
    std::vector<double> &aXYZ1,
    //
    const std::vector<int> &aEdgeFace0,
    const std::vector<unsigned int> &aQuad0,
    const std::vector<double> &aXYZ0);

DFM2_INLINE void SubdivisionPoints_Hex(
    std::vector<double> &aXYZ1,
    //
    const std::vector<unsigned int> &psupIndHex0,
    const std::vector<unsigned int> &psupHex0,
    const std::vector<unsigned int> &aQuadHex0,
    const unsigned int *aHex0,
    unsigned int nHex0,
    const double *aXYZ0,
    unsigned int nXYZ0);

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/mshsubdiv.cpp"
#endif

#endif /* DFM2_MSHSUBDIV_H */
