/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


/**
 * @details this file should not depends on anything except for  std::vector, std::functional
 */

#ifndef DFM2_MSH_NORMAL_H
#define DFM2_MSH_NORMAL_H

#include <vector>

#include "delfem2/dfm2_inline.h"

// -----------------
// work on points

namespace delfem2 {

/**
 * @brief Normal at the vertex of a triangle mesh.
 */
template<typename REAL>
DFM2_INLINE void Normal_MeshTri3D(
    REAL *vtx_normal,
    const REAL *vtx_xyz,
    size_t num_vtx,
    const unsigned int *tri_vtx,
    size_t num_tri);

/**
 * @brief Normal at the vertex of a triangle mesh.
 */
template<typename REAL>
std::array<REAL,3> Normal_TriInMeshTri3(
    unsigned int itri,
    const REAL *vtx_xyz,
    const unsigned int *tri_vtx);

/**
 * @brief Normal at the vertex of a quad mesh. Defined for "float" and "double"
 */
template<typename REAL>
DFM2_INLINE void Normal_MeshQuad3(
    std::vector<REAL> &aNorm,
    const std::vector<REAL> &aXYZ,
    const std::vector<unsigned int> &aQuad);

} // delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/msh_normal.cpp"
#endif

#endif /* DFM2_MSHMISC_H */
