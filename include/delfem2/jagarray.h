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
// DONE(2020/12/12): separated mshsubdiv
// TODO: change name mshuni.h
// TODO: separaete jarray.h

#ifndef DFM2_JAGARRAY_H
#define DFM2_JAGARRAY_H

#include <cstdio>
#include <vector>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

// ---------------------------------------------
// function related to jagged array

DFM2_INLINE void JArray_Sort(
    const std::vector<unsigned int> &index,
    std::vector<unsigned int> &array);

DFM2_INLINE void JArray_Sort(
    const unsigned int *index,
    unsigned int size,
    unsigned int *array);

DFM2_INLINE void JArray_AddDiagonal(
    std::vector<unsigned int> &psup_ind1,
    std::vector<unsigned int> &psup1,
    const unsigned int *psup_ind0,
    size_t npsup_ind0,
    const unsigned int *psup0,
    size_t npsup0);

DFM2_INLINE void JArray_Print(
    const std::vector<int> &index,
    const std::vector<int> &array);

/**
 * @details compute 2-ring neighborhood from 1-ring neighborhood
 */
DFM2_INLINE void JArray_Extend(
    std::vector<unsigned int> &psup_ind1,
    std::vector<unsigned int> &psup1,
    const unsigned int *psup_ind0,
    size_t npsup_ind0,
    const unsigned int *psup0);

DFM2_INLINE void JArrayEdgeUnidir_PointSurPoint(
    std::vector<unsigned int> &edge_ind,
    std::vector<unsigned int> &edge,
    //
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup);

} // end namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/jagarray.cpp"
#endif

#endif /* DFM2_JAGARRAY_H */
