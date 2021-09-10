/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @brief functions for mesh export/import
 */

#ifndef DFM2_MSH_IO_MISC_H
#define DFM2_MSH_IO_MISC_H

#include <cstdio>
#include <vector>
#include <string>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

// ---------------
// IO functions

DFM2_INLINE void Write_MeshTri3D(
    std::ofstream &fout,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri,
    int ndim = 3);

DFM2_INLINE void Read_MeshTri3D(
    std::ifstream &fin,
    std::vector<double> &aXYZ,
    std::vector<int> &aTri);

// -------
// TetGen

DFM2_INLINE void Read_MeshTet3D_TetGen(
    const std::string &fname,
    std::vector<double> &aXYZ,
    std::vector<int> &aTet,
    std::vector<int> &aTri);

// -------
// STL

DFM2_INLINE void Write_STL(
    const std::string &str,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri);

// -----------
// VTK

DFM2_INLINE void WriteVTK_Points(
    std::ofstream &fout,
    const std::string &name,
    const double *pXYZ,
    int nXYZ,
    int nDim);


// 5 : VTK_TRIANGLE
// 9 : VTK_QUAD
// 10: VTK_TETRA
// 12: VTK_HEXAHEDRON
// 13: VTK_WEDGE
// 14: VTK_PYRAMD
DFM2_INLINE void WriteVTK_Cells(
    std::ofstream &fout,
    int vtk_elem_type,
    const int *aElem,
    int nElem);

DFM2_INLINE void WriteVTK_Cells(
    std::ofstream &fout,
    const std::vector<int> &aTet,
    const std::vector<int> &aPyrm,
    const std::vector<int> &aPrsm);

DFM2_INLINE void WriteVTK_Data_PointVec(
    std::ofstream &fout,
    const double *aVal,
    int np,
    int nStrideVal,
    int ndim);

DFM2_INLINE void WriteVTK_Data_PointScalar(
    std::ofstream &fout,
    const double *aVal,
    int np,
    int nStrideVal = 1);

DFM2_INLINE void WriteVTK_MapTriScalar(
    const std::string &fpath,
    const std::string &name,
    const std::vector<double> &aXYZ,
    const std::vector<int> &map,
    const std::vector<int> &aTri,
    const std::vector<double> &aVal,
    int nStrideVal, int nOffset);

DFM2_INLINE void ReadVTK(
    std::vector<double> &aXYZ,
    int &type_elem,
    std::vector<int> &aElem,
    std::vector<double> &aPointVal,
    const std::string &fpath);

DFM2_INLINE void Read_MeshTri3D_Nas(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    const char *path);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/msh_io_misc.cpp"
#endif

#endif // DFM2_MSH_IOMISC_H
