/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @brief functions for mesh export/import
 */

#ifndef DFM2_MSHIO_H
#define DFM2_MSHIO_H

#include "delfem2/dfm2_inline.h"
#include <cstdio>
#include <vector>
#include <string>

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

// -----
// PLY

DFM2_INLINE void Write_Ply(
    const std::string &fname,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aTri);

DFM2_INLINE void Write_Ply(
    const std::string &fname,
    unsigned int nXYZ, double *paXYZ,
    unsigned int nTri, unsigned int *paTri);

DFM2_INLINE void Read_Ply(
    const std::string &fname,
    int &nnode_,
    double *&pXYZs_,
    int &ntri_,
    unsigned int *&aTri_);

DFM2_INLINE void Read_Ply(
    const std::string &fname,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri);

// ----------
// Obj


DFM2_INLINE void Write_Obj(
    const std::string &str,
    const double *aXYZ,
    int nXYZ,
    const unsigned int *aTri,
    int nTri);

DFM2_INLINE void Write_Obj_Quad(
    const std::string &str,
    const std::vector<double> &aXYZ,
    const std::vector<int> &aQuad);

/**
 * write obj file for the mesh the elemenet is a jagged array (tris and quads are mixed).
 * @param str
 * @param aXYZ
 * @param aElemInd
 * @param aElem
 */
DFM2_INLINE void Write_Obj_ElemJArray(
    const std::string &str, // mixed elem
    const std::vector<double> &aXYZ,
    const std::vector<int> &aElemInd,
    const std::vector<int> &aElem);


/**
 * to open the obj file with Blender, select the option "Split by Group".
 * @param pathf
 * @param aXYZ
 * @param aTri
 * @param aFlgTri
 */
DFM2_INLINE void Write_Obj_TriFlag(
    const std::string& pathf,
    std::vector<double>& aXYZ,
    std::vector<unsigned int>& aTri,
    std::vector<unsigned int>& aFlgTri);

DFM2_INLINE void Write_Obj(
    const std::string &str,
    const std::vector<std::pair<std::vector<double>, std::vector<int> > > &aMesh);

DFM2_INLINE void Read_Obj(
    const std::string &fname,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri);

/**
 * Read Obj file for quad-only mesh
 * @param fname
 * @param aXYZ
 * @param aQuad
 */
DFM2_INLINE void Read_Obj_MeshQuad3(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aQuad,
    const std::string &fname);

DFM2_INLINE void Read_Obj(
    std::stringstream &ssobj,
    std::vector<double> &aXYZ,
    std::vector<int> &aTri);

DFM2_INLINE void Read_Obj2(
    const std::string &fname,
    std::vector<double> &aXYZ,
    std::vector<int> &aTri);

DFM2_INLINE void Read_Obj3(
    const std::string &fname,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri);

class CTriGroup {
public:
  std::string name_group;
  std::string name_mtl;
  int imtl;
  std::vector<unsigned int> aTriVtx;
  std::vector<unsigned int> aTriNrm;
};

/**
 * Load wavefront Obj file with triangle group
 * @param[in] fname
 * @param[out] fname_mtl
 * @param[out] aXYZ
 * @param[out] aNorm
 * @param[out] aTriGroup
 */
DFM2_INLINE void Load_Obj(
    const std::string &fname,
    std::string &fname_mtl,
    std::vector<double> &aXYZ,
    std::vector<double> &aNorm,
    std::vector<CTriGroup> &aTriGroup);

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
    const int nElem);

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

}

#ifdef DFM2_HEADER_ONLY
#  include "delfem2/mshio.cpp"
#endif

#endif /* meshio_hpp */
