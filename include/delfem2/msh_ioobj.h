/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @brief functions for mesh export/import
 */

#ifndef DFM2_MSH_IOOBJ_H
#define DFM2_MSH_IOOBJ_H

#include <cstdio>
#include <vector>
#include <string>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {


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
    const std::string &pathf,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri,
    std::vector<unsigned int> &aFlgTri);

DFM2_INLINE void Write_Obj(
    const std::string &str,
    const std::vector<std::pair<std::vector<double>, std::vector<unsigned int> > > &aMesh);

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

DFM2_INLINE void Read_Obj2(
    const std::string &fname,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri);

DFM2_INLINE void Read_Obj3(
    const std::string &fname,
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTri);

class TriGroupWavefrontObj {
 public:
  std::string name_group;
  std::string name_mtl;
  int idx_material = -1;
  std::vector<unsigned int> vec_idx_vtx;
  std::vector<unsigned int> vec_idx_tex;
  std::vector<unsigned int> vec_idx_nrm;
};

/**
 * Load wavefront Obj file with triangle group
 * @param[in] fname
 * @param[out] fname_mtl
 * @param[out] vec_xyz
 * @param[out] vec_nrm
 * @param[out] vec_tri_group
 */
DFM2_INLINE void Load_Obj(
    const std::string &fname,
    std::string &fname_mtl,
    std::vector<double> &vec_xyz,
    std::vector<double> &vec_tex,
    std::vector<double> &vec_nrm,
    std::vector<TriGroupWavefrontObj> &vec_tri_group);

class CMaterial{
public:
  std::string name_mtl;
  float Kd[4];
  float Ka[4];
  float Ks[4];
  float Ke[4];
  float Ns;
  int illum;
  std::string map_Kd;
};

void Load_Mtl(const std::string& fname,
              std::vector<CMaterial>& aMtl);

class CMeshMultiElem{
public:
  void ReadObj(const std::string& fname);
//  void Draw() const;
  std::vector<double> AABB3_MinMax() const;
  void ScaleXYZ(double s);
  void TranslateXYZ(double x, double y, double z);
public:
  std::vector<double> aXYZ;
  std::vector<double> aTex;
  std::vector<double> aNorm;
  std::vector<delfem2::TriGroupWavefrontObj> aObjGroupTri;
  std::vector<CMaterial> aMaterial;
};



} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/msh_ioobj.cpp"
#endif

#endif // DFM2_MSH_IOOBJ_H
