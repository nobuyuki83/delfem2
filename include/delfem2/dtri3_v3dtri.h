/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_DTRIV3_H
#define DFM2_DTRIV3_H

#include <map>
#include <algorithm>
#include <stack>

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec3.h"
#include "delfem2/dtri.h"

namespace delfem2 {

DFM2_INLINE CVec3d UnitNormal_DTri3(
    unsigned int itri0,
    const std::vector<CDynTri>& aSTri,
    const std::vector<CVec3d>& aVec3);

DFM2_INLINE bool AssertMeshDTri2(
    const std::vector<CDynPntSur>& aPo3D,
    const std::vector<CDynTri>& aSTri,
    const std::vector<CVec3d>& aVec3);

DFM2_INLINE bool FindRayTriangleMeshIntersections(
    std::vector<CVec3d> &intersectionPoints,
    const CVec3d &line0,
    const CVec3d &line1,
    const std::vector<CDynTri>& aTri,
    const std::vector<CVec3d>& aVec3);

DFM2_INLINE bool DelaunayAroundPoint(
    const unsigned int ipo0,
    std::vector<CDynPntSur>& aPo,
    std::vector<CDynTri>& aTri,
    const std::vector<CVec3d>& aVec3);

// -------------------------------------------------------

class CMeshDynTri3D
{
public:
  void Initialize(const double* aPo, int nPo, int ndim,
                  const unsigned int* aTri, int nTri)
  {
    aVec3.resize(nPo);
    for(int ipo=0;ipo<nPo;ipo++){
      if( ndim == 3 ){
        aVec3[ipo].p[0] = aPo[ipo*3+0];
        aVec3[ipo].p[1] = aPo[ipo*3+1];
        aVec3[ipo].p[2] = aPo[ipo*3+2];
      }
      else if( ndim == 2 ){
        aVec3[ipo].p[0] = aPo[ipo*2+0];
        aVec3[ipo].p[1] = aPo[ipo*2+1];
        aVec3[ipo].p[2] = 0.0;
      }
    }
    InitializeMesh(aEPo, aETri,
                   aTri, nTri, nPo);
  }
  void Check()
  {
    AssertDTri(aETri);
    AssertMeshDTri(aEPo, aETri);
    AssertMeshDTri2(aEPo, aETri, aVec3);
  }
  std::vector<double> MinMax_XYZ() const {
    double x_min,x_max, y_min,y_max, z_min,z_max;
    x_min=x_max=aVec3[0].x;
    y_min=y_max=aVec3[0].y;
    z_min=z_max=aVec3[0].z;
    for(unsigned int ipo=0;ipo<aEPo.size();ipo++){
      const double x = aVec3[ipo].x;
      const double y = aVec3[ipo].y;
      const double z = aVec3[ipo].z;
      x_min = (x_min < x) ? x_min : x;
      x_max = (x_max > x) ? x_max : x;
      y_min = (y_min < y) ? y_min : y;
      y_max = (y_max > y) ? y_max : y;
      z_min = (z_min < z) ? z_min : z;
      z_max = (z_max > z) ? z_max : z;
    }
    std::vector<double> bb;
    bb.push_back(x_min);
    bb.push_back(x_max);
    bb.push_back(y_min);
    bb.push_back(y_max);
    bb.push_back(z_min);
    bb.push_back(z_max);
    return bb;
  }
  unsigned int insertPointElem(int itri0, double r0, double r1){
    const unsigned int ipo0 = static_cast<unsigned int>(aEPo.size());
    CVec3d v3;
    {
      int i0 = aETri[itri0].v[0];
      int i1 = aETri[itri0].v[1];
      int i2 = aETri[itri0].v[2];
      v3 = r0*aVec3[i0]+r1*aVec3[i1]+(1-r0-r1)*aVec3[i2];
    }
    aVec3.push_back(v3);
    aEPo.push_back(CDynPntSur());
    InsertPoint_Elem(ipo0, itri0, aEPo, aETri);
    return ipo0;
  }
  void DelaunayAroundPoint(int ipo){
    delfem2::DelaunayAroundPoint(ipo, aEPo, aETri, aVec3);
  }
  size_t nTri() const { return aETri.size(); }
  size_t nPoint() const { return aEPo.size(); }
  void DeleteTriEdge(int itri, int iedge){ CollapseEdge_MeshDTri(itri, iedge, aEPo, aETri); }
public:
  std::vector<CDynPntSur> aEPo;
  std::vector<CDynTri> aETri;
  std::vector<CVec3d> aVec3;
};
  
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/dtri3_v3dtri.cpp"
#endif

#endif // #endif SURFACE_MESH_H
