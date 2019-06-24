/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DYNTRI_V3_H
#define DYNTRI_V3_H

#include <map>
#include <algorithm>
#include <stack>

#include "delfem2/vec3.h"
#include "delfem2/dyntri.h"

CVector3 normalTri(int itri0,
                   const std::vector<ETri>& aSTri,
                   const std::vector<CVector3>& aVec3);

bool CheckTri(const std::vector<CEPo2>& aPo3D,
              const std::vector<ETri>& aSTri,
              const std::vector<CVector3>& aVec3);

bool FindRayTriangleMeshIntersections(std::vector<CVector3> &intersectionPoints,
                                      const CVector3 &line0,
                                      const CVector3 &line1,
                                      const std::vector<ETri>& aTri,
                                      const std::vector<CVector3>& aVec3);

bool DelaunayAroundPoint(int ipo0,
                         std::vector<CEPo2>& aPo,
                         std::vector<ETri>& aTri,
                         const std::vector<CVector3>& aVec3);

//////////////////////////////////////////////////////////////////////////////////////////////////


void DrawMeshDynTri_FaceNorm(const std::vector<ETri>& aSTri,
                             const std::vector<CVector3>& aVec3);

void DrawMeshDynTri_Edge(const std::vector<ETri>& aSTri,
                         const std::vector<CVector3>& aVec3);


class CMeshDynTri3D{
public:
  void Initialize(const double* aPo, int nPo, int ndim,
                  const unsigned int* aTri, int nTri)
  {
    aVec3.resize(nPo);
    for(int ipo=0;ipo<nPo;ipo++){
      if( ndim == 3 ){
        aVec3[ipo].x = aPo[ipo*3+0];
        aVec3[ipo].y = aPo[ipo*3+1];
        aVec3[ipo].z = aPo[ipo*3+2];
      }
      else if( ndim == 2 ){
        aVec3[ipo].x = aPo[ipo*2+0];
        aVec3[ipo].y = aPo[ipo*2+1];
        aVec3[ipo].z = 0.0;
      }
    }
    InitializeMesh(aEPo, aETri,
                   aTri, nTri, nPo);
  }
  void Check()
  {
    CheckTri(aETri);
    CheckTri(aEPo, aETri);
    CheckTri(aEPo, aETri, aVec3);
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
  int insertPointElem(int itri0, double r0, double r1){
    const int ipo0 = aEPo.size();
    CVector3 v3;
    {
      int i0 = aETri[itri0].v[0];
      int i1 = aETri[itri0].v[1];
      int i2 = aETri[itri0].v[2];
      v3 = r0*aVec3[i0]+r1*aVec3[i1]+(1-r0-r1)*aVec3[i2];
    }
    aVec3.push_back(v3);
    aEPo.push_back(CEPo2());
    InsertPoint_Elem(ipo0, itri0, aEPo, aETri);
    return ipo0;
  }
  void DelaunayAroundPoint(int ipo){
    ::DelaunayAroundPoint(ipo, aEPo, aETri, aVec3);
  }
  void Draw_FaceNorm()const { DrawMeshDynTri_FaceNorm(aETri,aVec3); }
  void Draw_Edge() const { DrawMeshDynTri_Edge(aETri,aVec3); }
  void draw() const { this->Draw_Edge(); }
  int nTri() const { return aETri.size(); }
  int nPoint() const { return aEPo.size(); }
  void DeleteTriEdge(int itri, int iedge){ Collapse_ElemEdge(itri, iedge, aEPo, aETri); }
public:
  std::vector<CEPo2> aEPo;
  std::vector<ETri> aETri;
  std::vector<CVector3> aVec3;
};


#endif // #endif SURFACE_MESH_H
