#ifndef DYNTRI_V2_H
#define DYNTRI_V2_H

#include <map>
#include <algorithm>
#include <stack>

#include "delfem2/vec2.h"
#include "delfem2/dyntri.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

class CMeshDensity
{
public:
  virtual double edgeLengthRatio(double px, double py) const = 0;
};

// TODO: there should be three optoions for adding point on edge, 0:none, 1:only between input points, 2: resample everything
bool GenerateTesselation2(std::vector<int>& aTri_out, // out
                          std::vector<double>& aXY_out, // out
                          std::vector<int>& aPtrVtxInd,
                          std::vector<int>& aVtxInd,
                          ////
                          double elen,
                          const CMeshDensity& mesh_density,
                          bool is_uniform_resample_loop, // good for polyline curve in
                          const std::vector< std::vector<double> >& aVecAry0); // in


bool GenerateTesselation2(std::vector<int>& aTri_out, // out
                          std::vector<double>& aXY_out, // out
                          std::vector<int>& aPtrVtxInd,
                          std::vector<int>& aVtxInd,
                          ////
                          double elen,
                          bool is_uniform_resample_loop, // good for polyline curve
                          const std::vector< std::vector<double> >& aVecAry0); // in

bool DelaunayAroundPoint(int ipo0,
                         std::vector<CEPo2>& aPo,
                         std::vector<ETri>& aTri,
                         const std::vector<CVector2>& aVec2);

void DrawMeshDynTri_Edge(const std::vector<ETri>& aSTri,
                         const std::vector<CVector2>& aVec2);

void DrawMeshDynTri_FaceNorm(const std::vector<ETri>& aSTri,
                             const std::vector<CVector2>& aVec2);

class CMeshDynTri2D{
public:
  void Initialize(const double* aXY, int nPo,
                  const int* aTri, int nTri)
  {
    aVec2.resize(nPo);
    for(int ipo=0;ipo<nPo;ipo++){
      aVec2[ipo].x = aXY[ipo*2+0];
      aVec2[ipo].y = aXY[ipo*2+1];
    }
    InitializeMesh(aEPo, aETri,
                   aTri, nTri, nPo);
  }
  void Check()
  {
    CheckTri(aETri);
    CheckTri(aEPo, aETri);
  }
  std::vector<double> MinMax_XYZ() const {
    double x_min,x_max, y_min,y_max;
    x_min=x_max=aVec2[0].x;
    y_min=y_max=aVec2[0].y;
    for(unsigned int ipo=0;ipo<aEPo.size();ipo++){
      const double x = aVec2[ipo].x;
      const double y = aVec2[ipo].y;
      x_min = (x_min < x) ? x_min : x;
      x_max = (x_max > x) ? x_max : x;
      y_min = (y_min < y) ? y_min : y;
      y_max = (y_max > y) ? y_max : y;
    }
    std::vector<double> bb;
    bb.push_back(x_min);
    bb.push_back(x_max);
    bb.push_back(y_min);
    bb.push_back(y_max);
    bb.push_back(0.0);
    bb.push_back(0.0);
    return bb;
  }
  int insertPointElem(int itri0, double r0, double r1){
    CVector2 v2;
    {
      int i0 = aETri[itri0].v[0];
      int i1 = aETri[itri0].v[1];
      int i2 = aETri[itri0].v[2];
      v2 = r0*aVec2[i0]+r1*aVec2[i1]+(1-r0-r1)*aVec2[i2];
    }
    const int ipo0 = aEPo.size();
    aVec2.push_back(v2);
    aEPo.push_back(CEPo2());
    InsertPoint_Elem(ipo0, itri0, aEPo, aETri);
    return ipo0;
  }
  void DelaunayAroundPoint(int ipo){
    ::DelaunayAroundPoint(ipo, aEPo, aETri, aVec2);
  }
  void Draw_FaceNorm()const { DrawMeshDynTri_FaceNorm(aETri,aVec2); }
  void Draw_Edge() const { DrawMeshDynTri_Edge(aETri,aVec2); }
  void draw() const { this->Draw_Edge(); }
  int nTri() const { return aETri.size(); }
  void DeleteTriEdge(int itri, int iedge){ Collapse_ElemEdge(itri, iedge, aEPo, aETri); }
public:
  std::vector<CEPo2> aEPo;
  std::vector<ETri> aETri;
  std::vector<CVector2> aVec2;
};

#endif // #endif SURFACE_MESH_H
