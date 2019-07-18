/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef CAD2D_h
#define CAD2D_h

#include "delfem2/vec2.h"
#include "delfem2/cadtopo.h"
#include "delfem2/dyntri_v2.h"

#include "delfem2/gl_funcs.h"


class CCad2D_VtxGeo{
public:
  CCad2D_VtxGeo(const CVector2& p) : pos(p){}
public:
  CVector2 pos;
};
class CCad2D_EdgeGeo{
public:
  CCad2D_EdgeGeo(){
    type_edge = 0;
  }
  void GenMesh(unsigned int iedge, const CCadTopo& topo,
               std::vector<CCad2D_VtxGeo>& aVtxGeo);
  void GetInternalPoints_ElemLen(std::vector<CVector2>& aV, double elen) const;
  double Distance(double x, double y) const;
  double Length() const;
public:
  CVector2 p0,p1;
  int type_edge; // 0: line, 1:Cubic Bezier
  std::vector<double> param;
  ///
  std::vector<CVector2> aP;
};
class CCad2D_FaceGeo{
public:
  std::vector<unsigned int> aTri;
  std::vector<double> aXY;
public:
  void GenMesh(unsigned int iface0, const CCadTopo& topo, 
               std::vector<CCad2D_EdgeGeo>& aEdgeGeo);
  bool IsInside(double x, double y) const {
    for(unsigned int it=0;it<aTri.size()/3;++it){
      const double q0[2] = {x,y};
      const int i0 = aTri[it*3+0];
      const int i1 = aTri[it*3+1];
      const int i2 = aTri[it*3+2];
      const double* p0 = aXY.data()+i0*2;
      const double* p1 = aXY.data()+i1*2;
      const double* p2 = aXY.data()+i2*2;
      double a0 = TriArea2D(q0, p1, p2);
      double a1 = TriArea2D(p0, q0, p2);
      double a2 = TriArea2D(p0, p1, q0);
      if( a0 > 0 && a1 > 0 && a2 > 0 ){ return true; }
    }
    return false;
  }
};

//////////////////

class CCad2D
{
public:
  CCad2D(){
    std::cout << "CCAD2D -- construct" << std::endl;
    ivtx_picked = -1;
    iedge_picked = -1;
    iface_picked = -1;
    is_draw_face = true;
  }
  void Clear(){
    aVtx.clear();
    aEdge.clear();
    aFace.clear();
    topo.Clear();
  }
  void Tessellation();
  void Pick(double x0, double y0,
            double view_height);
  void DragPicked(double p1x, double p1y, double p0x, double p0y);
  std::vector<double> MinMaxXYZ() const;
  void Check() const;
  void AddPolygon(const std::vector<double>& aXY);
  void AddVtxEdge(double x, double y, int ie_add);
  void Meshing(CMeshDynTri2D& dmesh,
               std::vector<int>& aFlgPnt,
               std::vector<int>& aFlgTri,
               double len) const;
  void GetPointsEdge(std::vector<int>& aIdP,
                     const double* pXY, int np,
                     const std::vector<int>& aIE,
                     double tolerance ) const;
  void SetEdgeType(int iedge, int itype, std::vector<double>& param){
    assert( iedge >= 0 && iedge< (int)aEdge.size() );
    aEdge[iedge].type_edge = itype;
    aEdge[iedge].param = param;
    this->Tessellation();
  }
  int GetEdgeType(int iedge) const {
    assert( iedge >= 0 && iedge< (int)aEdge.size() );
    return aEdge[iedge].type_edge;
  }
  int nFace() const { return aFace.size(); }
  int nVtx() const { return aVtx.size(); }
  int nEdge() const { return aEdge.size(); }
  std::vector<double> XY_VtxCtrl_Face(int iface) const;
  std::vector<std::pair<int,bool> >  Ind_Edge_Face(int iface) const;
  std::vector<int> Ind_Vtx_Face(int iface) const;
  std::vector<int> Ind_Vtx_Edge(int iedge) const;
public:
  CCadTopo topo;
  /////
  std::vector<CCad2D_VtxGeo> aVtx;
  std::vector<CCad2D_EdgeGeo> aEdge;
  std::vector<CCad2D_FaceGeo> aFace;
  
  int ivtx_picked;
  int iedge_picked;
  int iface_picked;
  int ipicked_elem;
  
  bool is_draw_face;
};



#endif
