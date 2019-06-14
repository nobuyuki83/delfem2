/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef CAD2D_h
#define CAD2D_h

#include "delfem2/vec2.h"
#include "delfem2/funcs_gl.h"
#include "delfem2/cadtopo.h"


class CCad2D_VtxGeo{
public:
  CCad2D_VtxGeo(const CVector2& p) : pos(p){}
public:
  CVector2 pos;
};
class CCad2D_EdgeGeo{
public:
  void GenMesh(unsigned int iedge, const CCadTopo& topo,
               std::vector<CCad2D_VtxGeo>& aVtxGeo);
  double Distance(double x, double y) const;
public:
  CVector2 p0,p1;
  std::vector<CVector2> aP;
};
class CCad2D_FaceGeo{
public:
  //    std::vector<CVector2> aP;
  std::vector<int> aTri;
  std::vector<double> aXY;
public:
  void GenMesh(unsigned int iface0, const CCadTopo& topo, 
               std::vector<CCad2D_EdgeGeo>& aEdgeGeo);
};

//////////////////

class CCad2D
{
public:
  CCad2D(){
    std::cout << "CCAD2D -- construct" << std::endl;
    ivtx_picked = -1;
    iedge_picked = -1;
    is_draw_face = true;
  }
  void Clear(){
    aVtx.clear();
    aEdge.clear();
    aFace.clear();
    topo.Clear();
  }
  void Draw() const;
  void Tessellation();
  void Pick(double x0, double y0,
            double view_height);
  void DragPicked(double p1x, double p1y, double p0x, double p0y);
  std::vector<double> MinMaxXYZ() const;
  void Check() const;
  void AddPolygon(const std::vector<double>& aXY);
  void AddPointEdge(double x, double y, int ie_add);
  void Meshing(std::vector<double>& aXY,
               std::vector<int>& aTri,
               double len) const;
  void GetPointsEdge(std::vector<int>& aIdP,
                     const double* pXY, int np,
                     const std::vector<int>& aIE,
                     double tolerance ) const;
  std::vector<double> GetVertexXY_Face(int iface) const;
public:
  CCadTopo topo;
  /////
  std::vector<CCad2D_VtxGeo> aVtx;
  std::vector<CCad2D_EdgeGeo> aEdge;
  std::vector<CCad2D_FaceGeo> aFace;
  int ivtx_picked;
  int iedge_picked;
  
  bool is_draw_face;
};



#endif
