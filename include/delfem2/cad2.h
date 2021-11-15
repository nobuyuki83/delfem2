/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_CAD2_H
#define DFM2_CAD2_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/vec2.h"
#include "delfem2/cadtopo.h"
#include "delfem2/srchbv2aabb.h"

namespace delfem2 {

class CCad2D_VtxGeo {
 public:
  CCad2D_VtxGeo(const CVec2d &p) : pos(p) {}
 public:
  CVec2d pos;
};

/**
 * @details this class should be independent from any other classes except for "CVector2" and "CBoundingBox2D"
 * std::vector<CCad2D_EdgeGeo> will stands for a loop of curves
 */
class CCad2D_EdgeGeo {
 public:
  enum EDGE_TYPE {
    LINE = 0,
    BEZIER_CUBIC = 1,
    BEZIER_QUADRATIC = 2,
  };
  CCad2D_EdgeGeo() {
    type_edge = LINE;
  }
  CCad2D_EdgeGeo(const CCad2D_EdgeGeo &e) {
    this->p0 = e.p0;
    this->p1 = e.p1;
    this->type_edge = e.type_edge;
    this->param = e.param;
  }

  void SetLine(){
    this->type_edge = EDGE_TYPE::LINE;
    this->param.clear();
  }

  void SetQuadraticBezierCurve(const CVec2d& pos0){
    this->type_edge = EDGE_TYPE::BEZIER_QUADRATIC;
    this->param = {
        pos0.x - p0.x,
        pos0.y - p0.y };
  }

  void SetCubicBezierCurve(const CVec2d& pos0, const CVec2d& pos1){
    this->type_edge = EDGE_TYPE::BEZIER_CUBIC;
    this->param = {
        pos0.x - p0.x,
        pos0.y - p0.y,
        pos1.x - p1.x,
        pos1.y - p1.y };
  }

  std::vector<double> GenMesh(unsigned int ndiv) const;

  double Distance(double x, double y) const;

  double LengthMesh() const;

  double LengthNDiv(unsigned int ndiv) const;

  CBoundingBox2<double> BB() const {
    CBoundingBox2<double> bb;
    bb.Add(p0.x, p0.y);
    bb.Add(p1.x, p1.y);
    // for (unsigned int ip = 0; ip < aP.size(); ++ip) { bb.Add(aP[ip].x, aP[ip].y); }
    return bb;
  }

  void Transform(double A[4]){
    this->MatVec2(p0.x,p0.y, A, p0.x, p0.y);
    this->MatVec2(p1.x,p1.y, A, p1.x, p1.y);
    if( type_edge == CCad2D_EdgeGeo::BEZIER_CUBIC ){
      assert( param.size() == 4);
      this->MatVec2(param[0],param[1],A,param[0],param[1]);
      this->MatVec2(param[2],param[3],A, param[2],param[3]);
    }
    if( type_edge == CCad2D_EdgeGeo::BEZIER_QUADRATIC ){
      assert( param.size() == 2);
      this->MatVec2(param[0],param[1],A,param[0],param[1]);
    }
  }

 public:
  CVec2d p0, p1;
  EDGE_TYPE type_edge; // 0: line, 1:Cubic Bezier 2: Quadratic Bezier
  std::vector<double> param;
 private:
  void MatVec2(double& x1, double& y1,
               const double A[4],
               double x0, double y0){
    x1 = A[0]*x0 + A[1]*y0;
    y1 = A[2]*x0 + A[3]*y0;
  }
};

double AreaLoop(
    const std::vector<CCad2D_EdgeGeo> &aEdge);

std::vector<CCad2D_EdgeGeo> InvertLoop(
    const std::vector<CCad2D_EdgeGeo> &aEdge);

std::vector<CCad2D_EdgeGeo> RemoveEdgeWithZeroLength(
    const std::vector<CCad2D_EdgeGeo> &aEdge);

CBoundingBox2<double> BB_LoopEdgeCad2D(
    const std::vector<CCad2D_EdgeGeo> &aEdge);

// ---------------------------------------------------------------------------------

/**
 * @brief class to define 2D shapes bounded by parametric curves
 */
class CCad2D {
 public:
  CCad2D() = default;

  void Clear() {
    aVtx.clear();
    aEdge.clear();
    topo.Clear();
  }
  // --------------------
  // const method here
  std::vector<double> MinMaxXYZ() const;

  CBoundingBox2<double> BB() const;

  bool Check() const;

  size_t nVtx() const { return aVtx.size(); }

  size_t nEdge() const { return aEdge.size(); }

  // ----------------------------------
  // geometric operations from here

  void AddPolygon(
      const std::vector<double> &aXY);

  void AddFace(
      const std::vector<CCad2D_EdgeGeo> &aEdge);

  void AddVtxFace(
      double x0, double y0, unsigned int ifc_add);

  void AddVtxEdge(
      double x, double y, unsigned int ie_add);

  void CopyVertexPositionsToEdges(){
    for(size_t ie=0;ie<topo.edges.size();++ie) {
      const int iv0 = topo.edges[ie].iv0;
      const int iv1 = topo.edges[ie].iv1;
      aEdge[ie].p0 = aVtx[iv0].pos;
      aEdge[ie].p1 = aVtx[iv1].pos;
    }
  }

 public:
  CadTopo topo;
  std::vector<CCad2D_VtxGeo> aVtx;
  std::vector<CCad2D_EdgeGeo> aEdge;
};

// ---------------------------------------

class Cad2_Ui{
 public:
  void Pick(
      double x0, double y0,
      double view_height,
      const CCad2D& cad);

  void DragPicked(
      CCad2D& cad,
      double p1x, double p1y,
      double p0x, double p0y);

 public:
  unsigned int ivtx_picked = UINT_MAX;
  unsigned int iedge_picked = UINT_MAX;
  unsigned int iface_picked = UINT_MAX;
  int ipicked_elem = 0;
  std::array<float, 2> picked_pos{0.f, 0.f};
};

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/cad2.cpp"
#endif

#endif /* DFM2_CAD2_H */
