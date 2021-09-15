/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_CAD2_DTRI_H
#define DFM2_CAD2_DTRI_H

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
    ip0 = -1;
  }
  CCad2D_EdgeGeo(const CCad2D_EdgeGeo &e) {
    this->p0 = e.p0;
    this->p1 = e.p1;
    this->type_edge = e.type_edge;
    this->param = e.param;
    this->aP = e.aP;
    this->ip0 = e.ip0;
  }
  void GenMeshNDiv(unsigned int ndiv);
  void GenMeshLength(double elen);
  double Distance(double x, double y) const;
  double LengthMesh() const;
  double LengthNDiv(unsigned int ndiv) const;
  CBoundingBox2<double> BB() const {
    CBoundingBox2<double> bb;
    bb.Add(p0.x, p0.y);
    bb.Add(p1.x, p1.y);
    for (unsigned int ip = 0; ip < aP.size(); ++ip) { bb.Add(aP[ip].x, aP[ip].y); }
    return bb;
  }
 public:
  CVec2d p0, p1;
  EDGE_TYPE type_edge; // 0: line, 1:Cubic Bezier 2: Quadratic Bezier
  std::vector<double> param;
  //
  std::vector<CVec2d> aP;
  int ip0; //!< ip0 is the p0's point index when mesh is generated
};

double AreaLoop(
    const std::vector<CCad2D_EdgeGeo> &aEdge);

std::vector<CCad2D_EdgeGeo> InvertLoop(
    const std::vector<CCad2D_EdgeGeo> &aEdge);

std::vector<CCad2D_EdgeGeo> RemoveEdgeWithZeroLength(
    const std::vector<CCad2D_EdgeGeo> &aEdge);

void Transform_LoopEdgeCad2D(
    std::vector<CCad2D_EdgeGeo> &aEdge,
    bool is_flip_holizontal,
    bool is_flip_vertical,
    double scale_x,
    double scale_y);

CBoundingBox2<double> BB_LoopEdgeCad2D(
    const std::vector<CCad2D_EdgeGeo> &aEdge);

/**
 * @details this class should be independent from everything
 */
class CCad2D_FaceGeo {
 public:
  std::vector<unsigned int> aTri;
 public:
  bool IsInside(
      double x, double y,
      const std::vector<CVec2d> &aVec2) const {
    for (unsigned int it = 0; it < aTri.size() / 3; ++it) {
      const CVec2d q0(x, y);
      const int i0 = aTri[it * 3 + 0];
      const int i1 = aTri[it * 3 + 1];
      const int i2 = aTri[it * 3 + 2];
      const CVec2d &p0 = aVec2[i0];
      const CVec2d &p1 = aVec2[i1];
      const CVec2d &p2 = aVec2[i2];
      double a0 = Area_Tri(q0, p1, p2);
      double a1 = Area_Tri(p0, q0, p2);
      double a2 = Area_Tri(p0, p1, q0);
      if (a0 > 0 && a1 > 0 && a2 > 0) { return true; }
    }
    return false;
  }
};


// ---------------------------------------------------------------------------------

/**
 * @brief class to define 2D shapes bounded by parametric curves
 */
class CCad2D {
 public:
  CCad2D() {
//    std::cout << "CCAD2D -- construct" << std::endl;
    ivtx_picked = -1;
    iedge_picked = -1;
    iface_picked = -1;
    is_draw_face = true;
  }
  void Clear() {
    aVtx.clear();
    aEdge.clear();
    aFace.clear();
    topo.Clear();
  }
  void Tessellation();
  void Pick(double x0, double y0,
            double view_height);
  void DragPicked(double p1x, double p1y, double p0x, double p0y);
  // --------------------
  // const method here
  std::vector<double> MinMaxXYZ() const;
  CBoundingBox2<double> BB() const;
  bool Check() const;
  int GetEdgeType(int iedge) const {
    assert(iedge >= 0 && iedge < (int) aEdge.size());
    return aEdge[iedge].type_edge;
  }
  size_t nFace() const { return aFace.size(); }
  size_t nVtx() const { return aVtx.size(); }
  size_t nEdge() const { return aEdge.size(); }
  /**
   * @brief return std::vector of XY that bounds the face with index iface
   */
  std::vector<double> XY_VtxCtrl_Face(int iface) const;
  std::vector<double> XY_Vtx(int ivtx) const;
  std::vector<std::pair<int, bool> > Ind_Edge_Face(int iface) const;
  std::vector<int> Ind_Vtx_Face(int iface) const;
  std::vector<int> Ind_Vtx_Edge(int iedge) const;
  /**
   * @brief add index to aIdP if a point in aXY is on the edge
   */
  void GetPointsEdge(std::vector<int> &aIdP,
                     const double *pXY, int np,
                     const std::vector<int> &aIE,
                     double tolerance) const;
  void MeshForVisualization(std::vector<float> &aXYf,
                            std::vector<std::vector<unsigned int> > &aaLine,
                            std::vector<unsigned int> &aTri) const;

  // ----------------------------------
  // geometric operations from here
  void AddPolygon(const std::vector<double> &aXY);
  void AddFace(const std::vector<CCad2D_EdgeGeo> &aEdge);
  void AddVtxFace(double x0, double y0, unsigned int ifc_add);
  void AddVtxEdge(double x, double y, unsigned int ie_add);
  void SetEdgeType(int iedge, CCad2D_EdgeGeo::EDGE_TYPE itype, std::vector<double> &param) {
    assert(iedge >= 0 && iedge < (int) aEdge.size());
    aEdge[iedge].type_edge = itype;
    aEdge[iedge].param = param;
    this->Tessellation();
  }
 public:
  CCadTopo topo;
  //
  std::vector<CCad2D_VtxGeo> aVtx;
  std::vector<CCad2D_EdgeGeo> aEdge;
  std::vector<CCad2D_FaceGeo> aFace;
  std::vector<CVec2d> aVec2_Tessellation;

  int ivtx_picked;
  int iedge_picked;
  int iface_picked;
  int ipicked_elem;

  bool is_draw_face;
};

/**
 * @brief mesher for 2 dimensional CAD
 */
class CMesher_Cad2D {
 public:
  CMesher_Cad2D() {
    edge_length = 0.1;
    nvtx = 0;
    nedge = 0;
    nface = 0;
  }
  void Meshing(CMeshDynTri2D &dmesh,
               const CCad2D &cad2d);
  std::vector<unsigned int> IndPoint_IndEdgeArray(
      const std::vector<int> &aIndEd,
      const CCad2D &cad2d);
  std::vector<int> IndPoint_IndFaceArray(
      const std::vector<int> &aIndFc,
      const CCad2D &cad2d);
  std::vector<unsigned int> IndPoint_IndEdge(
      const unsigned int ie,
      bool is_end_point,
      const CCad2D &cad2d);
 public:
  // inputs for meshing
  double edge_length;
  /**
   * @brief specifiation of how many divisions in the cad edge.
   * @details this specification has more priority than the this->edge_length
   */
  std::map<unsigned int, unsigned int> mapIdEd_NDiv;

  // --------------
  // output for meshing

  size_t nvtx;
  size_t nedge;
  size_t nface;
  std::vector<int> aFlgPnt;

  /**
   * @brief map triangle index to cad face index
   * @details after calling "this->Meshing()", the size of "this->aFlgTri" should be equal to the number of all the triangles
   */
  std::vector<unsigned int> aFlgTri;
};

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/cad2_dtri2.cpp"
#endif

#endif /* DFM2_CAD2_DTRI_H */
