/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_CAD2_MESHER_H
#define DFM2_CAD2_MESHER_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/vec2.h"
#include "delfem2/cadtopo.h"
#include "delfem2/cad2.h"
#include "delfem2/srchbv2aabb.h"

namespace delfem2 {

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
  void Meshing(
      CMeshDynTri2D &dmesh,
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

  /**
   * map point to the index to vertex, edge, and face
   */
  std::vector<unsigned int> aFlgPnt;

  /**
   * @brief map triangle index to cad face index
   * @details after calling "this->Meshing()", the size of "this->aFlgTri" should be equal to the number of all the triangles
   */
  std::vector<unsigned int> aFlgTri;

  std::vector< std::vector<std::pair<unsigned int, double> > > edge_point;
};

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/cad2_mesher.cpp"
#endif

#endif /* DFM2_CAD2_MESHER_H */
