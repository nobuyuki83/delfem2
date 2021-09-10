/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mshtopoio.h"

#include <string>
#include <vector>

#include "delfem2/msh_io_misc.h"
#include "delfem2/msh_io_obj.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/dtri3_v3dtri.h"

// --------------------------------------------------

void MeshTri3D_GeodesicPolyhedron(
    std::vector<double> &aXYZ1,
    std::vector<unsigned int> &aTri1) {
  namespace dfm2 = delfem2;
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri0;
  dfm2::MeshTri3D_Icosahedron(aXYZ0, aTri0);
  // -------
  const unsigned int np0 = aXYZ0.size() / 3;
  std::vector<unsigned int> elsup_ind, elsup;
  delfem2::JArray_ElSuP_MeshElem(
      elsup_ind, elsup,
      aTri0.data(), aTri0.size() / 3, 3, np0);
  // -------
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArrayPointSurPoint_MeshOneRingNeighborhood(
      psup_ind, psup,
      aTri0.data(),
      elsup_ind, elsup,
      3, np0);
  //  std::cout << "psup" << std::endl;
  //  Print_IndexedArray(psup_ind, psup);
  // ---------
  std::vector<unsigned int> edge_ind, edge;
  dfm2::JArrayEdgeUnidir_PointSurPoint(
      edge_ind, edge,
      psup_ind, psup);
  //  std::cout << "edge" << std::endl;
  //  Print_IndexedArray(edge_ind, edge);
  // ------------
  double r0 = sqrt((5 + sqrt(5)) * 0.5);
  aXYZ1 = aXYZ0;
  for (unsigned int ip = 0; ip < np0; ++ip) {
    for (unsigned int iedge = edge_ind[ip]; iedge < edge_ind[ip + 1]; ++iedge) {
      const unsigned int ip0 = edge[iedge];
      const double x1 = (aXYZ1[ip * 3 + 0] + aXYZ1[ip0 * 3 + 0]) * 0.5;
      const double y1 = (aXYZ1[ip * 3 + 1] + aXYZ1[ip0 * 3 + 1]) * 0.5;
      const double z1 = (aXYZ1[ip * 3 + 2] + aXYZ1[ip0 * 3 + 2]) * 0.5;
      double mag = r0 / sqrt(x1 * x1 + y1 * y1 + z1 * z1);
      aXYZ1.push_back(x1 * mag);
      aXYZ1.push_back(y1 * mag);
      aXYZ1.push_back(z1 * mag);
    }
  }
  aTri1.clear();
  aTri1.reserve(aTri0.size() * 3);
  for (unsigned int itri = 0; itri < aTri0.size() / 3; ++itri) {
    const unsigned int ip0 = aTri0[itri * 3 + 0];
    const unsigned int ip1 = aTri0[itri * 3 + 1];
    const unsigned int ip2 = aTri0[itri * 3 + 2];
    unsigned int iedge01, iedge12, iedge20;
    {
      if (ip0 < ip1) { iedge01 = dfm2::findEdge(ip0, ip1, edge_ind, edge); }
      else { iedge01 = dfm2::findEdge(ip1, ip0, edge_ind, edge); }
      if (ip1 < ip2) { iedge12 = dfm2::findEdge(ip1, ip2, edge_ind, edge); }
      else { iedge12 = dfm2::findEdge(ip2, ip1, edge_ind, edge); }
      if (ip2 < ip0) { iedge20 = dfm2::findEdge(ip2, ip0, edge_ind, edge); }
      else { iedge20 = dfm2::findEdge(ip0, ip2, edge_ind, edge); }
    }
    aTri1.push_back(ip0);
    aTri1.push_back(iedge01 + np0);
    aTri1.push_back(iedge20 + np0);
    aTri1.push_back(ip1);
    aTri1.push_back(iedge12 + np0);
    aTri1.push_back(iedge01 + np0);
    aTri1.push_back(ip2);
    aTri1.push_back(iedge20 + np0);
    aTri1.push_back(iedge12 + np0);
    aTri1.push_back(iedge01 + np0);
    aTri1.push_back(iedge12 + np0);
    aTri1.push_back(iedge20 + np0);
  }
}
