/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/mshtopoio.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "delfem2/msh_iomisc.h"
#include "delfem2/msh_ioobj.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshprimitive.h"
#include "delfem2/dtri3_v3dtri.h"

// ----------------------------------------------

// probably std::stroi is safer to use but it is only for C++11
static int myStoi(const std::string &str) {
  char *e;
  long d = std::strtol(str.c_str(), &e, 0);
  return (int) d;
}

static double myStof(const std::string &str) {
  char *e;
  float fval = std::strtof(str.c_str(), &e);
  return fval;
}

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
  dfm2::JArrayPointSurPoint_MeshOneRingNeighborhood(psup_ind, psup,
                                                    aTri0.data(),
                                                    elsup_ind, elsup,
                                                    3, np0);
  //  std::cout << "psup" << std::endl;
  //  Print_IndexedArray(psup_ind, psup);
  // ---------
  std::vector<unsigned int> edge_ind, edge;
  dfm2::JArrayEdgeUnidir_PointSurPoint(edge_ind, edge,
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
    int iedge01, iedge12, iedge20;
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

// ----------------------

void CMeshMultiElem::ReadObj(const std::string &path_obj) {
  std::string fname_mtl;
  Load_Obj(path_obj,
           fname_mtl, aXYZ, aTex, aNorm, aObjGroupTri);
  std::string path_dir = std::string(path_obj.begin(), path_obj.begin() + path_obj.rfind("/"));
  Load_Mtl(path_dir + "/" + fname_mtl,
           aMaterial);
//  std::cout << aObjGroupTri.size() << " " << aMaterial.size() << std::endl;
  { //
    std::map<std::string, int> mapMtlName2Ind;
    for (int imtl = 0; imtl < (int) aMaterial.size(); ++imtl) {
      mapMtlName2Ind.insert(std::make_pair(aMaterial[imtl].name_mtl, imtl));
    }
    for (auto &iogt : aObjGroupTri) {
      std::string name_mtl = iogt.name_mtl;
      auto itr = mapMtlName2Ind.find(name_mtl);
      if (name_mtl.empty() || itr == mapMtlName2Ind.end()) {
        iogt.idx_material = -1;
        continue;
      }
      iogt.idx_material = itr->second;
    }
  }
}

std::vector<double> CMeshMultiElem::AABB3_MinMax() const {
  double c[3], w[3];
  delfem2::CenterWidth_Points3(c, w,
                               aXYZ);
  std::vector<double> aabb(6);
  aabb[0] = c[0] - 0.5 * w[0];
  aabb[1] = c[0] + 0.5 * w[0];
  aabb[2] = c[1] - 0.5 * w[1];
  aabb[3] = c[1] + 0.5 * w[1];
  aabb[4] = c[2] - 0.5 * w[2];
  aabb[5] = c[2] + 0.5 * w[2];
  return aabb;
}

void CMeshMultiElem::ScaleXYZ(double s) {
  delfem2::Scale_PointsX(aXYZ,
                         s);
}

void CMeshMultiElem::TranslateXYZ(double x, double y, double z) {
  delfem2::Translate_Points3(aXYZ,
                             x, y, z);
}





