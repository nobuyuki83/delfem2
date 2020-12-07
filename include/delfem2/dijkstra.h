/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * @file geodesic & clustering on mesh using dijkstra method
 */

#ifndef DFM2_DIJKSTRA_H
#define DFM2_DIJKSTRA_H

#include "delfem2/dfm2_inline.h"
#include <queue>
#include <climits>
#include <random>

namespace delfem2 {

namespace dijkstra {

double Distance3(
    const double p0[3], const double p1[3])
{
  return sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]) );
}

template <typename DISTANCE>
class CNode
{
public:
  CNode(unsigned int ino, DISTANCE mesh_dist)
  : ind(ino), dist(mesh_dist){}
  bool operator < (const CNode& lhs) const {
    return this->dist > lhs.dist;
  }
public:
  unsigned int ind;
  DISTANCE dist;
};

}


void DijkstraElem_MeshElem(
    std::vector<unsigned int> &aDist,
    std::vector<unsigned int>& aOrder,
    //
    unsigned int ielm_ker,
    const std::vector<unsigned int> &aElSuEl)
{
  const unsigned int nelem = aDist.size();
  aOrder.assign(nelem,UINT_MAX);
  aDist.assign(nelem, UINT_MAX);
  aDist[ielm_ker] = 0;
  const unsigned int nedelm = aElSuEl.size() / nelem;
  std::priority_queue<dijkstra::CNode<unsigned int>> que;
  que.push(dijkstra::CNode<unsigned int>(ielm_ker, 0));
  unsigned int icnt = 0;
  while (!que.empty()) {
    const unsigned int ielm0 = que.top().ind;
    const unsigned int idist0 = que.top().dist;
    que.pop();
    if( aOrder[ielm0] != UINT_MAX ){ continue; } // already fixed so this is not the shortest path
    // found shortest path
    aOrder[ielm0] = icnt; icnt++;
    for (unsigned int iedelm = 0; iedelm < nedelm; ++iedelm) {
      const unsigned int ielm1 = aElSuEl[ielm0 * nedelm + iedelm];
      if (ielm1 == UINT_MAX) { continue; }
      const unsigned int idist1 = idist0+1;
      if (idist1 >= aDist[ielm1]) { continue; }
      // Found the shortest path ever examined.
      aDist[ielm1] = idist1;
      // put this in the que because this is a candidate for the shortest path
      que.push(dijkstra::CNode<unsigned int>(ielm1, idist1));
    }
  }
  assert(icnt==nelem);
}

void MeshClustering(
    std::vector<unsigned int> &aFlgElm,
    //
    unsigned int ncluster,
    const std::vector<unsigned int> &aTriSuTri,
    unsigned int ntri)
{
  std::vector<unsigned int> aDist0(ntri, UINT_MAX);
  std::random_device rd;
  std::mt19937 rdeng(rd());
  std::uniform_int_distribution<unsigned int> dist0(0, ntri - 1);
  const unsigned int itri_ker = dist0(rdeng);
  assert(itri_ker < ntri);
  aFlgElm.assign(ntri, 0);
  std::vector<unsigned int> aOrder;
  DijkstraElem_MeshElem(
      aDist0,aOrder,
      itri_ker, aTriSuTri);
  for (unsigned int icluster = 1; icluster < ncluster; ++icluster) {
    unsigned int itri_maxdist;
    { // find triangle with maximum distance
      double idist_max = 0;
      for (unsigned int it = 0; it < ntri; ++it) {
        if (aDist0[it] <= idist_max) { continue; }
        idist_max = aDist0[it];
        itri_maxdist = it;
      }
    }
    std::vector<unsigned int> aDist1(ntri, UINT_MAX);
    DijkstraElem_MeshElem(
        aDist1,aOrder,
        itri_maxdist, aTriSuTri);
    for (unsigned int it = 0; it < ntri; ++it) {
      if (aDist1[it] < aDist0[it]) {
        aDist0[it] = aDist1[it];
        aFlgElm[it] = icluster;
      }
    }
  }
}

template <typename PROC>
void DijkstraPoint_MeshTri3D(
    std::vector<double>& aDist,
    std::vector<unsigned int>& aOrder,
    PROC& proc,
    //
    unsigned int ip_ker,
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& psup_ind,
    const std::vector<unsigned int>& psup)
{
  const unsigned int np = aXYZ.size()/3;
  aOrder.assign(np,UINT_MAX);
  aDist.assign(np,-1.);
  std::priority_queue<dijkstra::CNode<double>> que;
  unsigned int icnt = 0;
  que.push(dijkstra::CNode<double>(ip_ker, 0.0));
  while(!que.empty()) {
    const unsigned int ip0 = que.top().ind; assert(ip0<np);
    const double dist0 = que.top().dist;
    que.pop();
    if( aOrder[ip0] != UINT_MAX ){ continue; } // already fixed (this wes necessary)
    aOrder[ip0] = icnt;
    proc.AddPoint(ip0,icnt,aOrder);
    icnt++;
    for(unsigned int ipsup=psup_ind[ip0];ipsup<psup_ind[ip0+1];++ipsup) {
      const unsigned int ip1 = psup[ipsup];
      const double len01 = dijkstra::Distance3(aXYZ.data()+ip0*3,aXYZ.data()+ip1*3);
      const double dist1 = dist0 + len01;
      if( aDist[ip1] >= -0.1 && aDist[ip1] < dist1 ){ continue; }
      aDist[ip1] = dist1; // shortest distance so far
      que.push(dijkstra::CNode<double>(ip1,dist1)); // add candidate for fix
    }
  }
  assert(icnt==np);
}


}

#endif
