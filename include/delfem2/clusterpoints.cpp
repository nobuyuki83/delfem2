/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <cmath>
#include <climits>
#include <cassert>
#include <vector>
#include <set>
#include "delfem2/clusterpoints.h"

namespace clusterpoints {

double Dot3(const double p0[3], const double p1[3]) {
  return p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2];
}

double Length3(const double p[3]) {
  return sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
}

template<typename VAL>
VAL myMin(VAL a, VAL b) {
  return (a < b) ? a : b;
}

template<typename VAL>
VAL myMax(VAL a, VAL b) {
  return (a > b) ? a : b;
}

class CData{
public:
  CData(double score0, unsigned int ipsup0): score(score0), ipsup(ipsup0) {}
  bool operator < (const CData& d) const {
    return score < d.score;
  }
public:
  double score;
  unsigned int ipsup;
};

}

// ----------------------------------

unsigned int Find_IndexPoint_From_IndexPsup(
    unsigned int ipsup,
    const unsigned int* psup_ind,
    const unsigned int* psup)
{
  const unsigned int jp = psup[ipsup];
  for(unsigned int jpsup=psup_ind[jp];jpsup<psup_ind[jp+1];++jpsup){
    const unsigned int ip = psup[jpsup];
    if( ipsup < psup_ind[ip] ) continue;
    if( psup_ind[ip+1] <= ipsup ) continue;
    return ip;
  }
  return UINT_MAX;
}

DFM2_INLINE void delfem2::BinaryClustering_Points3d(
    std::vector<double>& aXYZ1,
    std::vector<double>& aArea1,
    std::vector<double>& aNorm1,
    std::vector<unsigned int>& map01,
    //
    const std::vector<double>& aXYZ0,
    const std::vector<double>& aArea0,
    const std::vector<double>& aNorm0,
    const std::vector<unsigned int>& psup_ind0,
    const std::vector<unsigned int>& psup0)
{
  namespace cp = clusterpoints;
  const unsigned int np0 = static_cast<unsigned int>(aXYZ0.size()/3);
  assert( aArea0.size() == np0 );
  assert( aNorm0.size() == np0*3 );
  std::set<clusterpoints::CData> aData;
  for (unsigned int ip0 = 0; ip0 < np0; ++ip0) {
    for (unsigned int ipsup0 = psup_ind0[ip0]; ipsup0 < psup_ind0[ip0 + 1]; ++ipsup0) {
      const unsigned int jp0 = psup0[ipsup0];
      if (jp0 < ip0) { continue; }
      const double ai0 = aArea0[ip0];
      const double aj0 = aArea0[jp0];
//      const double dot0 = cp::Dot3(aNorm0.data() + ip0 * 3, aNorm0.data() + jp0 * 3);
//      const auto nvalencei = (double)(psup_ind0[ip0+1] - psup_ind0[ip0]);
//      const auto nvalencej = (double)(psup_ind0[jp0+1] - psup_ind0[jp0]);
//      const double score0 = dot0 * min(ai0 / aj0, aj0 / ai0)*nvalence;
      const double score0 = cp::myMin(ai0 / aj0, aj0 / ai0);// * cp::max(nvalencei,nvalencej);
//      const double score0 = min(ai0 / aj0, aj0 / ai0) / min(nvalencei/nvalencej, nvalencej/nvalencei);
      aData.insert(clusterpoints::CData(score0, ipsup0));
    }
  }
  {
    aArea1.resize(0);
    aXYZ1.resize(0);
    aNorm1.resize(0);
    const auto np1_guess = (unsigned int)(np0*0.7);
    aArea1.reserve(np1_guess);
    aXYZ1.reserve(np1_guess*3);
    aNorm1.reserve(np1_guess*3);
  }
  map01.assign(np0,UINT_MAX);
  for (const auto &data: aData) {
    const unsigned int ip0 = Find_IndexPoint_From_IndexPsup(
        data.ipsup,
        psup_ind0.data(),
        psup0.data());
    const unsigned int jp0 = psup0[data.ipsup];
    assert( ip0 < jp0 );
    const unsigned int ip1 = static_cast<unsigned int>(aXYZ1.size()/3); // index of new node
    assert( aArea1.size() == ip1 && aNorm1.size() == ip1*3 );
    if( map01[ip0] != UINT_MAX || map01[jp0] != UINT_MAX ) continue;
    map01[ip0] = ip1;
    map01[jp0] = ip1;
    const double ai0 = aArea0[ip0];
    const double aj0 = aArea0[jp0];
    aArea1.push_back(ai0 + aj0);
    { // average coordinates
      aXYZ1.push_back((ai0 * aXYZ0[ip0 * 3 + 0] + aj0 * aXYZ0[jp0 * 3 + 0]) / (ai0 + aj0) );
      aXYZ1.push_back((ai0 * aXYZ0[ip0 * 3 + 1] + aj0 * aXYZ0[jp0 * 3 + 1]) / (ai0 + aj0) );
      aXYZ1.push_back((ai0 * aXYZ0[ip0 * 3 + 2] + aj0 * aXYZ0[jp0 * 3 + 2]) / (ai0 + aj0) );
    }
    { // make new normal by averaging original normals
      const double n1[3] = {
          (ai0 * aNorm0[ip0 * 3 + 0] + aj0 * aNorm0[jp0 * 3 + 0]) / (ai0 + aj0),
          (ai0 * aNorm0[ip0 * 3 + 1] + aj0 * aNorm0[jp0 * 3 + 1]) / (ai0 + aj0),
          (ai0 * aNorm0[ip0 * 3 + 2] + aj0 * aNorm0[jp0 * 3 + 2]) / (ai0 + aj0)};
      const double ln1 = cp::Length3(n1);
      aNorm1.push_back(n1[0]/ln1);
      aNorm1.push_back(n1[1]/ln1);
      aNorm1.push_back(n1[2]/ln1);
    }
  }
  for(unsigned int ip0=0;ip0<np0;++ip0){ // points that are not clustered
    if( map01[ip0] != UINT_MAX ){ continue; }
    const unsigned int ip1 = static_cast<unsigned int>(aXYZ1.size()/3); // index of new node
    assert( aArea1.size() == ip1 && aNorm1.size() == ip1*3 );
    map01[ip0] = ip1;
    aXYZ1.push_back(aXYZ0[ip0*3+0]);
    aXYZ1.push_back(aXYZ0[ip0*3+1]);
    aXYZ1.push_back(aXYZ0[ip0*3+2]);
    aArea1.push_back(aArea0[ip0]);
    aNorm1.push_back(aNorm0[ip0*3+0]);
    aNorm1.push_back(aNorm0[ip0*3+1]);
    aNorm1.push_back(aNorm0[ip0*3+2]);
  }
}

DFM2_INLINE unsigned int
delfem2::BinaryClustering_Points2d(
    unsigned int* map01,
    //
    const unsigned int np0,
    const double* aArea0,
    const unsigned int* psup_ind0,
    const unsigned int* psup0)
{
  for(unsigned int ip=0;ip<np0;++ip){ map01[ip] = UINT_MAX; }
  namespace cp = clusterpoints;
  std::set<clusterpoints::CData> aData;
  for (unsigned int ip0 = 0; ip0 < np0; ++ip0) {
    for (unsigned int ipsup0 = psup_ind0[ip0]; ipsup0 < psup_ind0[ip0 + 1]; ++ipsup0) {
      const unsigned int jp0 = psup0[ipsup0];
      if (jp0 < ip0) { continue; }
      const double ai0 = aArea0[ip0];
      const double aj0 = aArea0[jp0];
      const double score0 = cp::myMin(ai0 / aj0, aj0 / ai0);
      aData.insert(clusterpoints::CData(score0, ipsup0));
    }
  }
  unsigned int np1 = 0;
  for (const auto &data: aData) {
    const unsigned int ip0 = Find_IndexPoint_From_IndexPsup(data.ipsup, psup_ind0, psup0);
    const unsigned int jp0 = psup0[data.ipsup];
    assert( ip0 < jp0 );
    if( map01[ip0] != UINT_MAX || map01[jp0] != UINT_MAX ) continue;
    map01[ip0] = np1;
    map01[jp0] = np1;
    np1++;
  }
  for(unsigned int ip0=0;ip0<np0;++ip0){ // points that are not clustered
    if( map01[ip0] != UINT_MAX ){ continue; }
    map01[ip0] = np1;
    np1++;
  }
  return np1;
}


DFM2_INLINE void
delfem2::Clustering_Psup(
    std::vector<unsigned int>& psup_ind1,
    std::vector<unsigned int>& psup1,
    //
    size_t np1,
    size_t np0,
    const unsigned int* map01,
    const unsigned int* psup_ind0,
    const unsigned int* psup0)
{
  std::vector<unsigned int> map10;
  unsigned int ncnt = 0;
  { // inverse map of map01.
    std::vector<unsigned int> cnt1(np1,0);
    for(unsigned int ip0=0;ip0<np0;++ip0) {
      const unsigned int ip1 = map01[ip0];
      assert(ip1 < np1);
      cnt1[ip1] += 1;
    }
    for(unsigned int ip1=0;ip1<np1;++ip1){
      ncnt = ( cnt1[ip1] > ncnt ) ? cnt1[ip1] : ncnt;
    }
    cnt1.assign(np1,0);
    map10.assign(np1*ncnt, UINT_MAX);
    for(unsigned int ip0=0;ip0<np0;++ip0){
      const unsigned int ip1 = map01[ip0];
      assert( ip1 < np1 );
      unsigned int icnt = cnt1[ip1];
      cnt1[ip1]++;
      assert(icnt<ncnt);
      assert( map10[ip1 * ncnt + icnt ] == UINT_MAX );
      map10[ip1 * ncnt + icnt] = ip0;
    }
  }
  psup_ind1.assign(np1+1,0);
  std::vector<unsigned int> aflag1(np1,UINT_MAX);
  for (unsigned int ip1 = 0; ip1 < np1; ++ip1) {
    aflag1[ip1] = ip1; // avoiding self-connection
    for(unsigned int icnt=0;icnt<ncnt;++icnt) {
      const unsigned int ip0 = map10[ip1 * ncnt + icnt];
      if ( ip0 == UINT_MAX) { break; } // non-clustered point
      assert(map01[ip0] == ip1);
      for (unsigned int ipsup0 = psup_ind0[ip0]; ipsup0 < psup_ind0[ip0 + 1]; ++ipsup0) {
        const unsigned int kp0 = psup0[ipsup0];
        const unsigned int kp1 = map01[kp0];
        if (aflag1[kp1] == ip1) { continue; }
        psup_ind1[ip1 + 1] += 1;
        aflag1[kp1] = ip1;
      }
    }
  }
  for(unsigned int ip1=0;ip1<np1;++ip1){
    psup_ind1[ip1+1] += psup_ind1[ip1];
  }
  const unsigned int npsup1 = psup_ind1[np1];
  psup1.resize(npsup1);
  aflag1.assign(np1,UINT_MAX);
  for (unsigned int ip1 = 0; ip1 < np1; ++ip1) {
    aflag1[ip1] = ip1; // avoiding self-connection
    for(unsigned int icnt=0;icnt<ncnt;++icnt) {
      const unsigned int ip0 = map10[ip1 * ncnt + icnt];
      if( ip0 == UINT_MAX ){ break; } // non-clustered point
      for(unsigned int ipsup0 = psup_ind0[ip0]; ipsup0 < psup_ind0[ip0 + 1]; ++ipsup0){
        const unsigned int kp0 = psup0[ipsup0];
        const unsigned int kp1 = map01[kp0];
        if( aflag1[kp1] == ip1 ){ continue; }
        aflag1[kp1] = ip1;
        psup1[ psup_ind1[ip1] ] = kp1;
        psup_ind1[ip1] += 1;
      }
    }
  }
  for(int ip1=(int)np1-1;ip1>0;--ip1){
    psup_ind1[ip1] = psup_ind1[ip1-1];
  }
  psup_ind1[0] = 0;
  assert( psup_ind1[np1] ==  psup1.size() );
}