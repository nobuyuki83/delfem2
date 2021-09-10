/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/jagarray.h"

#include <vector>
#include <cassert>
#include <iostream>
#include <climits>

// in the edge ip -> jp, it holds (ip < jp)
DFM2_INLINE void delfem2::JArrayEdgeUnidir_PointSurPoint(
    std::vector<unsigned int> &edge_ind,
    std::vector<unsigned int> &edge,
    //
    const std::vector<unsigned int> &psup_ind,
    const std::vector<unsigned int> &psup) {
  assert(psup_ind.size() >= 2);
  const std::size_t np = psup_ind.size() - 1;
  edge_ind.resize(np + 1);
  edge_ind[0] = 0;
  edge.clear();
  //
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      unsigned int ip0 = psup[ipsup];
      if (ip0 <= ip) continue;
      edge_ind[ip + 1]++;
    }
  }
  for (unsigned int ip = 0; ip < np; ip++) {
    edge_ind[ip + 1] += edge_ind[ip];
  }
  const unsigned int nedge = edge_ind[np];
  edge.resize(nedge);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      const unsigned int ip0 = psup[ipsup];
      if (ip0 <= ip) continue;
      const unsigned int iedge = edge_ind[ip];
      edge[iedge] = ip0;
      edge_ind[ip]++;
    }
  }
  for (int ip = (int) np; ip > 0; ip--) {
    edge_ind[ip] = edge_ind[ip - 1];
  }
  edge_ind[0] = 0;
}


// ---------------------------------------------

DFM2_INLINE void delfem2::JArray_Print(
    const std::vector<int> &index,
    const std::vector<int> &array) {
  assert(index.size() >= 2);
  const std::size_t np = index.size() - 1;
  for (unsigned int ip = 0; ip < np; ++ip) {
    std::cout << ip << " --> ";
    for (int ipsup = index[ip]; ipsup < index[ip + 1]; ++ipsup) {
      std::cout << array[ipsup] << " ";
    }
    std::cout << std::endl;
  }
}

DFM2_INLINE void delfem2::JArray_Sort(
    const std::vector<unsigned int> &index,
    std::vector<unsigned int> &array) {
  if (index.empty()) return;
  const int size = (int) index.size() - 1;
  for (int ipoin = 0; ipoin < size; ipoin++) {
    const unsigned int is = index[ipoin];
    const unsigned int ie = index[ipoin + 1];
    if (is == ie) continue;  // no element in the row
    assert(is < ie);
    for (unsigned int i = is; i < ie - 1; i++) {
      for (unsigned int j = ie - 1; j != i; j--) {
        if (array[j] < array[j - 1]) {
          const unsigned int itmp = array[j];
          array[j] = array[j - 1];
          array[j - 1] = itmp;
        }
      }
    }
  }
}

DFM2_INLINE void delfem2::JArray_Sort(
    const unsigned int *index,
    const unsigned int size,
    unsigned int *array) {
  if (size == 0) return;
//  if( index.size() == 0 ) return;
//  const int size = (int)index.size()-1;
  for (unsigned int ipoin = 0; ipoin < size; ipoin++) {
    const unsigned int is = index[ipoin];
    const unsigned int ie = index[ipoin + 1];
    if (is == ie) continue;
    assert(is < ie);
    for (unsigned int i = is; i < ie - 1; i++) {
      for (int j = (int) ie - 1; j > (int) i; j--) {
        if (array[j] < array[j - 1]) {
          unsigned int itmp = array[j];
          array[j] = array[j - 1];
          array[j - 1] = itmp;
        }
      }
    }
  }
}

DFM2_INLINE void delfem2::JArray_AddDiagonal(
    std::vector<unsigned int> &psup_ind1,
    std::vector<unsigned int> &psup1,
    const unsigned int *psup_ind0,
    size_t npsup_ind0,
    const unsigned int *psup0,
    [[maybe_unused]] size_t npsup0) {
  assert(npsup_ind0>0);
  const size_t np = npsup_ind0 - 1;
  std::vector<unsigned int> tmp(np, UINT_MAX);
  psup_ind1.assign(np + 1, 0);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (std::size_t ipsup = psup_ind0[ip]; ipsup < psup_ind0[ip + 1]; ++ipsup) {
      const unsigned int jp = psup0[ipsup];
      assert(tmp[jp] != ip);
      tmp[jp] = ip;
      psup_ind1[ip + 1] += 1;
    }
    if (tmp[ip] != ip) {
      tmp[ip] = ip;
      psup_ind1[ip + 1] += 1;
    }
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    psup_ind1[ip + 1] += psup_ind1[ip];
  }
  const unsigned int npsup = psup_ind1[np];
  psup1.resize(npsup);
  tmp.assign(np, UINT_MAX);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (std::size_t ipsup = psup_ind0[ip]; ipsup < psup_ind0[ip + 1]; ++ipsup) {
      const unsigned int jp = psup0[ipsup];
      assert(tmp[jp] != ip);
      tmp[jp] = ip;
      unsigned int iclstr = psup_ind1[ip];
      psup1[iclstr] = jp;
      psup_ind1[ip] += 1;
    }
    if (tmp[ip] != ip) {
      unsigned int iclstr = psup_ind1[ip];
      psup1[iclstr] = ip;
      psup_ind1[ip] += 1;
    }
  }
  for (int ip = static_cast<int>(np) - 1; ip >= 0; --ip) {
    psup_ind1[ip + 1] = psup_ind1[ip];
  }
  psup_ind1[0] = 0;
}

// -------------------------------------------


/**
 * @details compute 2-ring neighborhood from 1-ring neighborhood
 */
DFM2_INLINE void delfem2::JArray_Extend(
    std::vector<unsigned int> &psup_ind1,
    std::vector<unsigned int> &psup1,
    const unsigned int *psup_ind0,
    size_t npsup_ind0,
    const unsigned int *psup0) {
  const size_t np = npsup_ind0 - 1;
  psup_ind1.assign(np + 1, 0);
  std::vector<unsigned> aflg(np, UINT_MAX);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind0[ip]; ipsup < psup_ind0[ip + 1]; ++ipsup) {
      unsigned int jp0 = psup0[ipsup];
      for (unsigned int jpsup = psup_ind0[jp0]; jpsup < psup_ind0[jp0 + 1]; ++jpsup) {
        unsigned int kp0 = psup0[jpsup];
        if (aflg[kp0] == ip || kp0 == ip) { continue; }
        ++psup_ind1[ip + 1];
        aflg[kp0] = ip;
      }
    }
  }
  // ---------
  for (unsigned int ip = 0; ip < np; ++ip) {
    psup_ind1[ip + 1] += psup_ind1[ip];
  }
  psup1.resize(psup_ind1[np]);
  // ---------
  aflg.assign(np, UINT_MAX);
  for (unsigned int ip = 0; ip < np; ++ip) {
    for (unsigned int ipsup = psup_ind0[ip]; ipsup < psup_ind0[ip + 1]; ++ipsup) {
      unsigned int jp0 = psup0[ipsup];
      for (unsigned int jpsup = psup_ind0[jp0]; jpsup < psup_ind0[jp0 + 1]; ++jpsup) {
        unsigned int kp0 = psup0[jpsup];
        if (aflg[kp0] == ip || kp0 == ip) { continue; }
        unsigned int kpsup = psup_ind1[ip];
        ++psup_ind1[ip];
        psup1[kpsup] = kp0;
        aflg[kp0] = ip;
      }
    }
  }
  for (int ip = static_cast<int>(np); ip >= 1; --ip) {
    psup_ind1[ip] = psup_ind1[ip - 1];
  }
  psup_ind1[0] = 0;
}

