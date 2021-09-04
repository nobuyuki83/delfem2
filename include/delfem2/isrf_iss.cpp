/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/isrf_iss.h"

#include <cassert>
#include <cmath>
#include <iostream>

namespace delfem2 {
namespace iss {

DFM2_INLINE void makeElemSurroundingPoint(
    std::vector<int> &elsup_ind,
    std::vector<int> &elsup,
    //
    const std::vector<unsigned int> &aElem,
    int nPoEl,
    int nPo) {
  const int nElem = (int) aElem.size() / nPoEl;
  elsup_ind.assign(nPo + 1, 0);
  for (int ielem = 0; ielem < nElem; ielem++) {
    for (int inoel = 0; inoel < nPoEl; inoel++) {
      int ino1 = aElem[ielem * nPoEl + inoel];
      if (ino1 == -1) { break; }
      elsup_ind[ino1 + 1] += 1;
    }
  }
  for (int ino = 0; ino < nPo; ++ino) {
    elsup_ind[ino + 1] += elsup_ind[ino];
  }
  int nelsup = elsup_ind[nPo];
  elsup.resize(nelsup);
  for (int ielem = 0; ielem < nElem; ielem++) {
    for (int inoel = 0; inoel < nPoEl; inoel++) {
      int ino1 = aElem[ielem * nPoEl + inoel];
      if (ino1 == -1) { break; }
      int ind1 = elsup_ind[ino1];
      elsup[ind1] = ielem;
      elsup_ind[ino1] += 1;
    }
  }
  for (int ino = nPo; ino >= 1; ino--) {
    elsup_ind[ino] = elsup_ind[ino - 1];
  }
  elsup_ind[0] = 0;
}

DFM2_INLINE void makeOneRingNeighborhood(
    std::vector<int> &psup_ind,
    std::vector<int> &psup,
    //
    const std::vector<unsigned int> &aElem,
    const std::vector<int> &elsup_ind,
    const std::vector<int> &elsup,
    int nnoel,
    int nnode) {
  std::vector<int> aflg(nnode, -1);
  psup_ind.assign(nnode + 1, 0);
  for (int inode = 0; inode < nnode; inode++) {
    aflg[inode] = inode;
    for (int ielsup = elsup_ind[inode]; ielsup < elsup_ind[inode + 1]; ielsup++) {
      int jelem = elsup[ielsup];
      for (int jnoel = 0; jnoel < nnoel; jnoel++) {
        int jnode = aElem[jelem * nnoel + jnoel];
        if (aflg[jnode] != inode) {
          aflg[jnode] = inode;
          psup_ind[inode + 1]++;
        }
      }
    }
  }
  for (int ino = 0; ino < nnode; ino++) {
    psup_ind[ino + 1] += psup_ind[ino];
  }
  const int npsup = psup_ind[nnode];
  psup.resize(npsup);
  for (int ino = 0; ino < nnode; ino++) { aflg[ino] = -1; }
  for (int inode = 0; inode < nnode; inode++) {
    aflg[inode] = inode;
    for (int ielsup = elsup_ind[inode]; ielsup < elsup_ind[inode + 1]; ielsup++) {
      int jelem = elsup[ielsup];
      for (int jnoel = 0; jnoel < nnoel; jnoel++) {
        int jnode = aElem[jelem * nnoel + jnoel];
        if (aflg[jnode] != inode) {
          aflg[jnode] = inode;
          const int ind = psup_ind[inode];
          psup[ind] = jnode;
          psup_ind[inode]++;
        }
      }
    }
  }
  for (int inode = nnode; inode > 0; inode--) {
    psup_ind[inode] = psup_ind[inode - 1];
  }
  psup_ind[0] = 0;
}


// ---------------------------------------------------------

DFM2_INLINE int GetCutNode(
    unsigned int iln0,
    unsigned int iln1,
    const std::vector<int> &aCutInd,
    const std::vector<int> &aCut) {
  for (int ind0 = aCutInd[iln0]; ind0 < aCutInd[iln0 + 1]; ind0++) {
    if (aCut[ind0 * 2 + 0] != (int) iln1) continue;
    return aCut[ind0 * 2 + 1];
  }
  return -1;
}

DFM2_INLINE void FindCutNodeTet(
    int on[10],
    unsigned int iln0,
    unsigned int iln1,
    unsigned int iln2,
    unsigned int iln3,
    int f0,
    int f1,
    int f2,
    int f3,
    const std::vector<int> &mapLat2Out,
    const std::vector<int> &aCutInd,
    const std::vector<int> &aCut) {
  on[0] = mapLat2Out[iln0];
  on[1] = mapLat2Out[iln1];
  on[2] = mapLat2Out[iln2];
  on[3] = mapLat2Out[iln3];
  if (f0 != f1 && f0 != 0 && f1 != 0) {
    on[4] = GetCutNode(iln0, iln1, aCutInd, aCut);
    assert(on[4] != -1);
  }
  if (f0 != f2 && f0 != 0 && f2 != 0) {
    on[5] = GetCutNode(iln0, iln2, aCutInd, aCut);
    assert(on[5] != -1);
  }
  if (f0 != f3 && f0 != 0 && f3 != 0) {
    on[6] = GetCutNode(iln0, iln3, aCutInd, aCut);
    assert(on[6] != -1);
  }
  if (f1 != f2 && f1 != 0 && f2 != 0) {
    on[7] = GetCutNode(iln1, iln2, aCutInd, aCut);
    assert(on[7] != -1);
  }
  if (f1 != f3 && f1 != 0 && f3 != 0) {
    on[8] = GetCutNode(iln1, iln3, aCutInd, aCut);
    assert(on[8] != -1);
  }
  if (f2 != f3 && f2 != 0 && f3 != 0) {
    on[9] = GetCutNode(iln2, iln3, aCutInd, aCut);
    assert(on[9] != -1);
  }
}

DFM2_INLINE void GetClampTet(
    int tet[3][4],
    unsigned int &ntet,
    unsigned int iflg,
    int on[10]) {
  ntet = 0;
  if (iflg == 40) { return; }
  if (iflg == 13 || iflg == 37 || iflg == 31 || iflg == 39) { return; }
  // --00
  if (iflg == 4 || iflg == 10 || iflg == 28 || iflg == 12 || iflg == 30 || iflg == 36) { return; }
  if (iflg == 80) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[1];
    tet[0][2] = on[2];
    tet[0][3] = on[3];
    return;
  }
  if (iflg == 41) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[4];
    tet[0][2] = on[5];
    tet[0][3] = on[6];
    return;
  }
  if (iflg == 43) {
    ntet = 1;
    tet[0][0] = on[1];
    tet[0][1] = on[7];
    tet[0][2] = on[4];
    tet[0][3] = on[8];
    return;
  }
  if (iflg == 49) {
    ntet = 1;
    tet[0][0] = on[2];
    tet[0][1] = on[5];
    tet[0][2] = on[7];
    tet[0][3] = on[9];
    return;
  }
  if (iflg == 67) {
    ntet = 1;
    tet[0][0] = on[3];
    tet[0][1] = on[6];
    tet[0][2] = on[9];
    tet[0][3] = on[8];
    return;
  }
  {
    unsigned int k1 = 10, k2=0, k3=0, k4=0, k5=0, k6=0;
    if (iflg == 79) {
      k1 = 1;
      k2 = 2;
      k3 = 3;
      k4 = 4;
      k5 = 5;
      k6 = 6;
    }
    if (iflg == 77) {
      k1 = 0;
      k2 = 3;
      k3 = 2;
      k4 = 4;
      k5 = 8;
      k6 = 7;
    }
    if (iflg == 71) {
      k1 = 0;
      k2 = 1;
      k3 = 3;
      k4 = 5;
      k5 = 7;
      k6 = 9;
    }
    if (iflg == 53) {
      k1 = 0;
      k2 = 2;
      k3 = 1;
      k4 = 6;
      k5 = 9;
      k6 = 8;
    }
    if (k1 != 10) {
      ntet = 3;
      if (on[k1] > on[k2] && on[k1] > on[k3]) {
        tet[0][0] = on[k1];
        tet[0][1] = on[k5];
        tet[0][2] = on[k4];
        tet[0][3] = on[k6];
        if (on[k2] > on[k3]) {
          tet[1][0] = on[k1];
          tet[1][1] = on[k2];
          tet[1][2] = on[k6];
          tet[1][3] = on[k3];
          tet[2][0] = on[k1];
          tet[2][1] = on[k2];
          tet[2][2] = on[k5];
          tet[2][3] = on[k6];
        } else {
          tet[1][0] = on[k1];
          tet[1][1] = on[k2];
          tet[1][2] = on[k5];
          tet[1][3] = on[k3];
          tet[2][0] = on[k1];
          tet[2][1] = on[k5];
          tet[2][2] = on[k6];
          tet[2][3] = on[k3];
        }
      }
      if (on[k2] > on[k1] && on[k2] > on[k3]) {
        tet[0][0] = on[k2];
        tet[0][1] = on[k6];
        tet[0][2] = on[k5];
        tet[0][3] = on[k4];
        if (on[k1] > on[k3]) {
          tet[1][0] = on[k1];
          tet[1][1] = on[k2];
          tet[1][2] = on[k6];
          tet[1][3] = on[k3];
          tet[2][0] = on[k1];
          tet[2][1] = on[k6];
          tet[2][2] = on[k2];
          tet[2][3] = on[k4];
        } else {
          tet[1][0] = on[k1];
          tet[1][1] = on[k2];
          tet[1][2] = on[k4];
          tet[1][3] = on[k3];
          tet[2][0] = on[k4];
          tet[2][1] = on[k2];
          tet[2][2] = on[k6];
          tet[2][3] = on[k3];
        }
      }
      if (on[k3] > on[k1] && on[k3] > on[k2]) {
        tet[0][0] = on[k3];
        tet[0][1] = on[k5];
        tet[0][2] = on[k4];
        tet[0][3] = on[k6];
        if (on[k1] > on[k2]) {
          tet[1][0] = on[k1];
          tet[1][1] = on[k3];
          tet[1][2] = on[k2];
          tet[1][3] = on[k5];
          tet[2][0] = on[k1];
          tet[2][1] = on[k3];
          tet[2][2] = on[k5];
          tet[2][3] = on[k4];
        } else {
          tet[1][0] = on[k1];
          tet[1][1] = on[k3];
          tet[1][2] = on[k2];
          tet[1][3] = on[k4];
          tet[2][0] = on[k2];
          tet[2][1] = on[k3];
          tet[2][2] = on[k5];
          tet[2][3] = on[k4];
        }
      }
      return;
    }
  }
  if (iflg == 44 || iflg == 50 || iflg == 68 || iflg == 52 || iflg == 70 || iflg == 76) {
    ntet = 3;
    int k1 = -1, k2 = -1, k3 = -1, k4 = -1, k5 = -1, k6 = -1;
    if (iflg == 44) {
      k1 = 0;
      k2 = 1;
      k3 = 5;
      k4 = 6;
      k5 = 7;
      k6 = 8;
    }
    else if (iflg == 50) {
      k1 = 0;
      k2 = 2;
      k3 = 6;
      k4 = 4;
      k5 = 9;
      k6 = 7;
    }
    else if (iflg == 68) {
      k1 = 0;
      k2 = 3;
      k3 = 4;
      k4 = 5;
      k5 = 8;
      k6 = 9;
    }
    else if (iflg == 52) {
      k1 = 1;
      k2 = 2;
      k3 = 4;
      k4 = 8;
      k5 = 5;
      k6 = 9;
    }
    else if (iflg == 70) {
      k1 = 1;
      k2 = 3;
      k3 = 7;
      k4 = 4;
      k5 = 9;
      k6 = 6;
    }
    else if (iflg == 76) {
      k1 = 2;
      k2 = 3;
      k3 = 5;
      k4 = 7;
      k5 = 6;
      k6 = 8;
    }
    if (on[k1] > on[k2]) {
      tet[0][0] = on[k1];
      tet[0][1] = on[k2];
      tet[0][2] = on[k5];
      tet[0][3] = on[k6];
      tet[1][0] = on[k1];
      tet[1][1] = on[k6];
      tet[1][2] = on[k3];
      tet[1][3] = on[k4];
      tet[2][0] = on[k1];
      tet[2][1] = on[k3];
      tet[2][2] = on[k6];
      tet[2][3] = on[k5];
    } else {
      tet[0][0] = on[k1];
      tet[0][1] = on[k2];
      tet[0][2] = on[k3];
      tet[0][3] = on[k4];
      tet[1][0] = on[k2];
      tet[1][1] = on[k6];
      tet[1][2] = on[k3];
      tet[1][3] = on[k4];
      tet[2][0] = on[k2];
      tet[2][1] = on[k3];
      tet[2][2] = on[k6];
      tet[2][3] = on[k5];
    }
    return;
  }
  // +000
  if (iflg == 2 || iflg == 6 || iflg == 18 || iflg == 54) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[1];
    tet[0][2] = on[2];
    tet[0][3] = on[3];
    return;
  }
  // +++0
  if (iflg == 78 || iflg == 74 || iflg == 62 || iflg == 26) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[1];
    tet[0][2] = on[2];
    tet[0][3] = on[3];
    return;
  }
  // ++00
  if (iflg == 8 || iflg == 20 || iflg == 56 || iflg == 24 || iflg == 60 || iflg == 72) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[1];
    tet[0][2] = on[2];
    tet[0][3] = on[3];
    return;
  }

  // -+00
  if (iflg == 5) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[4];
    tet[0][2] = on[2];
    tet[0][3] = on[3];
    return;
  }
  if (iflg == 11) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[1];
    tet[0][2] = on[5];
    tet[0][3] = on[3];
    return;
  }
  if (iflg == 29) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[1];
    tet[0][2] = on[2];
    tet[0][3] = on[6];
    return;
  }
  if (iflg == 7) {
    ntet = 1;
    tet[0][0] = on[4];
    tet[0][1] = on[1];
    tet[0][2] = on[2];
    tet[0][3] = on[3];
    return;
  }
  if (iflg == 15) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[1];
    tet[0][2] = on[7];
    tet[0][3] = on[3];
    return;
  }
  if (iflg == 33) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[1];
    tet[0][2] = on[2];
    tet[0][3] = on[8];
    return;
  }
  if (iflg == 19) {
    ntet = 1;
    tet[0][0] = on[5];
    tet[0][1] = on[1];
    tet[0][2] = on[2];
    tet[0][3] = on[3];
    return;
  }
  if (iflg == 21) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[7];
    tet[0][2] = on[2];
    tet[0][3] = on[3];
    return;
  }
  if (iflg == 45) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[1];
    tet[0][2] = on[2];
    tet[0][3] = on[9];
    return;
  }
  if (iflg == 55) {
    ntet = 1;
    tet[0][0] = on[6];
    tet[0][1] = on[1];
    tet[0][2] = on[2];
    tet[0][3] = on[3];
    return;
  }
  if (iflg == 57) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[8];
    tet[0][2] = on[2];
    tet[0][3] = on[3];
    return;
  }
  if (iflg == 63) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[1];
    tet[0][2] = on[9];
    tet[0][3] = on[3];
    return;
  }

  // --+0
  if (iflg == 14) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[4];
    tet[0][2] = on[5];
    tet[0][3] = on[3];
    return;
  } // +--0
  if (iflg == 32) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[4];
    tet[0][2] = on[2];
    tet[0][3] = on[6];
    return;
  } // +-0-
  if (iflg == 38) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[1];
    tet[0][2] = on[5];
    tet[0][3] = on[6];
    return;
  } // +0--
  if (iflg == 16) {
    ntet = 1;
    tet[0][0] = on[1];
    tet[0][1] = on[7];
    tet[0][2] = on[4];
    tet[0][3] = on[3];
    return;
  } // -+-0
  if (iflg == 34) {
    ntet = 1;
    tet[0][0] = on[1];
    tet[0][1] = on[2];
    tet[0][2] = on[4];
    tet[0][3] = on[8];
    return;
  } // -+0-
  if (iflg == 42) {
    ntet = 1;
    tet[0][0] = on[1];
    tet[0][1] = on[7];
    tet[0][2] = on[0];
    tet[0][3] = on[8];
    return;
  } // 0+--
  if (iflg == 22) {
    ntet = 1;
    tet[0][0] = on[2];
    tet[0][1] = on[5];
    tet[0][2] = on[7];
    tet[0][3] = on[3];
    return;
  } // --+0
  if (iflg == 46) {
    ntet = 1;
    tet[0][0] = on[1];
    tet[0][1] = on[2];
    tet[0][2] = on[5];
    tet[0][3] = on[9];
    return;
  } // -0+-
  if (iflg == 48) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[7];
    tet[0][2] = on[2];
    tet[0][3] = on[9];
    return;
  } // 0-+-
  if (iflg == 58) {
    ntet = 1;
    tet[0][0] = on[6];
    tet[0][1] = on[8];
    tet[0][2] = on[2];
    tet[0][3] = on[3];
    return;
  } // --0+
  if (iflg == 64) {
    ntet = 1;
    tet[0][0] = on[1];
    tet[0][1] = on[9];
    tet[0][2] = on[6];
    tet[0][3] = on[3];
    return;
  } // -0-+
  if (iflg == 66) {
    ntet = 1;
    tet[0][0] = on[0];
    tet[0][1] = on[8];
    tet[0][2] = on[9];
    tet[0][3] = on[3];
    return;
  } // 0--+
  ////
  if (iflg == 73 || iflg == 75 || iflg == 25 || iflg == 23 || iflg == 59 || iflg == 61 || iflg == 17
      || iflg == 35 || iflg == 47 || iflg == 51 || iflg == 65 || iflg == 69) {
    ntet = 2;
    int k1 = -1, k2 = -1, k3 = -1, k4 = -1, k5 = -1;
    if (iflg == 75) {
      k1 = 0;
      k2 = 3;
      k3 = 2;
      k4 = 8;
      k5 = 7;
    } // 0-++
    else if (iflg == 69) {
      k1 = 0;
      k2 = 1;
      k3 = 3;
      k4 = 7;
      k5 = 9;
    } // 0+-+
    else if (iflg == 51) {
      k1 = 0;
      k2 = 2;
      k3 = 1;
      k4 = 9;
      k5 = 8;
    } // 0++-
    else if (iflg == 73) {
      k1 = 1;
      k2 = 2;
      k3 = 3;
      k4 = 5;
      k5 = 6;
    } // -0++
    else if (iflg == 65) {
      k1 = 1;
      k2 = 3;
      k3 = 0;
      k4 = 9;
      k5 = 5;
    } // +0-+
    else if (iflg == 47) {
      k1 = 1;
      k2 = 0;
      k3 = 2;
      k4 = 6;
      k5 = 9;
    } // +0+-
    else if (iflg == 61) {
      k1 = 2;
      k2 = 3;
      k3 = 1;
      k4 = 6;
      k5 = 4;
    } // -+0+
    else if (iflg == 59) {
      k1 = 2;
      k2 = 0;
      k3 = 3;
      k4 = 4;
      k5 = 8;
    } // +-0+
    else if (iflg == 35) {
      k1 = 2;
      k2 = 1;
      k3 = 0;
      k4 = 8;
      k5 = 6;
    } // ++0-
    else if (iflg == 25) {
      k1 = 3;
      k2 = 1;
      k3 = 2;
      k4 = 4;
      k5 = 5;
    } // -++0
    else if (iflg == 23) {
      k1 = 3;
      k2 = 2;
      k3 = 0;
      k4 = 7;
      k5 = 4;
    } // +-+0
    else if (iflg == 17) {
      k1 = 3;
      k2 = 0;
      k3 = 1;
      k4 = 5;
      k5 = 7;
    } // ++-0
    if (on[k2] > on[k3]) {
      tet[0][0] = on[k1];
      tet[0][1] = on[k3];
      tet[0][2] = on[k2];
      tet[0][3] = on[k5];
      tet[1][0] = on[k1];
      tet[1][1] = on[k5];
      tet[1][2] = on[k2];
      tet[1][3] = on[k4];
    } else {
      tet[0][0] = on[k1];
      tet[0][1] = on[k3];
      tet[0][2] = on[k2];
      tet[0][3] = on[k4];
      tet[1][0] = on[k1];
      tet[1][1] = on[k3];
      tet[1][2] = on[k4];
      tet[1][3] = on[k5];
    }
    return;
  }
}

// 0-7 corner
// 8-19 edge ---- 8:01, 9:23, 10:45, 11:67, 12:02, 13:13, 14:46, 15:57, 16:04, 17:15, 18:26, 19:37
// 20-25 face ----- 20:0462, 21:1375, 22:0154, 23:2673, 24:0231, 25:4576
// 26 center ------ 26:01234567
const double aCellPointDirection[27][3] = {
    // corner
    {-1.0, -1.0, -1.0},
    {+1.0, -1.0, -1.0},
    {-1.0, +1.0, -1.0},
    {+1.0, +1.0, -1.0},
    {-1.0, -1.0, +1.0},
    {+1.0, -1.0, +1.0},
    {-1.0, +1.0, +1.0},
    {+1.0, +1.0, +1.0},
    // edge
    {0.0, -1.0, -1.0}, // 8
    {0.0, +1.0, -1.0}, // 9
    {0.0, -1.0, +1.0}, // 10
    {0.0, +1.0, +1.0}, // 11
    {-1.0, 0.0, -1.0}, // 12
    {+1.0, 0.0, -1.0}, // 13
    {-1.0, 0.0, +1.0}, // 14
    {+1.0, 0.0, +1.0}, // 15
    {-1.0, -1.0, 0.0}, // 16
    {+1.0, -1.0, 0.0}, // 17
    {-1.0, +1.0, 0.0}, // 18
    {+1.0, +1.0, 0.0}, // 19
    // face
    {-1.0, 0.0, 0.0},
    {+1.0, 0.0, 0.0},
    {0.0, -1.0, 0.0},
    {0.0, +1.0, 0.0},
    {0.0, 0.0, -1.0},
    {0.0, 0.0, +1.0},
    // cntr
    {0.0, 0.0, 0.0}
};

class CCell {
 public:
  CCell()
      : size(1), ilevel(0), iparent(-1), iparent_pos(-1) {
    for (int i = 0; i < 27; ++i) { aIP[i] = -1; }
    for (int i = 0; i < 8; ++i) { aIC_Cld[i] = -1; }
    for (int i = 0; i < 6; ++i) { aIC_Adj[i] = -1; }
  }
  CCell(double size, int ilevel, int iparent, int iparent_pos)
      : size(size), ilevel(ilevel), iparent(iparent), iparent_pos(iparent_pos) {
    for (int i = 0; i < 27; ++i) { aIP[i] = -1; }
    for (int i = 0; i < 8; ++i) { aIC_Cld[i] = -1; }
    for (int i = 0; i < 6; ++i) { aIC_Adj[i] = -1; }
  }
  bool isHavingChild_Face(int iface) const {
    assert(iface >= 0 && iface < 6);
    const int faceHex[6][4] = {
        {0, 4, 6, 2},
        {1, 3, 7, 5},
        {0, 1, 5, 4},
        {2, 6, 7, 3},
        {0, 2, 3, 1},
        {4, 5, 7, 6}};
    if (aIC_Cld[faceHex[iface][0]] >= 0) return true;
    if (aIC_Cld[faceHex[iface][1]] >= 0) return true;
    if (aIC_Cld[faceHex[iface][2]] >= 0) return true;
    if (aIC_Cld[faceHex[iface][3]] >= 0) return true;
    return false;
  }
  bool isHavingChild_Edge(int iedge) const {
    assert(iedge >= 0 && iedge < 12);
    const int edgeHex[12][2] = {
        {0, 1}, {2, 3}, {4, 5}, {6, 7},
        {0, 2}, {1, 3}, {4, 6}, {5, 7},
        {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    if (aIC_Cld[edgeHex[iedge][0]] >= 0) return true;
    if (aIC_Cld[edgeHex[iedge][1]] >= 0) return true;
    return false;
  }
  bool isHavingChild() const {
    if (aIC_Cld[0] >= 0) return true;
    if (aIC_Cld[1] >= 0) return true;
    if (aIC_Cld[2] >= 0) return true;
    if (aIC_Cld[3] >= 0) return true;
    if (aIC_Cld[4] >= 0) return true;
    if (aIC_Cld[5] >= 0) return true;
    if (aIC_Cld[6] >= 0) return true;
    if (aIC_Cld[7] >= 0) return true;
    return false;
  }
  void setChildAdjRelation(std::vector<CCell> &aCell) {
    int icc0 = aIC_Cld[0];
    int icc1 = aIC_Cld[1];
    int icc2 = aIC_Cld[2];
    int icc3 = aIC_Cld[3];
    int icc4 = aIC_Cld[4];
    int icc5 = aIC_Cld[5];
    int icc6 = aIC_Cld[6];
    int icc7 = aIC_Cld[7];
    // x
    if (icc0 != -1 && icc1 != -1) {
      aCell[icc0].aIC_Adj[1] = icc1;
      aCell[icc1].aIC_Adj[0] = icc0;
    }
    if (icc2 != -1 && icc3 != -1) {
      aCell[icc2].aIC_Adj[1] = icc3;
      aCell[icc3].aIC_Adj[0] = icc2;
    }
    if (icc4 != -1 && icc5 != -1) {
      aCell[icc4].aIC_Adj[1] = icc5;
      aCell[icc5].aIC_Adj[0] = icc4;
    }
    if (icc6 != -1 && icc7 != -1) {
      aCell[icc6].aIC_Adj[1] = icc7;
      aCell[icc7].aIC_Adj[0] = icc6;
    }
    // y
    if (icc0 != -1 && icc2 != -1) {
      aCell[icc0].aIC_Adj[3] = icc2;
      aCell[icc2].aIC_Adj[2] = icc0;
    }
    if (icc1 != -1 && icc3 != -1) {
      aCell[icc1].aIC_Adj[3] = icc3;
      aCell[icc3].aIC_Adj[2] = icc1;
    }
    if (icc4 != -1 && icc6 != -1) {
      aCell[icc4].aIC_Adj[3] = icc6;
      aCell[icc6].aIC_Adj[2] = icc4;
    }
    if (icc5 != -1 && icc7 != -1) {
      aCell[icc5].aIC_Adj[3] = icc7;
      aCell[icc7].aIC_Adj[2] = icc5;
    }
    // z
    if (icc0 != -1 && icc4 != -1) {
      aCell[icc0].aIC_Adj[5] = icc4;
      aCell[icc4].aIC_Adj[4] = icc0;
    }
    if (icc1 != -1 && icc5 != -1) {
      aCell[icc1].aIC_Adj[5] = icc5;
      aCell[icc5].aIC_Adj[4] = icc1;
    }
    if (icc2 != -1 && icc6 != -1) {
      aCell[icc2].aIC_Adj[5] = icc6;
      aCell[icc6].aIC_Adj[4] = icc2;
    }
    if (icc3 != -1 && icc7 != -1) {
      aCell[icc3].aIC_Adj[5] = icc7;
      aCell[icc7].aIC_Adj[4] = icc3;
    }
  }
 public:
  double size;
  int ilevel;
  ////
  //  bool is_leaf;
  int iparent;
  int iparent_pos;
  ////
  int aIC_Cld[8]; // child cell
  int aIC_Adj[6]; // child cell
  int aIP[27]; // id of corner/center points
};


/*
bool isEdgeCenterPoint
(const std::vector<CCell>& aCell,
 int ic, int iedge)
{
  assert( ic >=0 && ic < aCell.size() );
  assert( iedge >= 0 && iedge < 12 );
  const int edgeHex[12][2] = {
    {0,1}, {2,3}, {4,5}, {6,7},
    {0,2}, {1,3}, {4,6}, {5,7},
    {0,4}, {1,5}, {2,6}, {3,7} };
  int icc0 = aCell[ic].aIC_Cld[ edgeHex[iedge][0] ];
  int icc1 = aCell[ic].aIC_Cld[ edgeHex[iedge][1] ];
  if( icc0 >= 0 || icc1 >= 0 ) return true;
  return false;
}
 */

DFM2_INLINE void makeLatticeCoasestLevel(
    std::vector<CPointLattice> &aPoint,
    std::vector<CCell> &aCell,
    //
    const CInput_IsosurfaceStuffing &input,
    double elen,
    int ndiv,
    const double org[3]) {
  aPoint.clear();
  aPoint.resize((ndiv + 1) * (ndiv + 1) * (ndiv + 1) + ndiv * ndiv * ndiv);
  for (int iz = 0; iz < ndiv + 1; ++iz) {
    for (int iy = 0; iy < ndiv + 1; ++iy) {
      for (int ix = 0; ix < ndiv + 1; ++ix) {
        const int icrnr0 = iz * (ndiv + 1) * (ndiv + 1) + iy * (ndiv + 1) + ix;
        double cx = ix * elen + org[0];
        double cy = iy * elen + org[1];
        double cz = iz * elen + org[2];
        double sdf = input.SignedDistance(cx, cy, cz);
        aPoint[icrnr0] = CPointLattice(cx, cy, cz, sdf);
      }
    }
  }
  for (int iz = 0; iz < ndiv; ++iz) {
    for (int iy = 0; iy < ndiv; ++iy) {
      for (int ix = 0; ix < ndiv; ++ix) {
        int ip0 = iz * ndiv * ndiv + iy * ndiv + ix + (ndiv + 1) * (ndiv + 1) * (ndiv + 1);
        double cx = (ix + 0.5) * elen + org[0];
        double cy = (iy + 0.5) * elen + org[1];
        double cz = (iz + 0.5) * elen + org[2];
        double sdf = input.SignedDistance(cx, cy, cz);
        aPoint[ip0] = CPointLattice(cx, cy, cz, sdf);
      }
    }
  }
  aCell.clear();
  aCell.reserve(ndiv * ndiv * ndiv * 2);
  for (int iz = 0; iz < ndiv; ++iz) {
    for (int iy = 0; iy < ndiv; ++iy) {
      for (int ix = 0; ix < ndiv; ++ix) {
        CCell c(elen, 0, -1, -1);
        c.aIP[0] = (iz + 0) * (ndiv + 1) * (ndiv + 1) + (iy + 0) * (ndiv + 1) + (ix + 0);
        c.aIP[1] = (iz + 0) * (ndiv + 1) * (ndiv + 1) + (iy + 0) * (ndiv + 1) + (ix + 1);
        c.aIP[2] = (iz + 0) * (ndiv + 1) * (ndiv + 1) + (iy + 1) * (ndiv + 1) + (ix + 0);
        c.aIP[3] = (iz + 0) * (ndiv + 1) * (ndiv + 1) + (iy + 1) * (ndiv + 1) + (ix + 1);
        c.aIP[4] = (iz + 1) * (ndiv + 1) * (ndiv + 1) + (iy + 0) * (ndiv + 1) + (ix + 0);
        c.aIP[5] = (iz + 1) * (ndiv + 1) * (ndiv + 1) + (iy + 0) * (ndiv + 1) + (ix + 1);
        c.aIP[6] = (iz + 1) * (ndiv + 1) * (ndiv + 1) + (iy + 1) * (ndiv + 1) + (ix + 0);
        c.aIP[7] = (iz + 1) * (ndiv + 1) * (ndiv + 1) + (iy + 1) * (ndiv + 1) + (ix + 1);
        c.aIP[26] = iz * ndiv * ndiv + iy * ndiv + ix + (ndiv + 1) * (ndiv + 1) * (ndiv + 1); // center
        if (ix != 0) { c.aIC_Adj[0] = (iz + 0) * ndiv * ndiv + (iy + 0) * ndiv + (ix - 1); } else { c.aIC_Adj[0] = -2; }
        if (ix != ndiv - 1) { c.aIC_Adj[1] = (iz + 0) * ndiv * ndiv + (iy + 0) * ndiv + (ix + 1); }
        else {
          c.aIC_Adj[1] = -2;
        }
        if (iy != 0) { c.aIC_Adj[2] = (iz + 0) * ndiv * ndiv + (iy - 1) * ndiv + (ix + 0); } else { c.aIC_Adj[2] = -2; }
        if (iy != ndiv - 1) { c.aIC_Adj[3] = (iz + 0) * ndiv * ndiv + (iy + 1) * ndiv + (ix + 0); }
        else {
          c.aIC_Adj[3] = -2;
        }
        if (iz != 0) { c.aIC_Adj[4] = (iz - 1) * ndiv * ndiv + (iy + 0) * ndiv + (ix + 0); } else { c.aIC_Adj[4] = -2; }
        if (iz != ndiv - 1) { c.aIC_Adj[5] = (iz + 1) * ndiv * ndiv + (iy + 0) * ndiv + (ix + 0); }
        else {
          c.aIC_Adj[5] = -2;
        }
        aCell.push_back(c);
      }
    }
  }
}

DFM2_INLINE void makeChild(
    std::vector<CCell> &aCell,
    std::vector<CPointLattice> &aPoint,
    const CInput_IsosurfaceStuffing &input,
    unsigned int icell,
    int ichild) {
  assert(icell < aCell.size());
  assert(ichild >= 0 && ichild < 8);
  if (aCell[icell].aIC_Cld[ichild] != -1) return; // already there
  int icntr0 = aCell[icell].aIP[26];
  const double dir[8][3] = {
      {-1.0, -1.0, -1.0},
      {+1.0, -1.0, -1.0},
      {-1.0, +1.0, -1.0},
      {+1.0, +1.0, -1.0},
      {-1.0, -1.0, +1.0},
      {+1.0, -1.0, +1.0},
      {-1.0, +1.0, +1.0},
      {+1.0, +1.0, +1.0}};
  double size0 = aCell[icell].size;
  double cx = aPoint[icntr0].pos[0];
  double cy = aPoint[icntr0].pos[1];
  double cz = aPoint[icntr0].pos[2];
  int ilevel0 = aCell[icell].ilevel;
  const double ccx = cx + dir[ichild][0] * size0 * 0.25;
  const double ccy = cy + dir[ichild][1] * size0 * 0.25;
  const double ccz = cz + dir[ichild][2] * size0 * 0.25;
  aCell[icell].aIC_Cld[ichild] = (int) aCell.size();
  CCell cc(size0 * 0.5, ilevel0 + 1, icell, ichild);
  {
    cc.aIP[26] = (int) aPoint.size();
    CPointLattice p(ccx, ccy, ccz, input.SignedDistance(ccx, ccy, ccz));
    aPoint.push_back(p);
  }
  aCell.push_back(cc);
  aCell[icell].setChildAdjRelation(aCell);
}

DFM2_INLINE void makeChild_Face(
    std::vector<CCell> &aCell,
    std::vector<CPointLattice> &aPoint,
    const CInput_IsosurfaceStuffing &input,
    int icell,
    int iface) {
  assert(icell >= 0 && icell < (int) aCell.size());
  const int faceHex[6][4] = {
      {0, 4, 6, 2},
      {1, 3, 7, 5},
      {0, 1, 5, 4},
      {2, 6, 7, 3},
      {0, 2, 3, 1},
      {4, 5, 7, 6}};
  for (int ifc = 0; ifc < 4; ++ifc) { // face child
    int ichild = faceHex[iface][ifc];
    if (aCell[icell].aIC_Cld[ichild] != -1) continue;
    makeChild(aCell, aPoint, input, icell, ichild);
  }
}

DFM2_INLINE void Continuation(
    std::vector<CPointLattice> &aPoint, std::vector<CCell> &aCell,
    const CInput_IsosurfaceStuffing &input) {
  const int faceHex[6][4] = {
      {0, 4, 6, 2},
      {1, 3, 7, 5},
      {0, 1, 5, 4},
      {2, 6, 7, 3},
      {0, 2, 3, 1},
      {4, 5, 7, 6}};
  const int oppFace[6] = {1, 0, 3, 2, 5, 4};
  const int adjFacePair[6][4][2] = {
      {{0, 1}, {2, 3}, {4, 5}, {6, 7}},
      {{1, 0}, {3, 2}, {5, 4}, {7, 6}},
      {{0, 2}, {1, 3}, {4, 6}, {5, 7}},
      {{2, 0}, {3, 1}, {6, 4}, {7, 5}},
      {{0, 4}, {1, 5}, {2, 6}, {3, 7}},
      {{4, 0}, {5, 1}, {6, 2}, {7, 3}}};
  const int adjChildInside[8][6] = {
      {-1, +1, -1, +2, -1, +4},
      {+0, -1, -1, +3, -1, +5},
      {-1, +3, +0, -1, -1, +6},
      {+2, -1, +1, -1, -1, +7},
      {-1, +5, -1, +6, +0, -1},
      {+4, -1, -1, +7, +1, -1},
      {-1, +7, +4, -1, +2, -1},
      {+6, -1, +5, -1, +3, -1}};
  const int edge2Face[12][4] = { // [f1,f2, e1,2], edge is blong to face f1,f2, for adjacent face, the edge is e1 and e2
      {2, 4, 1, 2}, // 0 8
      {3, 4, 0, 3}, // 1 9
      {2, 5, 3, 0}, // 2 10
      {3, 5, 2, 1}, // 3 11
      {0, 4, 5, 6}, // 4 12
      {1, 4, 4, 7}, // 5 13
      {0, 5, 7, 4}, // 6 14
      {1, 5, 6, 5}, // 7 15
      {0, 2, 9, 10}, // 8 16
      {1, 2, 8, 11}, // 9 17
      {0, 3, 11, 8}, //10 18
      {1, 3, 10, 9}};//11 19
  const int itr_continuation_max = 6; // this should be enough unless one cell is adjacent to 6 times finer cells
  int itr_continuation = 0;
  for (itr_continuation = 0; itr_continuation < itr_continuation_max; ++itr_continuation) {
    const int icnt0 = (int) aPoint.size();
    for (int icell = 0; icell < (int) aCell.size(); ++icell) {
      const int ipc = aCell[icell].iparent;
      if (ipc < 0) { continue; } // this cell is root (somebody will take care of this)
      for (int iedge = 0; iedge < 12; ++iedge) {
        if (!aCell[icell].isHavingChild_Edge(iedge)) { continue; }
        ////
        const int ich = aCell[icell].iparent_pos;
        {
          const int jface = edge2Face[iedge][0];
          if (aCell[icell].aIC_Adj[jface] < 0) { // find&make this cell
            const int jch = adjChildInside[ich][jface]; // inside neighbor
            const int jpca = aCell[ipc].aIC_Adj[jface]; // outside neighbor
            if (jch >= 0) { makeChild(aCell, aPoint, input, ipc, jch); }
            else if (jpca >= 0) { makeChild_Face(aCell, aPoint, input, jpca, oppFace[jface]); }
          }
        }
        {
          const int kface = edge2Face[iedge][1];
          if (aCell[icell].aIC_Adj[kface] < 0) { // find&make this cell
            const int kch = adjChildInside[ich][kface]; // inside neighbor
            const int kpca = aCell[ipc].aIC_Adj[kface]; // outside neighbor
            if (kch >= 0) { makeChild(aCell, aPoint, input, ipc, kch); }
            else if (kpca >= 0) { makeChild_Face(aCell, aPoint, input, kpca, oppFace[kface]); }
          }
        }
      }
    }
    // -------------
    for (int icell = 0; icell < (int) aCell.size(); ++icell) {
      const int ipc = aCell[icell].iparent;
      if (ipc < 0) { continue; } // this cell is root (somebody will take care of this)
      for (int iedge = 0; iedge < 12; ++iedge) {
        if (!aCell[icell].isHavingChild_Edge(iedge)) { continue; }
        const int ich = aCell[icell].iparent_pos;
        const int jface = edge2Face[iedge][0];
        const int kface = edge2Face[iedge][1];
        const int jca = aCell[icell].aIC_Adj[jface];
        const int kca = aCell[icell].aIC_Adj[kface];
        const int lkca = (kca >= 0) ? aCell[kca].aIC_Adj[jface] : -1;
        const int ljca = (jca >= 0) ? aCell[jca].aIC_Adj[kface] : -1;
        if (lkca != -1 || ljca != -1) { continue; }
        const int jch = adjChildInside[ich][jface];
        int lch = (jch >= 0) ? adjChildInside[jch][kface] : -1;
        if (lch != -1) {
#ifndef NDEBUG
          const int kch = adjChildInside[ich][kface];
          assert(kch != -1);
          assert(lch == adjChildInside[kch][jface]);
#endif
        }
        if (lch >= 0) { // diagonal child
          makeChild(aCell, aPoint, input, ipc, lch);
        }
        {
          const int jpca = aCell[ipc].aIC_Adj[jface];
          const int ljkpca = (jpca >= 0) ? aCell[jpca].aIC_Adj[kface] : -1;
          if (ljkpca >= 0) {
            makeChild_Face(aCell, aPoint, input, ljkpca, oppFace[kface]);
            makeChild_Face(aCell, aPoint, input, ljkpca, oppFace[jface]);
          }
          ////
          const int kpca = aCell[ipc].aIC_Adj[kface];
          const int lkjpca = (kpca >= 0) ? aCell[kpca].aIC_Adj[jface] : -1;
          if (lkjpca >= 0) {
            makeChild_Face(aCell, aPoint, input, lkjpca, oppFace[kface]);
            makeChild_Face(aCell, aPoint, input, lkjpca, oppFace[jface]);
          }
        }
      }
    }
    // -------
    for (int icell = 0; icell < (int) aCell.size(); ++icell) {
      for (int iface = 0; iface < 6; ++iface) {
        if (!aCell[icell].isHavingChild_Face(iface)) { continue; }
        const int jca0 = aCell[icell].aIC_Adj[iface];
        if (jca0 < 0) continue;
        const int jface = oppFace[iface];
        if (aCell[jca0].isHavingChild_Face(jface)) {
          for (int jfc = 0; jfc < 4; ++jfc) { // face child
            int jc_ch0 = aCell[jca0].aIC_Cld[faceHex[jface][jfc]];
            if (jc_ch0 < 0) { continue; }
//            if( aCell[jc_ch0].isHavingChild_Face(jface) ){
            if (aCell[jc_ch0].isHavingChild()) { // why not above?
              makeChild_Face(aCell, aPoint, input, icell, iface);
              break;
            }
          }
        } else {
          for (int ifc = 0; ifc < 4; ++ifc) { // face child
            int ic_ch0 = aCell[icell].aIC_Cld[faceHex[iface][ifc]];
            if (ic_ch0 < 0) { continue; }
//            if( aCell[ic_ch0].isHavingChild_Face(iface) ){
            if (aCell[ic_ch0].isHavingChild()) { // why not above?
              makeChild_Face(aCell, aPoint, input, jca0, jface);
              break;
            }
          }
        }
      }
    }
    // -------------
    // making adjacent to the neighboring cells
    for (int ic = 0; ic < (int) aCell.size(); ++ic) {
      for (int iface = 0; iface < 6; ++iface) {
        const int jface = oppFace[iface];
        int ica = aCell[ic].aIC_Adj[iface];
        if (ica < 0) { continue; }
        assert(aCell[ic].ilevel == aCell[ica].ilevel);
        for (int ifc = 0; ifc < 4; ++ifc) {
          const int ich0 = adjFacePair[iface][ifc][0];
          const int icha0 = adjFacePair[iface][ifc][1];
          const int icc0 = aCell[ic].aIC_Cld[ich0];
          const int icca0 = aCell[ica].aIC_Cld[icha0];
          if (icc0 != -1 && icca0 != -1) {
            aCell[icc0].aIC_Adj[iface] = icca0;
            aCell[icca0].aIC_Adj[jface] = icc0;
          }
        }
      }
    }
    std::cout << "point count: " << icnt0 << " " << aPoint.size() << std::endl;
    if (icnt0 == (int) aPoint.size()) break;
  }

  if (itr_continuation == itr_continuation_max) {
    std::cout << "too much continuation... we gave up" << std::endl;
  }
}
/*
void CheckContinuation(const std::vector<CCell>& aCell)
{
  const int faceHex[6][4] = {
    {0,4,6,2},
    {1,3,7,5},
    {0,1,5,4},
    {2,6,7,3},
    {0,2,3,1},
    {4,5,7,6} };
  const int oppDir[6] = {1,0,3,2,5,4};
  ////
  for(int icell=0;icell<(int)aCell.size();++icell){
    for(int iface=0;iface<6;++iface){
      if( !aCell[icell].isHavingChild_Face(iface) ) continue;
      int jca0 = aCell[icell].aIC_Adj[iface];
      if( jca0 < 0 ) continue;
      const int jface = oppDir[iface];
      if( aCell[jca0].isHavingChild_Face(jface) ){ continue; }
      for(int ifc=0;ifc<4;++ifc){
        int ic_ch0 = aCell[icell].aIC_Cld[ faceHex[iface][ifc] ];
        if( ic_ch0 < 0 ) continue;
        if( aCell[ic_ch0].isHavingChild_Face(iface) ){
          std::cout << "something is wrong " << icell << " " << iface << std::endl;
        }
      }
    }
  }
  //////
  const int edge2Face[12][4] = { // [f1,f2, e1,2], edge is blong to face f1,f2, for adjacent face, the edge is e1 and e2
    {2,4, 1,2}, // 0 8
    {3,4, 0,3}, // 1 9
    {2,5, 3,0}, // 2 10
    {3,5, 2,1}, // 3 11
    {0,4, 5,6}, // 4 12
    {1,4, 4,7}, // 5 13
    {0,5, 7,4}, // 6 14
    {1,5, 6,5}, // 7 15
    {0,2, 9,10}, // 8 16
    {1,2, 8,11}, // 9 17
    {0,3, 11,8}, //10 18
    {1,3, 10,9} };//11 19
  const int edgeHex[12][2] = {
    {0,1}, {2,3}, {4,5}, {6,7},
    {0,2}, {1,3}, {4,6}, {5,7},
    {0,4}, {1,5}, {2,6}, {3,7} };
  for(int icell=0;icell<(int)aCell.size();++icell){
    for(int iedge=0;iedge<12;++iedge){
      if( !aCell[icell].isHavingChild_Edge(iedge) ){ continue; }
      { // jface
        const int jca0 = aCell[icell].aIC_Adj[ edge2Face[iedge][0] ];
        if( jca0 >= 0 ){
          const int jedge = edge2Face[iedge][2];
          if( aCell[jca0].isHavingChild_Edge(jedge) ){ continue; }
          for(int ien=0;ien<2;++ien){
            int icch0 = aCell[icell].aIC_Cld[ edgeHex[iedge][ien] ];
            if( icch0 < 0 ) continue;
            if( aCell[icch0].isHavingChild_Edge(iedge) ){
              std::cout << "somehting is wrong" << std::endl;
            }
          }
        }
      }
      { // kface
        const int kca0 = aCell[icell].aIC_Adj[ edge2Face[iedge][1] ];
        if( kca0 >= 0 ){
          const int kedge = edge2Face[iedge][3];
          if( aCell[kca0].isHavingChild_Edge(kedge) ){ continue; }
          for(int ien=0;ien<2;++ien){
            int icch0 = aCell[icell].aIC_Cld[ edgeHex[iedge][ien] ];
            if( icch0 < 0 ) continue;
            if( aCell[icch0].isHavingChild_Edge(iedge) ){
              std::cout << "somehting is wrong" << std::endl;
            }
          }
        }
      }
    }
  }
  std::cout << "check finished " << aCell.size() << std::endl;
}
 */

DFM2_INLINE void addEdgeFacePoints(
    std::vector<CPointLattice> &aPoint,
    std::vector<CCell> &aCell,
    //
    const CInput_IsosurfaceStuffing &input) {
  std::vector<int> orderCell;
  orderCell.reserve(aCell.size());
  {
    for (int ilevel = 0;; ++ilevel) {
      int icnt = 0;
      for (int icell = 0; icell < (int) aCell.size(); ++icell) {
        if (aCell[icell].ilevel == ilevel) {
          orderCell.push_back(icell);
          icnt++;
        }
      }
      if (icnt == 0) break;
    }
  }

  ////////////////////////////////////////////////////////////////////////////

  // this need to run coarse to fine
  for (int iicell = 0; iicell < (int) aCell.size(); ++iicell) {
    const int icell = orderCell[iicell];
    //  for(int icell=0;icell<aCell.size();++icell){
    CCell &c = aCell[icell];
    const int icc0 = c.aIC_Cld[0];
    const int icc1 = c.aIC_Cld[1];
    const int icc2 = c.aIC_Cld[2];
    const int icc3 = c.aIC_Cld[3];
    const int icc4 = c.aIC_Cld[4];
    const int icc5 = c.aIC_Cld[5];
    const int icc6 = c.aIC_Cld[6];
    const int icc7 = c.aIC_Cld[7];
    //// set children's corner point that shared with parents
    if (icc0 >= 0) {
      aCell[icc0].aIP[0] = c.aIP[0];
      aCell[icc0].aIP[7] = c.aIP[26];
    }
    if (icc1 >= 0) {
      aCell[icc1].aIP[1] = c.aIP[1];
      aCell[icc1].aIP[6] = c.aIP[26];
    }
    if (icc2 >= 0) {
      aCell[icc2].aIP[2] = c.aIP[2];
      aCell[icc2].aIP[5] = c.aIP[26];
    }
    if (icc3 >= 0) {
      aCell[icc3].aIP[3] = c.aIP[3];
      aCell[icc3].aIP[4] = c.aIP[26];
    }
    if (icc4 >= 0) {
      aCell[icc4].aIP[4] = c.aIP[4];
      aCell[icc4].aIP[3] = c.aIP[26];
    }
    if (icc5 >= 0) {
      aCell[icc5].aIP[5] = c.aIP[5];
      aCell[icc5].aIP[2] = c.aIP[26];
    }
    if (icc6 >= 0) {
      aCell[icc6].aIP[6] = c.aIP[6];
      aCell[icc6].aIP[1] = c.aIP[26];
    }
    if (icc7 >= 0) {
      aCell[icc7].aIP[7] = c.aIP[7];
      aCell[icc7].aIP[0] = c.aIP[26];
    }
    ////
    const int icax = c.aIC_Adj[0];
    const int icaX = c.aIC_Adj[1];
    const int icay = c.aIC_Adj[2];
    const int icaY = c.aIC_Adj[3];
    const int icaz = c.aIC_Adj[4];
    const int icaZ = c.aIC_Adj[5];
    ////
    if (c.aIP[8] == -1) { // edge point (8)
      const int icayc2 = (icay >= 0) ? aCell[icay].aIC_Cld[2] : -1;
      const int icayc3 = (icay >= 0) ? aCell[icay].aIC_Cld[3] : -1;
      const int icazc4 = (icaz >= 0) ? aCell[icaz].aIC_Cld[4] : -1;
      const int icazc5 = (icaz >= 0) ? aCell[icaz].aIC_Cld[5] : -1;
      int icayz = -1, icayzc6 = -1, icayzc7 = -1;
      if (icay >= 0 && icaz >= 0) {
        icayz = aCell[icay].aIC_Adj[4];
        if (icayz >= 0) {
          assert(icayz == aCell[icaz].aIC_Adj[2]);
          icayzc6 = aCell[icayz].aIC_Cld[6];
          icayzc7 = aCell[icayz].aIC_Cld[7];
        }
      }
      if (icc0 >= 0 || icc1 >= 0 || icayc2 >= 0 || icayc3 >= 0 || icazc4 >= 0 || icazc5 >= 0 || icayzc6 >= 0
          || icayzc7 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0];
        const double y0 = aPoint[c.aIP[26]].pos[1] - c.size * 0.5;
        const double z0 = aPoint[c.aIP[26]].pos[2] - c.size * 0.5;
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[8] = ip0;
        if (icay >= 0) { aCell[icay].aIP[9] = ip0; }
        if (icaz >= 0) { aCell[icaz].aIP[10] = ip0; }
        if (icayz >= 0) { aCell[icayz].aIP[11] = ip0; }
        if (icc0 >= 0) { aCell[icc0].aIP[1] = ip0; }
        if (icc1 >= 0) { aCell[icc1].aIP[0] = ip0; }
        if (icayc2 >= 0) { aCell[icayc2].aIP[3] = ip0; }
        if (icayc3 >= 0) { aCell[icayc3].aIP[2] = ip0; }
        if (icazc4 >= 0) { aCell[icazc4].aIP[5] = ip0; }
        if (icazc5 >= 0) { aCell[icazc5].aIP[4] = ip0; }
        if (icayzc6 >= 0) { aCell[icayzc6].aIP[7] = ip0; }
        if (icayzc7 >= 0) { aCell[icayzc7].aIP[6] = ip0; }
      }
    }
    if (c.aIP[9] == -1) { // edge point (20)
      int icaYc0 = (icaY >= 0) ? aCell[icaY].aIC_Cld[0] : -1;
      int icaYc1 = (icaY >= 0) ? aCell[icaY].aIC_Cld[1] : -1;
      int icazc6 = (icaz >= 0) ? aCell[icaz].aIC_Cld[6] : -1;
      int icazc7 = (icaz >= 0) ? aCell[icaz].aIC_Cld[7] : -1;
      int icaYz = -1, icaYzc4 = -1, icaYzc5 = -1;
      if (icaY >= 0 && icaz >= 0) {
        icaYz = aCell[icaY].aIC_Adj[4];
        if (icaYz >= 0) {
          assert(icaYz == aCell[icaz].aIC_Adj[3]);
          icaYzc4 = aCell[icaYz].aIC_Cld[4];
          icaYzc5 = aCell[icaYz].aIC_Cld[5];
        }
      }
      if (icc2 >= 0 || icc3 >= 0 || icaYc0 >= 0 || icaYc1 >= 0 || icazc6 >= 0 || icazc7 >= 0 || icaYzc4 >= 0
          || icaYzc5 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0];
        const double y0 = aPoint[c.aIP[26]].pos[1] + c.size * 0.5;
        const double z0 = aPoint[c.aIP[26]].pos[2] - c.size * 0.5;
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[9] = ip0;
        if (icaY >= 0) { aCell[icaY].aIP[8] = ip0; }
        if (icaz >= 0) { aCell[icaz].aIP[11] = ip0; }
        if (icaYz >= 0) { aCell[icaYz].aIP[10] = ip0; }
        if (icc2 >= 0) { aCell[icc2].aIP[3] = ip0; }
        if (icc3 >= 0) { aCell[icc3].aIP[2] = ip0; }
        if (icaYc0 >= 0) { aCell[icaYc0].aIP[1] = ip0; }
        if (icaYc1 >= 0) { aCell[icaYc1].aIP[0] = ip0; }
        if (icazc6 >= 0) { aCell[icazc6].aIP[7] = ip0; }
        if (icazc7 >= 0) { aCell[icazc7].aIP[6] = ip0; }
        if (icaYzc4 >= 0) { aCell[icaYzc4].aIP[5] = ip0; }
        if (icaYzc5 >= 0) { aCell[icaYzc5].aIP[4] = ip0; }
      }
    }
    if (c.aIP[10] == -1) { // edge point (20)
      int icayc6 = (icay >= 0) ? aCell[icay].aIC_Cld[6] : -1;
      int icayc7 = (icay >= 0) ? aCell[icay].aIC_Cld[7] : -1;
      int icaZc0 = (icaZ >= 0) ? aCell[icaZ].aIC_Cld[0] : -1;
      int icaZc1 = (icaZ >= 0) ? aCell[icaZ].aIC_Cld[1] : -1;
      int icayZ = -1, icayZc2 = -1, icayZc3 = -1;
      if (icay >= 0 && icaZ >= 0) {
        icayZ = aCell[icay].aIC_Adj[5];
        if (icayZ >= 0) {
          assert(icayZ == aCell[icaZ].aIC_Adj[2]);
          icayZc2 = aCell[icayZ].aIC_Cld[2];
          icayZc3 = aCell[icayZ].aIC_Cld[3];
        }
      }
      if (icc4 >= 0 || icc5 >= 0 || icayc6 >= 0 || icayc7 >= 0 || icaZc0 >= 0 || icaZc1 >= 0 || icayZc2 >= 0
          || icayZc3 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0];
        const double y0 = aPoint[c.aIP[26]].pos[1] - c.size * 0.5;
        const double z0 = aPoint[c.aIP[26]].pos[2] + c.size * 0.5;
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[10] = ip0;
        if (icay >= 0) { aCell[icay].aIP[11] = ip0; }
        if (icaZ >= 0) { aCell[icaZ].aIP[8] = ip0; }
        if (icayZ >= 0) { aCell[icayZ].aIP[9] = ip0; }
        if (icc4 >= 0) { aCell[icc4].aIP[5] = ip0; }
        if (icc5 >= 0) { aCell[icc5].aIP[4] = ip0; }
        if (icayc6 >= 0) { aCell[icayc6].aIP[7] = ip0; }
        if (icayc7 >= 0) { aCell[icayc7].aIP[6] = ip0; }
        if (icaZc0 >= 0) { aCell[icaZc0].aIP[1] = ip0; }
        if (icaZc1 >= 0) { aCell[icaZc1].aIP[0] = ip0; }
        if (icayZc2 >= 0) { aCell[icayZc2].aIP[3] = ip0; }
        if (icayZc3 >= 0) { aCell[icayZc3].aIP[2] = ip0; }
      }
    }
    if (c.aIP[11] == -1) { // edge point (20)
      int icaYc4 = (icaY >= 0) ? aCell[icaY].aIC_Cld[4] : -1;
      int icaYc5 = (icaY >= 0) ? aCell[icaY].aIC_Cld[5] : -1;
      int icaZc2 = (icaZ >= 0) ? aCell[icaZ].aIC_Cld[2] : -1;
      int icaZc3 = (icaZ >= 0) ? aCell[icaZ].aIC_Cld[3] : -1;
      int icaYZ = -1, icaYZc0 = -1, icaYZc1 = -1;
      if (icaY >= 0 && icaZ >= 0) {
        icaYZ = aCell[icaY].aIC_Adj[5];
        if (icaYZ >= 0) {
          assert(icaYZ == aCell[icaZ].aIC_Adj[3]);
          icaYZc0 = aCell[icaYZ].aIC_Cld[0];
          icaYZc1 = aCell[icaYZ].aIC_Cld[1];
        }
      }
      if (icc6 >= 0 || icc7 >= 0 || icaYc4 >= 0 || icaYc5 >= 0 || icaZc2 >= 0 || icaZc3 >= 0 || icaYZc0 >= 0
          || icaYZc1 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0];
        const double y0 = aPoint[c.aIP[26]].pos[1] + c.size * 0.5;
        const double z0 = aPoint[c.aIP[26]].pos[2] + c.size * 0.5;
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[11] = ip0;
        if (icaY >= 0) { aCell[icaY].aIP[10] = ip0; }
        if (icaZ >= 0) { aCell[icaZ].aIP[9] = ip0; }
        if (icaYZ >= 0) { aCell[icaYZ].aIP[8] = ip0; }
        if (icc6 >= 0) { aCell[icc6].aIP[7] = ip0; }
        if (icc7 >= 0) { aCell[icc7].aIP[6] = ip0; }
        if (icaYc4 >= 0) { aCell[icaYc4].aIP[5] = ip0; }
        if (icaYc5 >= 0) { aCell[icaYc5].aIP[4] = ip0; }
        if (icaZc2 >= 0) { aCell[icaZc2].aIP[3] = ip0; }
        if (icaZc3 >= 0) { aCell[icaZc3].aIP[2] = ip0; }
        if (icaYZc0 >= 0) { aCell[icaYZc0].aIP[1] = ip0; }
        if (icaYZc1 >= 0) { aCell[icaYZc1].aIP[0] = ip0; }
      }
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (c.aIP[12] == -1) { // edge point (12)
      int icaxc1 = (icax >= 0) ? aCell[icax].aIC_Cld[1] : -1;
      int icaxc3 = (icax >= 0) ? aCell[icax].aIC_Cld[3] : -1;
      int icazc4 = (icaz >= 0) ? aCell[icaz].aIC_Cld[4] : -1;
      int icazc6 = (icaz >= 0) ? aCell[icaz].aIC_Cld[6] : -1;
      int icaxz = -1, icaxzc5 = -1, icaxzc7 = -1;
      if (icax >= 0 && icaz >= 0) {
        icaxz = aCell[icax].aIC_Adj[4];
        if (icaxz >= 0) {
          assert(icaxz == aCell[icaz].aIC_Adj[0]);
          icaxzc5 = aCell[icaxz].aIC_Cld[5];
          icaxzc7 = aCell[icaxz].aIC_Cld[7];
        }
      }
      if (icc0 >= 0 || icc2 >= 0 || icaxc1 >= 0 || icaxc3 >= 0 || icazc4 >= 0 || icazc6 >= 0 || icaxzc5 >= 0
          || icaxzc7 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0] - c.size * 0.5;
        const double y0 = aPoint[c.aIP[26]].pos[1];
        const double z0 = aPoint[c.aIP[26]].pos[2] - c.size * 0.5;
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[12] = ip0;
        if (icax >= 0) { aCell[icax].aIP[13] = ip0; }
        if (icaz >= 0) { aCell[icaz].aIP[14] = ip0; }
        if (icaxz >= 0) { aCell[icaxz].aIP[15] = ip0; }
        if (icc0 >= 0) { aCell[icc0].aIP[2] = ip0; }
        if (icc2 >= 0) { aCell[icc2].aIP[0] = ip0; }
        if (icaxc1 >= 0) { aCell[icaxc1].aIP[3] = ip0; }
        if (icaxc3 >= 0) { aCell[icaxc3].aIP[1] = ip0; }
        if (icazc4 >= 0) { aCell[icazc4].aIP[6] = ip0; }
        if (icazc6 >= 0) { aCell[icazc6].aIP[4] = ip0; }
        if (icaxzc5 >= 0) { aCell[icaxzc5].aIP[7] = ip0; }
        if (icaxzc7 >= 0) { aCell[icaxzc7].aIP[5] = ip0; }
      }
    }
    if (c.aIP[13] == -1) { // edge point (12)
      int icaXc0 = (icaX >= 0) ? aCell[icaX].aIC_Cld[0] : -1;
      int icaXc2 = (icaX >= 0) ? aCell[icaX].aIC_Cld[2] : -1;
      int icazc5 = (icaz >= 0) ? aCell[icaz].aIC_Cld[5] : -1;
      int icazc7 = (icaz >= 0) ? aCell[icaz].aIC_Cld[7] : -1;
      int icaXz = -1, icaXzc4 = -1, icaXzc6 = -1;
      if (icaX >= 0 && icaz >= 0) {
        icaXz = aCell[icaX].aIC_Adj[4];
        if (icaXz >= 0) {
          assert(icaXz == aCell[icaz].aIC_Adj[1]);
          icaXzc4 = aCell[icaXz].aIC_Cld[4];
          icaXzc6 = aCell[icaXz].aIC_Cld[6];
        }
      }
      if (icc1 >= 0 || icc3 >= 0 || icaXc0 >= 0 || icaXc2 >= 0 || icazc5 >= 0 || icazc7 >= 0 || icaXzc4 >= 0
          || icaXzc6 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0] + c.size * 0.5;
        const double y0 = aPoint[c.aIP[26]].pos[1];
        const double z0 = aPoint[c.aIP[26]].pos[2] - c.size * 0.5;
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[13] = ip0;
        if (icaX >= 0) { aCell[icaX].aIP[12] = ip0; }
        if (icaz >= 0) { aCell[icaz].aIP[15] = ip0; }
        if (icaXz >= 0) { aCell[icaXz].aIP[14] = ip0; }
        if (icc1 >= 0) { aCell[icc1].aIP[3] = ip0; }
        if (icc3 >= 0) { aCell[icc3].aIP[1] = ip0; }
        if (icaXc0 >= 0) { aCell[icaXc0].aIP[2] = ip0; }
        if (icaXc2 >= 0) { aCell[icaXc2].aIP[0] = ip0; }
        if (icazc5 >= 0) { aCell[icazc5].aIP[7] = ip0; }
        if (icazc7 >= 0) { aCell[icazc7].aIP[5] = ip0; }
        if (icaXzc4 >= 0) { aCell[icaXzc4].aIP[6] = ip0; }
        if (icaXzc6 >= 0) { aCell[icaXzc6].aIP[4] = ip0; }
      }
    }
    if (c.aIP[14] == -1) { // edge point (14)
      int icaxc5 = (icax >= 0) ? aCell[icax].aIC_Cld[5] : -1;
      int icaxc7 = (icax >= 0) ? aCell[icax].aIC_Cld[7] : -1;
      int icaZc0 = (icaZ >= 0) ? aCell[icaZ].aIC_Cld[0] : -1;
      int icaZc2 = (icaZ >= 0) ? aCell[icaZ].aIC_Cld[2] : -1;
      int icaxZ = -1, icaxZc1 = -1, icaxZc3 = -1;
      if (icax >= 0 && icaZ >= 0) {
        icaxZ = aCell[icax].aIC_Adj[5];
        if (icaxZ >= 0) {
          assert(icaxZ == aCell[icaZ].aIC_Adj[0]);
          icaxZc1 = aCell[icaxZ].aIC_Cld[1];
          icaxZc3 = aCell[icaxZ].aIC_Cld[3];
        }
      }
      if (icc4 >= 0 || icc6 >= 0 || icaxc5 >= 0 || icaxc7 >= 0 || icaZc0 >= 0 || icaZc2 >= 0 || icaxZc1 >= 0
          || icaxZc3 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0] - c.size * 0.5;
        const double y0 = aPoint[c.aIP[26]].pos[1];
        const double z0 = aPoint[c.aIP[26]].pos[2] + c.size * 0.5;
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[14] = ip0;
        if (icax >= 0) { aCell[icax].aIP[15] = ip0; }
        if (icaZ >= 0) { aCell[icaZ].aIP[12] = ip0; }
        if (icaxZ >= 0) { aCell[icaxZ].aIP[13] = ip0; }
        if (icc4 >= 0) { aCell[icc4].aIP[6] = ip0; }
        if (icc6 >= 0) { aCell[icc6].aIP[4] = ip0; }
        if (icaxc5 >= 0) { aCell[icaxc5].aIP[7] = ip0; }
        if (icaxc7 >= 0) { aCell[icaxc7].aIP[5] = ip0; }
        if (icaZc0 >= 0) { aCell[icaZc0].aIP[2] = ip0; }
        if (icaZc2 >= 0) { aCell[icaZc2].aIP[0] = ip0; }
        if (icaxZc1 >= 0) { aCell[icaxZc1].aIP[3] = ip0; }
        if (icaxZc3 >= 0) { aCell[icaxZc3].aIP[1] = ip0; }
      }
    }
    if (c.aIP[15] == -1) { // edge point (15)
      int icaXc4 = (icaX >= 0) ? aCell[icaX].aIC_Cld[4] : -1;
      int icaXc6 = (icaX >= 0) ? aCell[icaX].aIC_Cld[6] : -1;
      int icaZc1 = (icaZ >= 0) ? aCell[icaZ].aIC_Cld[1] : -1;
      int icaZc3 = (icaZ >= 0) ? aCell[icaZ].aIC_Cld[3] : -1;
      int icaXZ = -1, icaXZc0 = -1, icaXZc2 = -1;
      if (icaX >= 0 && icaZ >= 0) {
        icaXZ = aCell[icaX].aIC_Adj[5];
        if (icaXZ >= 0) {
          assert(icaXZ == aCell[icaZ].aIC_Adj[1]);
          icaXZc0 = aCell[icaXZ].aIC_Cld[0];
          icaXZc2 = aCell[icaXZ].aIC_Cld[2];
        }
      }
      if (icc5 >= 0 || icc7 >= 0 || icaXc4 >= 0 || icaXc6 >= 0 || icaZc1 >= 0 || icaZc3 >= 0 || icaXZc0 >= 0
          || icaXZc2 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0] + c.size * 0.5;
        const double y0 = aPoint[c.aIP[26]].pos[1];
        const double z0 = aPoint[c.aIP[26]].pos[2] + c.size * 0.5;
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[15] = ip0;
        if (icaX >= 0) { aCell[icaX].aIP[14] = ip0; }
        if (icaZ >= 0) { aCell[icaZ].aIP[13] = ip0; }
        if (icaXZ >= 0) { aCell[icaXZ].aIP[12] = ip0; }
        if (icc5 >= 0) { aCell[icc5].aIP[7] = ip0; }
        if (icc7 >= 0) { aCell[icc7].aIP[5] = ip0; }
        if (icaXc4 >= 0) { aCell[icaXc4].aIP[6] = ip0; }
        if (icaXc6 >= 0) { aCell[icaXc6].aIP[4] = ip0; }
        if (icaZc1 >= 0) { aCell[icaZc1].aIP[3] = ip0; }
        if (icaZc3 >= 0) { aCell[icaZc3].aIP[1] = ip0; }
        if (icaXZc0 >= 0) { aCell[icaXZc0].aIP[2] = ip0; }
        if (icaXZc2 >= 0) { aCell[icaXZc2].aIP[0] = ip0; }
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (c.aIP[16] == -1) { // edge point (16)
      int icaxc1 = (icax >= 0) ? aCell[icax].aIC_Cld[1] : -1;
      int icaxc5 = (icax >= 0) ? aCell[icax].aIC_Cld[5] : -1;
      int icayc2 = (icay >= 0) ? aCell[icay].aIC_Cld[2] : -1;
      int icayc6 = (icay >= 0) ? aCell[icay].aIC_Cld[6] : -1;
      int icaxy = -1, icaxyc3 = -1, icaxyc7 = -1;
      if (icax >= 0 && icay >= 0) {
        icaxy = aCell[icax].aIC_Adj[2];
        if (icaxy >= 0) {
          assert(icaxy == aCell[icay].aIC_Adj[0]);
          icaxyc3 = aCell[icaxy].aIC_Cld[3];
          icaxyc7 = aCell[icaxy].aIC_Cld[7];
        }
      }
      if (icc0 >= 0 || icc4 >= 0 || icaxc1 >= 0 || icaxc5 >= 0 || icayc2 >= 0 || icayc6 >= 0 || icaxyc3 >= 0
          || icaxyc7 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0] - c.size * 0.5;
        const double y0 = aPoint[c.aIP[26]].pos[1] - c.size * 0.5;
        const double z0 = aPoint[c.aIP[26]].pos[2];
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[16] = ip0;
        if (icax >= 0) { aCell[icax].aIP[17] = ip0; }
        if (icay >= 0) { aCell[icay].aIP[18] = ip0; }
        if (icaxy >= 0) { aCell[icaxy].aIP[19] = ip0; }
        if (icc0 >= 0) { aCell[icc0].aIP[4] = ip0; }
        if (icc4 >= 0) { aCell[icc4].aIP[0] = ip0; }
        if (icaxc1 >= 0) { aCell[icaxc1].aIP[5] = ip0; }
        if (icaxc5 >= 0) { aCell[icaxc5].aIP[1] = ip0; }
        if (icayc2 >= 0) { aCell[icayc2].aIP[6] = ip0; }
        if (icayc6 >= 0) { aCell[icayc6].aIP[2] = ip0; }
        if (icaxyc3 >= 0) { aCell[icaxyc3].aIP[7] = ip0; }
        if (icaxyc7 >= 0) { aCell[icaxyc7].aIP[3] = ip0; }
      }
    }
    if (c.aIP[17] == -1) { // edge point (17)
      int icaXc0 = (icaX >= 0) ? aCell[icaX].aIC_Cld[0] : -1;
      int icaXc4 = (icaX >= 0) ? aCell[icaX].aIC_Cld[4] : -1;
      int icayc3 = (icay >= 0) ? aCell[icay].aIC_Cld[3] : -1;
      int icayc7 = (icay >= 0) ? aCell[icay].aIC_Cld[7] : -1;
      int icaXy = -1, icaXyc2 = -1, icaXyc6 = -1;
      if (icaX >= 0 && icay >= 0) {
        icaXy = aCell[icaX].aIC_Adj[2];
        if (icaXy >= 0) {
          assert(icaXy == aCell[icay].aIC_Adj[1]);
          icaXyc2 = aCell[icaXy].aIC_Cld[2];
          icaXyc6 = aCell[icaXy].aIC_Cld[6];
        }
      }
      if (icc1 >= 0 || icc5 >= 0 || icaXc0 >= 0 || icaXc4 >= 0 || icayc3 >= 0 || icayc7 >= 0 || icaXyc2 >= 0
          || icaXyc6 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0] + c.size * 0.5;
        const double y0 = aPoint[c.aIP[26]].pos[1] - c.size * 0.5;
        const double z0 = aPoint[c.aIP[26]].pos[2];
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[17] = ip0;
        if (icaX >= 0) { aCell[icaX].aIP[16] = ip0; }
        if (icay >= 0) { aCell[icay].aIP[19] = ip0; }
        if (icaXy >= 0) { aCell[icaXy].aIP[18] = ip0; }
        if (icc1 >= 0) { aCell[icc1].aIP[5] = ip0; }
        if (icc5 >= 0) { aCell[icc5].aIP[1] = ip0; }
        if (icaXc0 >= 0) { aCell[icaXc0].aIP[4] = ip0; }
        if (icaXc4 >= 0) { aCell[icaXc4].aIP[0] = ip0; }
        if (icayc3 >= 0) { aCell[icayc3].aIP[7] = ip0; }
        if (icayc7 >= 0) { aCell[icayc7].aIP[3] = ip0; }
        if (icaXyc2 >= 0) { aCell[icaXyc2].aIP[6] = ip0; }
        if (icaXyc6 >= 0) { aCell[icaXyc6].aIP[2] = ip0; }
      }
    }
    if (c.aIP[18] == -1) { // edge point (16)
      int icaxc3 = (icax >= 0) ? aCell[icax].aIC_Cld[3] : -1;
      int icaxc7 = (icax >= 0) ? aCell[icax].aIC_Cld[7] : -1;
      int icaYc0 = (icaY >= 0) ? aCell[icaY].aIC_Cld[0] : -1;
      int icaYc4 = (icaY >= 0) ? aCell[icaY].aIC_Cld[4] : -1;
      int icaxY = -1, icaxYc1 = -1, icaxYc5 = -1;
      if (icax >= 0 && icaY >= 0) {
        icaxY = aCell[icax].aIC_Adj[3];
        if (icaxY >= 0) {
          assert(icaxY == aCell[icaY].aIC_Adj[0]);
          icaxYc1 = aCell[icaxY].aIC_Cld[1];
          icaxYc5 = aCell[icaxY].aIC_Cld[5];
        }
      }
      if (icc2 >= 0 || icc6 >= 0 || icaxc3 >= 0 || icaxc7 >= 0 || icaYc0 >= 0 || icaYc4 >= 0 || icaxYc1 >= 0
          || icaxYc5 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0] - c.size * 0.5;
        const double y0 = aPoint[c.aIP[26]].pos[1] + c.size * 0.5;
        const double z0 = aPoint[c.aIP[26]].pos[2];
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[18] = ip0;
        if (icax >= 0) { aCell[icax].aIP[19] = ip0; }
        if (icaY >= 0) { aCell[icaY].aIP[16] = ip0; }
        if (icaxY >= 0) { aCell[icaxY].aIP[17] = ip0; }
        if (icc2 >= 0) { aCell[icc2].aIP[6] = ip0; }
        if (icc6 >= 0) { aCell[icc6].aIP[2] = ip0; }
        if (icaxc3 >= 0) { aCell[icaxc3].aIP[7] = ip0; }
        if (icaxc7 >= 0) { aCell[icaxc7].aIP[3] = ip0; }
        if (icaYc0 >= 0) { aCell[icaYc0].aIP[4] = ip0; }
        if (icaYc4 >= 0) { aCell[icaYc4].aIP[0] = ip0; }
        if (icaxYc1 >= 0) { aCell[icaxYc1].aIP[5] = ip0; }
        if (icaxYc5 >= 0) { aCell[icaxYc5].aIP[1] = ip0; }
      }
    }

    if (c.aIP[19] == -1) { // edge point (19)
      int icaXc2 = (icaX >= 0) ? aCell[icaX].aIC_Cld[2] : -1;
      int icaXc6 = (icaX >= 0) ? aCell[icaX].aIC_Cld[6] : -1;
      int icaYc1 = (icaY >= 0) ? aCell[icaY].aIC_Cld[1] : -1;
      int icaYc5 = (icaY >= 0) ? aCell[icaY].aIC_Cld[5] : -1;
      int icaXY = -1, icaXYc0 = -1, icaXYc4 = -1;
      if (icaX >= 0 && icaY >= 0) {
        icaXY = aCell[icaX].aIC_Adj[3];
        if (icaXY >= 0) {
          assert(icaXY == aCell[icaY].aIC_Adj[1]);
          icaXYc0 = aCell[icaXY].aIC_Cld[0];
          icaXYc4 = aCell[icaXY].aIC_Cld[4];
        }
      }
      if (icc3 >= 0 || icc7 >= 0 || icaXc2 >= 0 || icaXc6 >= 0 || icaYc1 >= 0 || icaYc5 >= 0 || icaXYc0 >= 0
          || icaXYc4 >= 0) {
        const int ip0 = (int) aPoint.size();
        const double x0 = aPoint[c.aIP[26]].pos[0] + c.size * 0.5;
        const double y0 = aPoint[c.aIP[26]].pos[1] + c.size * 0.5;
        const double z0 = aPoint[c.aIP[26]].pos[2];
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[19] = ip0;
        if (icaX >= 0) { aCell[icaX].aIP[18] = ip0; }
        if (icaY >= 0) { aCell[icaY].aIP[17] = ip0; }
        if (icaXY >= 0) { aCell[icaXY].aIP[16] = ip0; }
        if (icc3 >= 0) { aCell[icc3].aIP[7] = ip0; }
        if (icc7 >= 0) { aCell[icc7].aIP[3] = ip0; }
        if (icaXc2 >= 0) { aCell[icaXc2].aIP[6] = ip0; }
        if (icaXc6 >= 0) { aCell[icaXc6].aIP[2] = ip0; }
        if (icaYc1 >= 0) { aCell[icaYc1].aIP[5] = ip0; }
        if (icaYc5 >= 0) { aCell[icaYc5].aIP[1] = ip0; }
        if (icaXYc0 >= 0) { aCell[icaXYc0].aIP[4] = ip0; }
        if (icaXYc4 >= 0) { aCell[icaXYc4].aIP[0] = ip0; }
      }
    }
    // ------------------------------------
    if (c.aIP[20] == -1) { // face point (20)
      int icaxc1 = -1, icaxc3 = -1, icaxc5 = -1, icaxc7 = -1;
      if (icax >= 0) {
        CCell &ca = aCell[icax];
        icaxc1 = ca.aIC_Cld[1];
        icaxc3 = ca.aIC_Cld[3];
        icaxc5 = ca.aIC_Cld[5];
        icaxc7 = ca.aIC_Cld[7];
      }
      if (icc0 >= 0 || icc2 >= 0 || icc4 >= 0 || icc6 >= 0 || icaxc1 >= 0 || icaxc3 >= 0 || icaxc5 >= 0
          || icaxc7 >= 0) {
        int ip0 = (int) aPoint.size();
        double x0 = aPoint[c.aIP[26]].pos[0] - c.size * 0.5;
        double y0 = aPoint[c.aIP[26]].pos[1];
        double z0 = aPoint[c.aIP[26]].pos[2];
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[20] = ip0;
        if (icax >= 0) { aCell[icax].aIP[21] = ip0; }
        if (icc0 >= 0) { aCell[icc0].aIP[6] = ip0; }
        if (icc2 >= 0) { aCell[icc2].aIP[4] = ip0; }
        if (icc4 >= 0) { aCell[icc4].aIP[2] = ip0; }
        if (icc6 >= 0) { aCell[icc6].aIP[0] = ip0; }
        if (icaxc1 >= 0) { aCell[icaxc1].aIP[7] = ip0; }
        if (icaxc3 >= 0) { aCell[icaxc3].aIP[5] = ip0; }
        if (icaxc5 >= 0) { aCell[icaxc5].aIP[3] = ip0; }
        if (icaxc7 >= 0) { aCell[icaxc7].aIP[1] = ip0; }
      }
    }
    if (c.aIP[21] == -1) { // face point (21)
      int icaXc0 = -1, icaXc2 = -1, icaXc4 = -1, icaXc6 = -1;
      if (icaX >= 0) {
        CCell &ca = aCell[icaX];
        icaXc0 = ca.aIC_Cld[0];
        icaXc2 = ca.aIC_Cld[2];
        icaXc4 = ca.aIC_Cld[4];
        icaXc6 = ca.aIC_Cld[6];
      }
      if (icc0 >= 0 || icc2 >= 0 || icc4 >= 0 || icc6 >= 0 || icaXc0 >= 0 || icaXc2 >= 0 || icaXc4 >= 0
          || icaXc6 >= 0) {
        int ip0 = (int) aPoint.size();
        double x0 = aPoint[c.aIP[26]].pos[0] + c.size * 0.5;
        double y0 = aPoint[c.aIP[26]].pos[1];
        double z0 = aPoint[c.aIP[26]].pos[2];
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[21] = ip0;
        if (icaX >= 0) { aCell[icaX].aIP[20] = ip0; }
        if (icc1 >= 0) { aCell[icc1].aIP[7] = ip0; }
        if (icc3 >= 0) { aCell[icc3].aIP[5] = ip0; }
        if (icc5 >= 0) { aCell[icc5].aIP[3] = ip0; }
        if (icc7 >= 0) { aCell[icc7].aIP[1] = ip0; }
        if (icaXc0 >= 0) { aCell[icaXc0].aIP[6] = ip0; }
        if (icaXc2 >= 0) { aCell[icaXc2].aIP[4] = ip0; }
        if (icaXc4 >= 0) { aCell[icaXc4].aIP[2] = ip0; }
        if (icaXc6 >= 0) { aCell[icaXc6].aIP[0] = ip0; }
      }
    }
    if (c.aIP[22] == -1) { // face point (22)
      int icayc2 = -1, icayc3 = -1, icayc6 = -1, icayc7 = -1;
      if (icay >= 0) {
        CCell &ca = aCell[icay];
        icayc2 = ca.aIC_Cld[2];
        icayc3 = ca.aIC_Cld[3];
        icayc6 = ca.aIC_Cld[6];
        icayc7 = ca.aIC_Cld[7];
      }
      if (icc0 >= 0 || icc1 >= 0 || icc4 >= 0 || icc5 >= 0 || icayc2 >= 0 || icayc3 >= 0 || icayc6 >= 0
          || icayc7 >= 0) {
        int ip0 = (int) aPoint.size();
        double x0 = aPoint[c.aIP[26]].pos[0];
        double y0 = aPoint[c.aIP[26]].pos[1] - c.size * 0.5;
        double z0 = aPoint[c.aIP[26]].pos[2];
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[22] = ip0;
        if (icay >= 0) { aCell[icay].aIP[23] = ip0; }
        if (icc0 >= 0) { aCell[icc0].aIP[5] = ip0; }
        if (icc1 >= 0) { aCell[icc1].aIP[4] = ip0; }
        if (icc4 >= 0) { aCell[icc4].aIP[1] = ip0; }
        if (icc5 >= 0) { aCell[icc5].aIP[0] = ip0; }
        if (icayc2 >= 0) { aCell[icayc2].aIP[7] = ip0; }
        if (icayc3 >= 0) { aCell[icayc3].aIP[6] = ip0; }
        if (icayc6 >= 0) { aCell[icayc6].aIP[3] = ip0; }
        if (icayc7 >= 0) { aCell[icayc7].aIP[2] = ip0; }
      }
    }
    if (c.aIP[23] == -1) { // face point (23)
      int icaYc0 = -1, icaYc1 = -1, icaYc4 = -1, icaYc5 = -1;
      if (icaY >= 0) {
        CCell &ca = aCell[icaY];
        icaYc0 = ca.aIC_Cld[0];
        icaYc1 = ca.aIC_Cld[1];
        icaYc4 = ca.aIC_Cld[4];
        icaYc5 = ca.aIC_Cld[5];
      }
      if (icc2 >= 0 || icc3 >= 0 || icc6 >= 0 || icc7 >= 0 || icaYc0 >= 0 || icaYc1 >= 0 || icaYc4 >= 0
          || icaYc5 >= 0) {
        int ip0 = (int) aPoint.size();
        double x0 = aPoint[c.aIP[26]].pos[0];
        double y0 = aPoint[c.aIP[26]].pos[1] + c.size * 0.5;
        double z0 = aPoint[c.aIP[26]].pos[2];
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[23] = ip0;
        if (icaY >= 0) { aCell[icaY].aIP[22] = ip0; }
        if (icc2 >= 0) { aCell[icc2].aIP[7] = ip0; }
        if (icc3 >= 0) { aCell[icc3].aIP[6] = ip0; }
        if (icc6 >= 0) { aCell[icc6].aIP[3] = ip0; }
        if (icc7 >= 0) { aCell[icc7].aIP[2] = ip0; }
        if (icaYc0 >= 0) { aCell[icaYc0].aIP[5] = ip0; }
        if (icaYc1 >= 0) { aCell[icaYc1].aIP[4] = ip0; }
        if (icaYc4 >= 0) { aCell[icaYc4].aIP[1] = ip0; }
        if (icaYc5 >= 0) { aCell[icaYc5].aIP[0] = ip0; }
      }
    }
    if (c.aIP[24] == -1) { // face point (24)
      int icac4 = -1, icac5 = -1, icac6 = -1, icac7 = -1;
      if (icaz >= 0) {
        CCell &ca = aCell[icaz];
        icac4 = ca.aIC_Cld[4];
        icac5 = ca.aIC_Cld[5];
        icac6 = ca.aIC_Cld[6];
        icac7 = ca.aIC_Cld[7];
      }
      if (icc0 >= 0 || icc1 >= 0 || icc2 >= 0 || icc3 >= 0 || icac4 >= 0 || icac5 >= 0 || icac6 >= 0 || icac7 >= 0) {
        int ip0 = (int) aPoint.size();
        double x0 = aPoint[c.aIP[26]].pos[0];
        double y0 = aPoint[c.aIP[26]].pos[1];
        double z0 = aPoint[c.aIP[26]].pos[2] - c.size * 0.5;
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[24] = ip0;
        if (icaz >= 0) { aCell[icaz].aIP[25] = ip0; }
        if (icc0 >= 0) { aCell[icc0].aIP[3] = ip0; }
        if (icc1 >= 0) { aCell[icc1].aIP[2] = ip0; }
        if (icc2 >= 0) { aCell[icc2].aIP[1] = ip0; }
        if (icc3 >= 0) { aCell[icc3].aIP[0] = ip0; }
        if (icac4 >= 0) { aCell[icac4].aIP[7] = ip0; }
        if (icac5 >= 0) { aCell[icac5].aIP[6] = ip0; }
        if (icac6 >= 0) { aCell[icac6].aIP[5] = ip0; }
        if (icac7 >= 0) { aCell[icac7].aIP[4] = ip0; }
      }
    }
    if (c.aIP[25] == -1) { // face point (24)
      int icac0 = -1, icac1 = -1, icac2 = -1, icac3 = -1;
      if (icaZ >= 0) {
        CCell &ca = aCell[icaZ];
        icac0 = ca.aIC_Cld[0];
        icac1 = ca.aIC_Cld[1];
        icac2 = ca.aIC_Cld[2];
        icac3 = ca.aIC_Cld[3];
      }
      if (icc4 >= 0 || icc5 >= 0 || icc6 >= 0 || icc7 >= 0 || icac0 >= 0 || icac1 >= 0 || icac2 >= 0 || icac3 >= 0) {
        int ip0 = (int) aPoint.size();
        double x0 = aPoint[c.aIP[26]].pos[0];
        double y0 = aPoint[c.aIP[26]].pos[1];
        double z0 = aPoint[c.aIP[26]].pos[2] + c.size * 0.5;
        aPoint.emplace_back(x0, y0, z0, input.SignedDistance(x0, y0, z0));
        c.aIP[25] = ip0;
        if (icaZ >= 0) { aCell[icaZ].aIP[24] = ip0; }
        if (icc4 >= 0) { aCell[icc4].aIP[7] = ip0; }
        if (icc5 >= 0) { aCell[icc5].aIP[6] = ip0; }
        if (icc6 >= 0) { aCell[icc6].aIP[5] = ip0; }
        if (icc7 >= 0) { aCell[icc7].aIP[4] = ip0; }
        if (icac0 >= 0) { aCell[icac0].aIP[3] = ip0; }
        if (icac1 >= 0) { aCell[icac1].aIP[2] = ip0; }
        if (icac2 >= 0) { aCell[icac2].aIP[1] = ip0; }
        if (icac3 >= 0) { aCell[icac3].aIP[0] = ip0; }
      }
    }
  }
}

DFM2_INLINE void makeTetLattice(
    std::vector<unsigned int> &aTet,
    const std::vector<CCell> &aCell) {
  aTet.clear();
  const int faceHex[6][4] = {
      {0, 4, 6, 2},
      {1, 3, 7, 5},
      {0, 1, 5, 4},
      {2, 6, 7, 3},
      {0, 2, 3, 1},
      {4, 5, 7, 6}};
  const int halfPyramidPtn[8][6] = {
      {0, 0, 0, 0, 0, 0},
      {0, 0, 1, 1, 1, 1},
      {1, 1, 0, 0, 1, 1},
      {1, 1, 1, 1, 0, 0},
      {1, 1, 1, 1, 0, 0},
      {1, 1, 0, 0, 1, 1},
      {0, 0, 1, 1, 1, 1},
      {0, 0, 0, 0, 0, 0}};
  const int edgePointFace[6][4] = {
      {16, 14, 18, 12},
      {13, 19, 15, 17},
      {8, 17, 10, 16},
      {18, 11, 19, 9},
      {12, 9, 13, 8},
      {10, 15, 11, 14}};
  for (int ic = 0; ic < (int) aCell.size(); ++ic) {
    const CCell &c = aCell[ic];
    const int ip_c = c.aIP[26]; // center
    for (int iface = 0; iface < 6; ++iface) {
      int ip_d = c.aIP[20 + iface]; // face point
      if (ip_d == -1) { // no face point
        int ica0 = c.aIC_Adj[iface];
        if (ica0 < 0) { // facing outer or intternal parent hex
          int ip_f0 = c.aIP[faceHex[iface][0]];
          int ip_f1 = c.aIP[faceHex[iface][1]];
          int ip_f2 = c.aIP[faceHex[iface][2]];
          int ip_f3 = c.aIP[faceHex[iface][3]];
          int iptn = 0;
          if (c.iparent >= 0) {
            iptn = halfPyramidPtn[c.iparent_pos][iface];
          }
          if (iptn == 0) { // connect 0-2
            aTet.push_back(ip_c);
            aTet.push_back(ip_f0);
            aTet.push_back(ip_f1);
            aTet.push_back(ip_f2); // half pyramid
            aTet.push_back(ip_c);
            aTet.push_back(ip_f2);
            aTet.push_back(ip_f3);
            aTet.push_back(ip_f0); // half pyramid
          } else { // connect 1-3
            aTet.push_back(ip_c);
            aTet.push_back(ip_f0);
            aTet.push_back(ip_f1);
            aTet.push_back(ip_f3); // half pyramid
            aTet.push_back(ip_c);
            aTet.push_back(ip_f1);
            aTet.push_back(ip_f2);
            aTet.push_back(ip_f3); // half pyramid
          }
          continue;
        }
        const CCell &ca = aCell[ica0];
        if (ca.ilevel != c.ilevel) {
          std::cout << ca.ilevel << " " << c.ilevel << std::endl;
          continue;
        }
        int ip_ca = ca.aIP[26]; // adjacent center
        if (ip_c < ip_ca) continue; // avoid duplicate
        for (int ieface = 0; ieface < 4; ++ieface) {
          int ip_e0 = c.aIP[faceHex[iface][(ieface + 0) % 4]];
          int ip_em = c.aIP[edgePointFace[iface][ieface]];
          int ip_e1 = c.aIP[faceHex[iface][(ieface + 1) % 4]];
          if (ip_em == -1) {
            aTet.push_back(ip_c);
            aTet.push_back(ip_ca);
            aTet.push_back(ip_e0);
            aTet.push_back(ip_e1); // bcc
          } else {
            aTet.push_back(ip_c);
            aTet.push_back(ip_ca);
            aTet.push_back(ip_e0);
            aTet.push_back(ip_em); // bcc bisect
            aTet.push_back(ip_c);
            aTet.push_back(ip_ca);
            aTet.push_back(ip_em);
            aTet.push_back(ip_e1); // bcc bisect
          }
        }
      } else { // ip_d != -1
        for (int ieface = 0; ieface < 4; ++ieface) {
          const int ihp0 = faceHex[iface][(ieface + 0) % 4];
          const int ihp1 = faceHex[iface][(ieface + 1) % 4];
          const int ihpm = edgePointFace[iface][ieface];
          int ip_e0 = c.aIP[ihp0];
          int ip_em = c.aIP[ihpm];
          int ip_e1 = c.aIP[ihp1];
          if (ip_em < 0) {
            aTet.push_back(ip_c);
            aTet.push_back(ip_d);
            aTet.push_back(ip_e0);
            aTet.push_back(ip_e1); // bisect
          } else {
            if (c.aIC_Cld[ihp0] < 0) {
              aTet.push_back(ip_c);
              aTet.push_back(ip_d);
              aTet.push_back(ip_e0);
              aTet.push_back(ip_em); // quadrisect
            }
            if (c.aIC_Cld[ihp1] < 0) {
              aTet.push_back(ip_c);
              aTet.push_back(ip_d);
              aTet.push_back(ip_em);
              aTet.push_back(ip_e1); // quadrisect
            }
          }
        }
      }
    }
  }
}

DFM2_INLINE void WarpLattice(
    std::vector<CPointLattice> &aPointLattice,
    const std::vector<int> &psup_ind,
    const std::vector<int> &psup) {
  for (unsigned int ip = 0; ip < aPointLattice.size(); ++ip) {
    if (aPointLattice[ip].sdf < 0) { aPointLattice[ip].iflg = 1; }  // outside
    else { aPointLattice[ip].iflg = 2; } // inside
  }
  // warp background mesh
  for (int ip = 0; ip < (int) aPointLattice.size(); ++ip) { // lattice points
    double min_len = -1.0;
    double min_dist = -1.0;
    double min_cut[3] = {0.0, 0.0, 0.0};
    for (int ipsup = psup_ind[ip]; ipsup < psup_ind[ip + 1]; ++ipsup) {
      const int jp = psup[ipsup];
      assert(jp >= 0 && jp < (int) aPointLattice.size());
      if (aPointLattice[ip].iflg == aPointLattice[jp].iflg) continue;
      if (aPointLattice[ip].iflg == 0) continue;
      if (aPointLattice[jp].iflg == 0) continue;
      const double sdfi = aPointLattice[ip].sdf;
      const double sdfj = aPointLattice[jp].sdf;
      double ratioi = sdfj / (sdfj - sdfi);
      double ratioj = 1.0 - ratioi;
      double vi2j[3] = {
          (aPointLattice[ip].pos[0] - aPointLattice[jp].pos[0]),
          (aPointLattice[ip].pos[1] - aPointLattice[jp].pos[1]),
          (aPointLattice[ip].pos[2] - aPointLattice[jp].pos[2])};

      double len = sqrt(vi2j[0] * vi2j[0] + vi2j[1] * vi2j[1] + vi2j[2] * vi2j[2]);
      double dist = len * (1 - ratioi);
      if (min_dist < 0 || dist < min_dist) {
        min_dist = dist;
        min_cut[0] = aPointLattice[ip].pos[0] * ratioi + aPointLattice[jp].pos[0] * ratioj;
        min_cut[1] = aPointLattice[ip].pos[1] * ratioi + aPointLattice[jp].pos[1] * ratioj;
        min_cut[2] = aPointLattice[ip].pos[2] * ratioi + aPointLattice[jp].pos[2] * ratioj;
      }
      if (min_len < 0 || len < min_len) {
        min_len = len;
      }
    }
    if (min_dist < 0) continue;
    if (min_len < 0) continue;
    if (min_dist > min_len * 0.3) continue;
    aPointLattice[ip].pos[0] = min_cut[0];
    aPointLattice[ip].pos[1] = min_cut[1];
    aPointLattice[ip].pos[2] = min_cut[2];
    aPointLattice[ip].iflg = 0;
  }
}

DFM2_INLINE void MakeCutPoint(
    std::vector<double> &aXYZ,
    std::vector<int> &mapLat2Out,
    std::vector<int> &lat2cut_ind,
    std::vector<int> &lat2cut,
    const std::vector<CPointLattice> &aPointLattice,
    const std::vector<int> &psup_ind,
    const std::vector<int> &psup) {
  const int nno_lat = (int) aPointLattice.size();
  aXYZ.clear();
  aXYZ.reserve(nno_lat * 6);
  //
  mapLat2Out.assign(nno_lat, -1);
  int nno_out_lat = 0;
  for (int ino = 0; ino < nno_lat; ++ino) {
    if (aPointLattice[ino].iflg == 1) continue;
    mapLat2Out[ino] = nno_out_lat;
    aXYZ.push_back(aPointLattice[ino].pos[0]);
    aXYZ.push_back(aPointLattice[ino].pos[1]);
    aXYZ.push_back(aPointLattice[ino].pos[2]);
    nno_out_lat++;
  }
  ////
  lat2cut_ind.assign(nno_lat + 1, 0);
  for (int ino0 = 0; ino0 < nno_lat; ++ino0) {
    for (int ipsup = psup_ind[ino0]; ipsup < psup_ind[ino0 + 1]; ++ipsup) {
      const int ino1 = psup[ipsup];
      if (aPointLattice[ino0].iflg == aPointLattice[ino1].iflg) continue;
      if (aPointLattice[ino0].iflg == 0) continue;
      if (aPointLattice[ino1].iflg == 0) continue;
      lat2cut_ind[ino0 + 1]++; // register cut point
    }
  }
  //
  for (int ino = 0; ino < nno_lat; ++ino) { lat2cut_ind[ino + 1] += lat2cut_ind[ino]; }
  int icut_cur = nno_out_lat;
  int nno_cut = lat2cut_ind[nno_lat] / 2;
  aXYZ.resize(aXYZ.size() + nno_cut * 3);
  lat2cut.resize(nno_cut * 4, -1);
  for (int ino0 = 0; ino0 < nno_lat; ++ino0) {
    for (int ipsup = psup_ind[ino0]; ipsup < psup_ind[ino0 + 1]; ++ipsup) {
      const int ino1 = psup[ipsup];
      if (aPointLattice[ino0].iflg == aPointLattice[ino1].iflg) continue;
      if (aPointLattice[ino0].iflg == 0) continue;
      if (aPointLattice[ino1].iflg == 0) continue;
      const int ilat2cut0 = lat2cut_ind[ino0];
      lat2cut[ilat2cut0 * 2 + 0] = ino1; // oposite side
      if (ino0 > ino1) {
        lat2cut_ind[ino0]++;
        continue;
      }
      const double sdf0 = aPointLattice[ino0].sdf;
      const double sdf1 = aPointLattice[ino1].sdf;
      const double ratio0 = sdf1 / (sdf1 - sdf0);
      aXYZ[icut_cur * 3 + 0] = aPointLattice[ino0].pos[0] * ratio0 + aPointLattice[ino1].pos[0] * (1 - ratio0);
      aXYZ[icut_cur * 3 + 1] = aPointLattice[ino0].pos[1] * ratio0 + aPointLattice[ino1].pos[1] * (1 - ratio0);
      aXYZ[icut_cur * 3 + 2] = aPointLattice[ino0].pos[2] * ratio0 + aPointLattice[ino1].pos[2] * (1 - ratio0);
      lat2cut[ilat2cut0 * 2 + 1] = icut_cur; // cut
      icut_cur++;
      lat2cut_ind[ino0]++;
    }
  }
  assert(icut_cur == nno_cut + nno_out_lat);
  for (int ino = nno_lat; ino >= 1; ino--) { lat2cut_ind[ino] = lat2cut_ind[ino - 1]; }
  lat2cut_ind[0] = 0;
  ////
  for (int ino0 = 0; ino0 < nno_lat; ++ino0) {
    for (int ind0 = lat2cut_ind[ino0]; ind0 < lat2cut_ind[ino0 + 1]; ++ind0) {
      int ino1 = lat2cut[ind0 * 2 + 0];
      if (ino0 < ino1) {
        assert(lat2cut[ind0 * 2 + 1] != -1);
        continue;
      }
      for (int ind1 = lat2cut_ind[ino1]; ind1 < lat2cut_ind[ino1 + 1]; ind1++) {
        int ino2 = lat2cut[ind1 * 2 + 0];
        if (ino2 != ino0) { continue; }
        assert(lat2cut[ind1 * 2 + 1] != -1);
        lat2cut[ind0 * 2 + 1] = lat2cut[ind1 * 2 + 1];
        assert(lat2cut[ind0 * 2 + 1] != -1);
        break;
      }
    }
  }
}

DFM2_INLINE void cutoutTetFromLattice(
    std::vector<unsigned int> &aTet,
    const std::vector<CPointLattice> &aPointLattice,
    const std::vector<unsigned int> &aTetLattice,
    const std::vector<int> &mapLat2Out,
    const std::vector<int> &lat2cut_ind,
    const std::vector<int> &lat2cut) {
  aTet.clear();
  aTet.reserve(aTetLattice.size());
  for (unsigned int it = 0; it < aTetLattice.size() / 4; it++) {
    const int iln0 = aTetLattice[it * 4 + 0];
    const int iln1 = aTetLattice[it * 4 + 1];
    const int iln2 = aTetLattice[it * 4 + 2];
    const int iln3 = aTetLattice[it * 4 + 3];
    const int f0 = aPointLattice[iln0].iflg;
    const int f1 = aPointLattice[iln1].iflg;
    const int f2 = aPointLattice[iln2].iflg;
    const int f3 = aPointLattice[iln3].iflg;
    int on[10] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    FindCutNodeTet(on,
                   iln0, iln1, iln2, iln3, f0, f1, f2, f3,
                   mapLat2Out, lat2cut_ind, lat2cut);
    const int iflg = aPointLattice[iln0].iflg
        + aPointLattice[iln1].iflg * 3
        + aPointLattice[iln2].iflg * 9
        + aPointLattice[iln3].iflg * 27;
    int tet[3][4];
    unsigned int ntet;
    GetClampTet(tet, ntet, iflg, on);
    for (unsigned int itet = 0; itet < ntet; itet++) {
      aTet.push_back(tet[itet][0]);
      aTet.push_back(tet[itet][1]);
      aTet.push_back(tet[itet][2]);
      aTet.push_back(tet[itet][3]);
    }
  }
}

} // namespace iss
} // namespace delfem2

/**
 * @brief internal function for debug
 */
DFM2_INLINE void delfem2::makeBackgroundLattice
    (std::vector<CPointLattice> &aPointLattice,
     std::vector<unsigned int> &aTetLattice,
     const CInput_IsosurfaceStuffing &input,
     double elen, int ndiv, const double org[3]) {
  std::vector<iss::CCell> aCell;
  iss::makeLatticeCoasestLevel(aPointLattice, aCell,
                               input, elen, ndiv, org);
  for (int icell = 0; icell < (int) aCell.size(); ++icell) {
    int icntr0 = aCell[icell].aIP[26];
    double size_parent = aCell[icell].size;
    double sdf_parent = aPointLattice[icntr0].sdf;
    if (-sdf_parent
        > size_parent * 1.5) { continue; } // this cell is completely outside the region no need to look into children
    double cx = aPointLattice[icntr0].pos[0];
    double cy = aPointLattice[icntr0].pos[1];
    double cz = aPointLattice[icntr0].pos[2];
    const int ilevel_parent = aCell[icell].ilevel;
    for (int ichild = 0; ichild < 8; ++ichild) {
      const double ccx = cx + iss::aCellPointDirection[ichild][0] * size_parent * 0.25;
      const double ccy = cy + iss::aCellPointDirection[ichild][1] * size_parent * 0.25;
      const double ccz = cz + iss::aCellPointDirection[ichild][2] * size_parent * 0.25;
      double sdf0;
      int level_vol_goal, level_srf, nlayer;
      input.Level(level_vol_goal, level_srf, nlayer, sdf0,
                  ccx, ccy, ccz);
      int level_srf_goal = 0;
      if (level_srf > -0.1) {
        double elen_srf = elen / pow(2, level_srf);
        double dist = sdf0 - elen_srf * nlayer;
        if (dist < elen_srf) {
          level_srf_goal = level_srf;
        } else {
          level_srf_goal = static_cast<int>(level_srf - log2(dist / elen_srf));
        }
      }
//      =
      if (level_vol_goal >= (ilevel_parent + 1) ||
//         (level_srf_goal >= (ilevel_parent+1) && level_srf > level_srf_goal) )
          level_srf_goal >= (ilevel_parent + 1)) {
        aCell[icell].aIC_Cld[ichild] = (int) aCell.size();
        iss::CCell cc(size_parent * 0.5, ilevel_parent + 1, icell, ichild);
        {
          cc.aIP[26] = (int) aPointLattice.size();
          CPointLattice p(ccx, ccy, ccz, sdf0);
          aPointLattice.push_back(p);
        }
        aCell.push_back(cc);
      } else {
        aCell[icell].aIC_Cld[ichild] = -1;
      }
    }
    aCell[icell].setChildAdjRelation(aCell); // make relation ship inside children
  }
  Continuation(aPointLattice, aCell,
               input);
//  CheckContinuation(aCell);
  addEdgeFacePoints(aPointLattice, aCell,
                    input);
  makeTetLattice(aTetLattice,
                 aCell);
}

DFM2_INLINE bool delfem2::IsoSurfaceStuffing(
    std::vector<double> &aXYZ,
    std::vector<unsigned int> &aTet,
    std::vector<int> &aIsOnSurfXYZ,
    const CInput_IsosurfaceStuffing &input,
    double elen_in,
    double width,
    const double center[3]) {
  if (elen_in <= 0) return false;

  int ndiv = (int) (width / elen_in);
  if (ndiv == 0) { ndiv = 1; }
  double elen = width / (double) ndiv;
  const double org[3] = {center[0] - 0.5 * elen * ndiv, center[1] - 0.5 * elen * ndiv, center[2] - 0.5 * elen * ndiv};

  std::vector<CPointLattice> aPointLattice;
  std::vector<unsigned int> aTetLattice;
  makeBackgroundLattice(aPointLattice, aTetLattice, input, elen, ndiv, org);

  std::vector<int> mapLat2Out;
  std::vector<int> lat2cut_ind, lat2cut;
  {
    std::vector<int> psup_ind, psup;
    {
      std::vector<int> elsup_ind, elsup;
      iss::makeElemSurroundingPoint(elsup_ind, elsup, aTetLattice, 4, (int) aPointLattice.size());
      iss::makeOneRingNeighborhood(psup_ind, psup, aTetLattice, elsup_ind, elsup, 4, (int) aPointLattice.size());
    }
    iss::WarpLattice(aPointLattice,
                     psup_ind, psup);
    iss::MakeCutPoint(aXYZ, mapLat2Out, lat2cut_ind, lat2cut,
                      aPointLattice, psup_ind, psup);
  }

  iss::cutoutTetFromLattice(aTet,
                            aPointLattice, aTetLattice, mapLat2Out, lat2cut_ind, lat2cut);

  {
    aIsOnSurfXYZ.assign((int) aXYZ.size() / 3, 1);
    const int np = (int) aPointLattice.size();
    for (int ip = 0; ip < np; ++ip) {
      if (aPointLattice[ip].iflg == 0) {
        int ion0 = mapLat2Out[ip];
        if (ion0 == -1) continue;
        aIsOnSurfXYZ[ion0] = 1;
      }
      if (aPointLattice[ip].iflg == 2) {
        int ion0 = mapLat2Out[ip];
        if (ion0 == -1) continue;
        aIsOnSurfXYZ[ion0] = 0;
      }
    }
  }

  return true;
}

// ------------------------------------------------
// Legacy code for single level isosurface stuffing

#if 0

                                                                                                                        class CInputIsosurfaceStuffing
{
public:
  virtual double Projection(double px, double py, double pz) const = 0;
};

void GetCoord_BBCLatticePoint
(double co[3], unsigned int ipo,
 unsigned int ndiv, double elen, const double org[3])
{
  if( ipo < (ndiv+1)*(ndiv+1)*(ndiv+1) ){
    const int ipoz = ipo/((ndiv+1)*(ndiv+1));
    int ipoxy = ipo-ipoz*(ndiv+1)*(ndiv+1);
    const int ipoy = ipoxy/(ndiv+1);
    const int ipox = ipoxy-ipoy*(ndiv+1);
    co[0] = org[0]+ipox*elen;
    co[1] = org[1]+ipoy*elen;
    co[2] = org[2]+ipoz*elen;
  }
  else{
    int ipo0 = ipo-(ndiv+1)*(ndiv+1)*(ndiv+1);
    const int ipoz = ipo0/(ndiv*ndiv);
    unsigned int ipoxy = ipo0-ipoz*ndiv*ndiv;
    const int ipoy = ipoxy/ndiv;
    const int ipox = ipoxy-ipoy*ndiv;
    co[0] = org[0]+(ipox+0.5)*elen;
    co[1] = org[1]+(ipoy+0.5)*elen;
    co[2] = org[2]+(ipoz+0.5)*elen;
  }
}

void EvaluateSignedDistanceField_BBCLatticePoint
(std::vector<double>& sdf_lat, std::vector<double>& co_lat,
 const CInputIsosurfaceStuffing& sdf,
 unsigned int ndiv, double elen, const double org[3])
{
  const int nno_lat = (ndiv+1)*(ndiv+1)*(ndiv+1)+ndiv*ndiv*ndiv;
  sdf_lat.assign(nno_lat,-1);
  co_lat.clear();
  co_lat.resize(nno_lat*3);
  for(int ino=0;ino<nno_lat;++ino){
    double co[3]; GetCoord_BBCLatticePoint(co,ino, ndiv,elen,org);
    co_lat[ino*3+0] = co[0];
    co_lat[ino*3+1] = co[1];
    co_lat[ino*3+2] = co[2];
    double dist = sdf.Projection(co[0],co[1],co[2]);
    sdf_lat[ino] = dist;
  }
}


void GetAdjacentNode
(int adj[14],
 unsigned int ix, unsigned int iy, unsigned int iz, bool is_corner, unsigned int ndiv)
{
  const int nno_cor = (ndiv+1)*(ndiv+1)*(ndiv+1);
  const int nno_lat = (ndiv+1)*(ndiv+1)*(ndiv+1)+ndiv*ndiv*ndiv;
  adj[ 0]=-1; adj[ 1]=-1; adj[ 2]=-1; adj[ 3]=-1; adj[ 4]=-1; adj[ 5]=-1;
  adj[ 6]=-1; adj[ 7]=-1; adj[ 8]=-1; adj[ 9]=-1; adj[10]=-1; adj[11]=-1; adj[12]=-1; adj[13]=-1;
  if( is_corner ){
    const unsigned int ino0 = iz*(ndiv+1)*(ndiv+1)+iy*(ndiv+1)+ix;
    if( ix!=0    ){ adj[0] = ino0-1;                 assert( adj[0] >= 0       ); }
    if( ix!=ndiv ){ adj[1] = ino0+1;                 assert( adj[1] <  nno_cor ); }
    if( iy!=0    ){ adj[2] = ino0-(ndiv+1);          assert( adj[2] >= 0       ); }
    if( iy!=ndiv ){ adj[3] = ino0+(ndiv+1);          assert( adj[3] <  nno_cor ); }
    if( iz!=0    ){ adj[4] = ino0-(ndiv+1)*(ndiv+1); assert( adj[4] >= 0       ); }
    if( iz!=ndiv ){ adj[5] = ino0+(ndiv+1)*(ndiv+1); assert( adj[5] <  nno_cor ); }
    ////
    if( ix!=0    && iy!=0    && iz!=0    ){ adj[ 6] = (ix-1)+(iy-1)*ndiv+(iz-1)*ndiv*ndiv + nno_cor; }
    if( ix!=ndiv && iy!=0    && iz!=0    ){ adj[ 7] = (ix-0)+(iy-1)*ndiv+(iz-1)*ndiv*ndiv + nno_cor; }
    if( ix!=0    && iy!=ndiv && iz!=0    ){ adj[ 8] = (ix-1)+(iy-0)*ndiv+(iz-1)*ndiv*ndiv + nno_cor; }
    if( ix!=ndiv && iy!=ndiv && iz!=0    ){ adj[ 9] = (ix-0)+(iy-0)*ndiv+(iz-1)*ndiv*ndiv + nno_cor; }
    if( ix!=0    && iy!=0    && iz!=ndiv ){ adj[10] = (ix-1)+(iy-1)*ndiv+(iz-0)*ndiv*ndiv + nno_cor; }
    if( ix!=ndiv && iy!=0    && iz!=ndiv ){ adj[11] = (ix-0)+(iy-1)*ndiv+(iz-0)*ndiv*ndiv + nno_cor; }
    if( ix!=0    && iy!=ndiv && iz!=ndiv ){ adj[12] = (ix-1)+(iy-0)*ndiv+(iz-0)*ndiv*ndiv + nno_cor; }
    if( ix!=ndiv && iy!=ndiv && iz!=ndiv ){ adj[13] = (ix-0)+(iy-0)*ndiv+(iz-0)*ndiv*ndiv + nno_cor; }
  }
  else{
    unsigned int ino0 = iz*ndiv*ndiv+iy*ndiv+ix + nno_cor;
    if( ix!=0      ){ adj[0] = ino0-1;         assert( adj[0] >= nno_cor ); }
    if( ix!=ndiv-1 ){ adj[1] = ino0+1;         assert( adj[1] <  nno_lat ); }
    if( iy!=0      ){ adj[2] = ino0-ndiv;      assert( adj[2] >= nno_cor ); }
    if( iy!=ndiv-1 ){ adj[3] = ino0+ndiv;      assert( adj[3] <  nno_lat ); }
    if( iz!=0      ){ adj[4] = ino0-ndiv*ndiv; assert( adj[4] >= nno_cor ); }
    if( iz!=ndiv-1 ){ adj[5] = ino0+ndiv*ndiv; assert( adj[5] <  nno_lat ); }
    ////
    adj[ 6] = (ix+0)+(iy+0)*(ndiv+1)+(iz+0)*(ndiv+1)*(ndiv+1);
    adj[ 7] = (ix+1)+(iy+0)*(ndiv+1)+(iz+0)*(ndiv+1)*(ndiv+1);
    adj[ 8] = (ix+0)+(iy+1)*(ndiv+1)+(iz+0)*(ndiv+1)*(ndiv+1);
    adj[ 9] = (ix+1)+(iy+1)*(ndiv+1)+(iz+0)*(ndiv+1)*(ndiv+1);
    adj[10] = (ix+0)+(iy+0)*(ndiv+1)+(iz+1)*(ndiv+1)*(ndiv+1);
    adj[11] = (ix+1)+(iy+0)*(ndiv+1)+(iz+1)*(ndiv+1)*(ndiv+1);
    adj[12] = (ix+0)+(iy+1)*(ndiv+1)+(iz+1)*(ndiv+1)*(ndiv+1);
    adj[13] = (ix+1)+(iy+1)*(ndiv+1)+(iz+1)*(ndiv+1)*(ndiv+1);
  }
}

void GetAdjacentNode
(int adj[14],
 const int ino0, int ndiv)
{
  assert( ino0 >= 0 );
  const int nno_cor = (ndiv+1)*(ndiv+1)*(ndiv+1);
  if( ino0 < nno_cor ){
    const int iz = ino0/((ndiv+1)*(ndiv+1));
    const int ino1 = ino0-iz*(ndiv+1)*(ndiv+1);
    const int iy = ino1/(ndiv+1);
    const int ix = ino1-iy*(ndiv+1);
    GetAdjacentNode(adj, ix, iy, iz, true, ndiv);
  }
  else{
    assert( ino0 < (ndiv+1)*(ndiv+1)*(ndiv+1)+ndiv*ndiv*ndiv );
    const int ino1 = ino0-nno_cor;
    const int iz = ino1/(ndiv*ndiv);
    const int ino2 = ino1-ndiv*ndiv*iz;
    const int iy = ino2/ndiv;
    const int ix = ino2-iy*ndiv;
    GetAdjacentNode(adj, ix, iy, iz, false, ndiv);
  }
}


void MakeCoordCutPoints
(std::vector<double>& aXYZ, std::vector<int>& mapLat2Out,
 std::vector<int>& lat2cut_ind, std::vector<int>& lat2cut,
 const std::vector<double>& aSDF_lat, const std::vector<double>& aXYZ_lat, const std::vector<unsigned int>& aflg_lat,
 unsigned int ndiv)
{
  const int nno_lat = (ndiv+1)*(ndiv+1)*(ndiv+1)+ndiv*ndiv*ndiv; // lattice point
  assert( aSDF_lat.size()  == nno_lat );
  assert( aXYZ_lat.size()   == nno_lat*3 );
  assert( aflg_lat.size() == nno_lat );
  aXYZ.clear();
  aXYZ.reserve(nno_lat*6);
  /////
  mapLat2Out.assign(nno_lat,-1);
  int nno_out_lat = 0;
  for(int ino=0;ino<nno_lat;++ino){
    if( aflg_lat[ino] == 1 ) continue; // this is outside of mesh
    mapLat2Out[ino] = nno_out_lat;
    aXYZ.push_back( aXYZ_lat[ino*3+0] );
    aXYZ.push_back( aXYZ_lat[ino*3+1] );
    aXYZ.push_back( aXYZ_lat[ino*3+2] );
    nno_out_lat++;
  }
  ////
  lat2cut_ind.assign(nno_lat+1,0);
  for(int ino0=0;ino0<nno_lat;++ino0){
    int adj[14]; GetAdjacentNode(adj, ino0,ndiv);
    for(int iadj=0;iadj<14;++iadj){
      if( adj[iadj] == -1 ) continue;
      const int ino1 = adj[iadj];
      if( aflg_lat[ino0] == aflg_lat[ino1] ) continue;
      if( aflg_lat[ino0] == 0 ) continue;
      if( aflg_lat[ino1] == 0 ) continue;
      lat2cut_ind[ino0+1]++; // register cut point
    }
  }
  /////
  for(int ino=0;ino<nno_lat;++ino){ lat2cut_ind[ino+1] += lat2cut_ind[ino]; }
  int icut_cur = nno_out_lat;
  int nno_cut = lat2cut_ind[nno_lat]/2;
  aXYZ.resize(aXYZ.size()+nno_cut*3);
  lat2cut.resize(nno_cut*4,-1);
  for(int ino0=0;ino0<nno_lat;++ino0){
    int adj[14];   GetAdjacentNode(adj, ino0,ndiv);
    for(int iadj=0;iadj<14;++iadj){
      if( adj[iadj] == -1 ) continue;
      const int ino1 = adj[iadj];
      if( aflg_lat[ino0] == aflg_lat[ino1] ) continue;
      if( aflg_lat[ino0] == 0 ) continue;
      if( aflg_lat[ino1] == 0 ) continue;
      const int ilat2cut0 = lat2cut_ind[ino0];
      lat2cut[ilat2cut0*2+0] = ino1; // oposite side
      if( ino0 > ino1 ){
        lat2cut_ind[ino0]++;
        continue;
      }
      const double sdf0 = aSDF_lat[ino0];
      const double sdf1 = aSDF_lat[ino1];
      const double ratio0 = sdf1/(sdf1-sdf0);
      aXYZ[icut_cur*3+0] = aXYZ_lat[ino0*3+0]*ratio0 + aXYZ_lat[ino1*3+0]*(1-ratio0);
      aXYZ[icut_cur*3+1] = aXYZ_lat[ino0*3+1]*ratio0 + aXYZ_lat[ino1*3+1]*(1-ratio0);
      aXYZ[icut_cur*3+2] = aXYZ_lat[ino0*3+2]*ratio0 + aXYZ_lat[ino1*3+2]*(1-ratio0);
      lat2cut[ilat2cut0*2+1] = icut_cur; // cut
      icut_cur++;
      lat2cut_ind[ino0]++;
    }
  }
  assert( icut_cur == nno_cut+nno_out_lat );
  for(int ino=nno_lat;ino>=1;ino--){ lat2cut_ind[ino] = lat2cut_ind[ino-1]; }
  lat2cut_ind[0] = 0;
  ////
  for(int ino0=0;ino0<nno_lat;++ino0){
    for(int ind0=lat2cut_ind[ino0];ind0<lat2cut_ind[ino0+1];++ind0){
      int ino1 = lat2cut[ind0*2+0];
      if( ino0 < ino1 ){
        int iout = lat2cut[ind0*2+1];
        assert( iout != -1 );
        continue;
      }
      for(int ind1=lat2cut_ind[ino1];ind1<lat2cut_ind[ino1+1];ind1++){
        int ino2 = lat2cut[ind1*2+0];
        if( ino2 != ino0 ){ continue; }
        assert( lat2cut[ind1*2+1] != -1 );
        lat2cut[ind0*2+1] = lat2cut[ind1*2+1];
        assert( lat2cut[ind0*2+1] != -1 );
        break;
      }
    }
  }
}

// 12 tetrahedrons around the body point
// call by "MarchingTet"
void GetBBCTet
(int ix, int iy, int iz,
 int lnods[12][4],
 int ndiv)
{
  for(unsigned int i=0;i<12*4;i++){ (&lnods[0][0])[i] = -1; }
  const int ib0 = ix+iy*ndiv+iz*ndiv*ndiv + (ndiv+1)*(ndiv+1)*(ndiv+1); // body point
  const int ic000 = (ix+0)+(iy+0)*(ndiv+1)+(iz+0)*(ndiv+1)*(ndiv+1);
  const int ic100 = (ix+1)+(iy+0)*(ndiv+1)+(iz+0)*(ndiv+1)*(ndiv+1);
  const int ic010 = (ix+0)+(iy+1)*(ndiv+1)+(iz+0)*(ndiv+1)*(ndiv+1);
  const int ic110 = (ix+1)+(iy+1)*(ndiv+1)+(iz+0)*(ndiv+1)*(ndiv+1);
  const int ic001 = (ix+0)+(iy+0)*(ndiv+1)+(iz+1)*(ndiv+1)*(ndiv+1);
  const int ic101 = (ix+1)+(iy+0)*(ndiv+1)+(iz+1)*(ndiv+1)*(ndiv+1);
  const int ic011 = (ix+0)+(iy+1)*(ndiv+1)+(iz+1)*(ndiv+1)*(ndiv+1);
  const int ic111 = (ix+1)+(iy+1)*(ndiv+1)+(iz+1)*(ndiv+1)*(ndiv+1);
  if( ix != ndiv-1 ){ // +x
    lnods[ 0][0]=ib0;  lnods[ 0][2]=ic101;  lnods[ 0][3]=ic100;  lnods[ 0][1]=ib0+1;
    lnods[ 1][0]=ib0;  lnods[ 1][2]=ic100;  lnods[ 1][3]=ic110;  lnods[ 1][1]=ib0+1;
    lnods[ 2][0]=ib0;  lnods[ 2][2]=ic110;  lnods[ 2][3]=ic111;  lnods[ 2][1]=ib0+1;
    lnods[ 3][0]=ib0;  lnods[ 3][2]=ic111;  lnods[ 3][3]=ic101;  lnods[ 3][1]=ib0+1;
  }
  if( iy != ndiv-1 ){ // +y
    lnods[ 4][0]=ib0;  lnods[ 4][2]=ic111;  lnods[ 4][3]=ic110;  lnods[ 4][1]=ib0+ndiv;
    lnods[ 5][0]=ib0;  lnods[ 5][2]=ic110;  lnods[ 5][3]=ic010;  lnods[ 5][1]=ib0+ndiv;
    lnods[ 6][0]=ib0;  lnods[ 6][2]=ic010;  lnods[ 6][3]=ic011;  lnods[ 6][1]=ib0+ndiv;
    lnods[ 7][0]=ib0;  lnods[ 7][2]=ic011;  lnods[ 7][3]=ic111;  lnods[ 7][1]=ib0+ndiv;
  }
  if( iz != ndiv-1 ){ // +z
    lnods[ 8][0]=ib0;  lnods[ 8][2]=ic001;  lnods[ 8][3]=ic101;  lnods[ 8][1]=ib0+ndiv*ndiv;
    lnods[ 9][0]=ib0;  lnods[ 9][2]=ic101;  lnods[ 9][3]=ic111;  lnods[ 9][1]=ib0+ndiv*ndiv;
    lnods[10][0]=ib0;  lnods[10][2]=ic111;  lnods[10][3]=ic011;  lnods[10][1]=ib0+ndiv*ndiv;
    lnods[11][0]=ib0;  lnods[11][2]=ic011;  lnods[11][3]=ic001;  lnods[11][1]=ib0+ndiv*ndiv;
  }
}


void Warp_BBCLatticePoint
(std::vector<double>& aXYZ_lat,
 std::vector<unsigned int>& flg_lat,
 const std::vector<double>& sdf_lat,
 unsigned int ndiv, double elen)
{
  const int nno_lat = (ndiv+1)*(ndiv+1)*(ndiv+1)+ndiv*ndiv*ndiv; // lattice point
  assert( sdf_lat.size()  == nno_lat );
  assert( aXYZ_lat.size()   == nno_lat*3 );
  /////
  flg_lat.resize(nno_lat);
  for(int ino=0;ino<nno_lat;++ino){
    if( sdf_lat[ino] < 0 ){ flg_lat[ino] = 1; } // outside
    else{                   flg_lat[ino] = 2; } // inside
  }
  for(int ino0=0;ino0<nno_lat;++ino0){ // lattice points
    int adj[14]; GetAdjacentNode(adj, ino0,ndiv);
    double min_dist = -1.0;
    double min_cut[3];
    for(int iadj=0;iadj<14;iadj++){ // find out minium distance to the surface
      if( adj[iadj] == -1 ) continue;
      const int ino1 = adj[iadj];
      if( flg_lat[ino0] == flg_lat[ino1] ) continue;
      if( flg_lat[ino0] == 0 ) continue;
      if( flg_lat[ino1] == 0 ) continue;
      const double sdf0 = sdf_lat[ino0];
      const double sdf1 = sdf_lat[ino1];
      double ratio0 = sdf1/(sdf1-sdf0);
      double disp[3] = {
        (aXYZ_lat[ino0*3+0]-aXYZ_lat[ino1*3+0])*(1-ratio0),
        (aXYZ_lat[ino0*3+1]-aXYZ_lat[ino1*3+1])*(1-ratio0),
        (aXYZ_lat[ino0*3+2]-aXYZ_lat[ino1*3+2])*(1-ratio0) };
      double dist = SqLength3D(disp);
      if( min_dist > 0 && dist > min_dist  ) continue;
      min_dist = dist;
      min_cut[0] = aXYZ_lat[ino0*3+0]*ratio0+aXYZ_lat[ino1*3+0]*(1-ratio0);
      min_cut[1] = aXYZ_lat[ino0*3+1]*ratio0+aXYZ_lat[ino1*3+1]*(1-ratio0);
      min_cut[2] = aXYZ_lat[ino0*3+2]*ratio0+aXYZ_lat[ino1*3+2]*(1-ratio0);
    }
    if( min_dist < 0 ) continue;
    min_dist = sqrt(min_dist);
    if( min_dist > elen*0.3 ) continue;
    aXYZ_lat[ino0*3+0] = min_cut[0];
    aXYZ_lat[ino0*3+1] = min_cut[1];
    aXYZ_lat[ino0*3+2] = min_cut[2];
    flg_lat[ino0] = 0;
  }
}


void MarchingTet
(std::vector<int>& aTet,
 std::vector<int>& mapLat2Out,
 std::vector<int>& aCutInd, std::vector<int>& aCut,
 const std::vector<unsigned int>& flg_lat,
 unsigned int ndiv)
{
  const int nno_lat = (ndiv+1)*(ndiv+1)*(ndiv+1)+ndiv*ndiv*ndiv;
  assert( flg_lat.size() == nno_lat );
  aTet.clear();
  aTet.reserve(ndiv*ndiv*ndiv*24*4);
  for(int ix=0;ix<ndiv;ix++){
    for(int iy=0;iy<ndiv;iy++){
      for(int iz=0;iz<ndiv;iz++){
        int lnods[12][4]; GetBBCTet(ix, iy, iz, lnods, ndiv); // 12 tetrahedrons around the body point
        for(int i0=0;i0<12;++i0){
          if( lnods[i0][0] == -1 ) continue;
          const int iln0 = lnods[i0][0];
          const int iln1 = lnods[i0][1];
          const int iln2 = lnods[i0][2];
          const int iln3 = lnods[i0][3];
          int on[10] = { -1,-1,-1,-1,-1, -1,-1,-1,-1,-1 };
          const int f0 = flg_lat[iln0];
          const int f1 = flg_lat[iln1];
          const int f2 = flg_lat[iln2];
          const int f3 = flg_lat[iln3];
          FindCutNodeTet(on,
                         iln0,iln1,iln2,iln3,
                         f0,f1,f2,f3,
                         mapLat2Out,aCutInd,aCut);
          unsigned int iflg = flg_lat[iln0] + flg_lat[iln1]*3 + flg_lat[iln2]*9 + flg_lat[iln3]*27;
          int tet[3][4];
          unsigned int ntet;
          GetClampTet(tet, ntet, iflg, on);
          for(int itet=0;itet<ntet;itet++){
            aTet.push_back(tet[itet][0]);
            aTet.push_back(tet[itet][1]);
            aTet.push_back(tet[itet][2]);
            aTet.push_back(tet[itet][3]);
          }            
        }
      }    
    }
  }
}

bool IsoSurfaceStuffing
(std::vector<double>& aXYZ, std::vector<int>& aTet, std::vector<int>& aIsSurfNode,
 const CInputIsosurfaceStuffing& sdf, double elen_in, double width, const double center[3])
{
  if( elen_in <= 0 ) return false;
  ////
  aXYZ.clear();
  aTet.clear();
  aIsSurfNode.clear();
  ////
  unsigned int ndiv = (int)(width/elen_in);
  if( ndiv == 0 ){ ndiv = 1; }
  double elen = width/(double)ndiv;
  const double org[3] = { center[0]-0.5*width, center[1]-0.5*width, center[2]-0.5*width };
  ////
  std::vector<double> aSDF_lat;
  std::vector<double> aXYZ_lat;
  EvaluateSignedDistanceField_BBCLatticePoint(aSDF_lat,aXYZ_lat,
                                              sdf, ndiv,elen,org);
  ////
  std::vector<unsigned int> flg_lat;
  Warp_BBCLatticePoint(aXYZ_lat,flg_lat,
                       aSDF_lat,ndiv,elen);
  ////
  std::vector<int> mapLat2Out;
  std::vector<int> lat2cut_ind, lat2cut;
  MakeCoordCutPoints(aXYZ, mapLat2Out,
                     lat2cut_ind,lat2cut,
                     aSDF_lat,aXYZ_lat,flg_lat, ndiv);
  ////
  MarchingTet(aTet,
              mapLat2Out,lat2cut_ind,lat2cut,flg_lat,ndiv);
  ////
  {
    aIsSurfNode.assign(aXYZ.size()/3,1);
    const int nno_lat = (ndiv+1)*(ndiv+1)*(ndiv+1)+ndiv*ndiv*ndiv;
    for(int iln=0;iln<nno_lat;++iln){
      if( flg_lat[iln] == 0 ){
        int ion0 = mapLat2Out[iln];
        if( ion0 == -1 ) continue;
        aIsSurfNode[ion0] = 1;
      }
      if( flg_lat[iln] == 2 ){
        int ion0 = mapLat2Out[iln];
        if( ion0 == -1 ) continue;
        aIsSurfNode[ion0] = 0;
      }
    }
  }
  for(int itet=0;itet<aTet.size()/4;itet++){
    const int i0 = aTet[itet*4+0];
    const int i1 = aTet[itet*4+1];
    const int i2 = aTet[itet*4+2];
    const int i3 = aTet[itet*4+3];
    const double p0[3] = { aXYZ[i0*3+0], aXYZ[i0*3+1], aXYZ[i0*3+2] };
    const double p1[3] = { aXYZ[i1*3+0], aXYZ[i1*3+1], aXYZ[i1*3+2] };
    const double p2[3] = { aXYZ[i2*3+0], aXYZ[i2*3+1], aXYZ[i2*3+2] };
    const double p3[3] = { aXYZ[i3*3+0], aXYZ[i3*3+1], aXYZ[i3*3+2] };
    double vol = TetVolume3D(p0, p1, p2, p3);
    if( vol < 0 ){
      std::cout << "nega : " << itet << " " << vol << std::endl;
    }
  }
  
  /*
   /////
   nXYZ_ = (int)aXYZ.size()/3;
   if( pXYZ_ != 0 ){ delete[] pXYZ_; }
   pXYZ_ = new double [aXYZ.size()];
   for(unsigned int i=0;i<aXYZ.size();i++){ pXYZ_[i] = aXYZ[i]; }
   //
   nTet_ = (int)aTet.size()/4;
   if( pTet_ != 0 ){ delete[] pTet_; }
   pTet_ = new unsigned int [aTet.size()];
   for(unsigned int i=0;i<aTet.size();i++){ pTet_[i] = aTet[i]; }
   for(unsigned int itet=0;itet<nTet_;itet++){
   unsigned int* tet = pTet_+itet*4;
   double vol = TetVolume(pXYZ_+tet[0]*3, pXYZ_+tet[1]*3, pXYZ_+tet[2]*3, pXYZ_+tet[3]*3);
   if( vol < 0 ){
   std::cout << "nega : " << itet << " " << vol << std::endl;
   }
   }
   */
  return true;
}


#endif







