//
// Created by Nobuyuki Umetani on 2022/01/01.
//

#ifndef POLYLINE_ELASTIC_EDIT_ROD2_H_
#define POLYLINE_ELASTIC_EDIT_ROD2_H_

#include <array>

#include "delfem2/ls_pentadiagonal.h"
#include "delfem2/fem_rod2.h"

template<class ARRAY2D>
void DragPolylineElastic_Rod2(
    ARRAY2D &vtx_xy,
    delfem2::LinearSystemSolver_BlockPentaDiagonal<2> &sparse,
    unsigned int idx_vtx0,
    const std::array<double,2> &pos_vtx0,
    double edge_length = 0.05,
    double stiff_bend = 0.001) {
  const double stiff_stretch = 1.0;
  const size_t np = vtx_xy.size();
  if( np < 3 ){
    return;
  }
  if( sparse.nblk() != np ){
    sparse.Initialize(np);
  }
  {
    sparse.dof_bcflag.assign(np * 2, 0);
    sparse.dof_bcflag[idx_vtx0 * 2 + 0] = 1;
    sparse.dof_bcflag[idx_vtx0 * 2 + 1] = 1;
  }
  sparse.BeginMerge();
  double W = 0.0;
  for (unsigned int ihinge = 0; ihinge < np - 2; ++ihinge) {
    const unsigned int aIP[3] = {ihinge, ihinge + 1, ihinge + 2};
    double ap[3][2] = {
        {vtx_xy[aIP[0]][0], vtx_xy[aIP[0]][1]},
        {vtx_xy[aIP[1]][0], vtx_xy[aIP[1]][1]},
        {vtx_xy[aIP[2]][0], vtx_xy[aIP[2]][1]}};
    for (int ino = 0; ino < 3; ++ino) {
      if (aIP[ino] != idx_vtx0) { continue; }
      ap[ino][0] = pos_vtx0[0];
      ap[ino][1] = pos_vtx0[1];
    }
    const double aL[2] = {edge_length, edge_length};
    double We, dWe[3][2], ddWe[3][3][2][2];
    delfem2::WdWddW_Rod2(
        We, dWe, ddWe,
        ap, aL, stiff_stretch, stiff_stretch, stiff_bend);
    W += We;
    for (int ino = 0; ino < 3; ++ino) {
      sparse.vec_r[aIP[ino] * 2 + 0] -= dWe[ino][0];
      sparse.vec_r[aIP[ino] * 2 + 1] -= dWe[ino][1];
    }
    sparse.template Merge<3, 3, 2, 2>(aIP, aIP, ddWe);
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    sparse.AddValueToDiagonal(ip, 0, 1.);
    sparse.AddValueToDiagonal(ip, 1, 1.);
  }
  sparse.Solve();
  // std::cout << W << std::endl;
  for (unsigned int ip = 0; ip < np; ++ip) {
    vtx_xy[ip][0] += static_cast<float>(sparse.vec_x[ip * 2 + 0]);
    vtx_xy[ip][1] += static_cast<float>(sparse.vec_x[ip * 2 + 1]);
  }
  { // fixed boundary condition
    vtx_xy[idx_vtx0][0] = pos_vtx0[0];
    vtx_xy[idx_vtx0][1] = pos_vtx0[1];
  }
}

#endif //POLYLINE_ELASTIC_EDIT_ROD2_H_
