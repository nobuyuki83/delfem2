/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_HAIR_DARBOUX_COLLISION_H
#define DFM2_HAIR_DARBOUX_COLLISION_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/hair_darboux_solver.h"
#include "delfem2/ls_solver_block_sparse.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/lsmats.h"

namespace delfem2 {

class CContactHair {
 public:
  unsigned int ip0, ip1;
  double s;
  unsigned int iq0, iq1;
  double t;
  CVec3d norm;
 public:
  [[nodiscard]] CVec3d Direction(
      const std::vector<CVec3d> &aP) const {
    return (1 - s) * aP[ip0] + s * aP[ip1] - (1 - t) * aP[iq0] - t * aP[iq1];
  }
};

DFM2_INLINE void Solve_RodHairContact(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    delfem2::LinearSystemSolver_BlockSparse &mats,
    double stiff_stretch,
    const double stiff_bendtwist[3],
    double mdtt,
    const std::vector<CVec3d> &aPt0,
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aS0,
    const std::vector<unsigned int> &aIP_HairRoot,
    double clearance,
    double stiff_contact,
    const std::vector<CContactHair> &aContact);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/hair_darboux_collision.cpp"
#endif

#endif  /* DFM2_HAIR_DARBOUX_COLLISION_H */
