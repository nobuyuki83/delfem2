/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_FEMROD_H
#define DFM2_FEMROD_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/lsmats.h"

namespace delfem2 {

// --------------
// below: rod hair

DFM2_INLINE void ParallelTransport_RodHair(
    std::vector<CVec3d> &aP0,
    std::vector<CVec3d> &aS0,
    const std::vector<unsigned int> &aIP_HairRoot);

DFM2_INLINE void MakeBCFlag_RodHair(
    std::vector<int> &aBCFlag,
    const std::vector<unsigned int> &aIP_HairRoot);

DFM2_INLINE void MakeSparseMatrix_RodHair(
    CMatrixSparse<double> &mats,
    const std::vector<unsigned int> &aIP_HairRoot,
    unsigned int ndof_par_node);

DFM2_INLINE void MakeDirectorOrthogonal_RodHair(
    std::vector<CVec3d> &aS,
    const std::vector<CVec3d> &aP);

/**
 * Merge linear system for hair
 * @param vec_r
 * @param[in,out] mats sparse matrix
 * @param[in] stiff_stretch
 * @param[in] stiff_bendtwist
 * @param[in] aIP_HairRoot
 * @param[in] aP
 * @param[in] aS
 * @param[in] aP0 vertex position of rest shape
 * @param[in] aS0
 * @return
 */
DFM2_INLINE double MergeLinSys_Hair(
    std::vector<double> &vec_r,
    CMatrixSparse<double> &mats,
    double stiff_stretch,
    const double stiff_bendtwist[3],
    const std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<CVec3d> &aP,
    const std::vector<CVec3d> &aS,
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aS0);

DFM2_INLINE void UpdateSolutionHair(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    const std::vector<double> &vec_x,
    const std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<int> &aBCFlag);

/**
 * @brief static minimization of the deformation energy
 * @param aP (in&out) position of the vertices of the rods
 * @param aS (in&out) director vectors
 * @param mats (in&out) sparse matrix
 * @param mdtt (in) mass divided by square of timestep (mass/dt/dt)
 * @param aP0 (in) initial positions of the vertices of the rods
 * @param aS0 (in) initial darboux vectors
 * @param aBCFlag (in) boundary condition flag. Non zero value means fixed value
 * @param aIP_HairRoot (in) indeces of the root points
 */
DFM2_INLINE void Solve_RodHair(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    CMatrixSparse<double> &mats,
    double stiff_stretch,
    const double stiff_bendtwist[3],
    double mdtt,
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aS0,
    const std::vector<int> &aBCFlag,
    const std::vector<unsigned int> &aIP_HairRoot);

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
    CMatrixSparse<double> &mats,
    double stiff_stretch,
    const double stiff_bendtwist[3],
    double mdtt,
    const std::vector<CVec3d> &aPt0,
    const std::vector<CVec3d> &aP0,
    const std::vector<CVec3d> &aS0,
    const std::vector<int> &aBCFlag,
    const std::vector<unsigned int> &aIP_HairRoot,
    double clearance,
    double stiff_contact,
    const std::vector<CContactHair> &aContact);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/femrod.cpp"
#endif

#endif  /* DFM2_FEMROD_H */
