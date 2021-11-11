/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_HAIR_DARBOUX_H
#define DFM2_HAIR_DARBOUX_H

#include "delfem2/dfm2_inline.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"

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

DFM2_INLINE void MakeDirectorOrthogonal_RodHair(
    std::vector<CVec3d> &aS,
    const std::vector<CVec3d> &aP);

DFM2_INLINE void UpdateSolutionHair(
    std::vector<CVec3d> &aP,
    std::vector<CVec3d> &aS,
    const std::vector<double> &vec_x,
    const std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<int> &aBCFlag);

class CHairShape {
 public:
  //! number of points
  unsigned int np;

  //! axis position increment for the helix
  double pitch;

  //! radius of helix
  double rad0;

  //! angle increment for the helix
  double dangle;

  //! root shape
  double p0[3];
};

void MakeProblemSetting_Spiral(
    std::vector<CVec3d> &aP0,
    std::vector<CVec3d> &aS0,
    std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<CHairShape> &aHairShape);

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/hair_darboux_util.cpp"
#endif

#endif  /* DFM2_HAIR_DARBOUX_H */
