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
#include "delfem2/mat3_funcs.h"

#ifndef M_PI
#  define M_PI 3.14159265359
#endif

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

class HairDarbouxShape {
 public:
  void MakeConfigDaboux(
    std::vector<CVec3d> & aP,
    std::vector<CVec3d> & aS) const {
    const unsigned int np = num_points;
    const double elen = edge_length;
    const double rad0 = radius;
    double dangle;
    {
      double tmp = pitch / (2 * M_PI);
      double tmp2 = tmp * tmp;
      dangle = elen / std::sqrt(rad0*rad0 + tmp2);
      for (unsigned int itr = 0; itr < 10; ++itr) {
        double elen0 = std::sqrt(2 * rad0 * rad0 * (1 - cos(dangle)) + tmp2 * dangle * dangle);
        double df = (rad0 * rad0 * sin(dangle) + dangle * tmp2) / elen0;
        double f = elen0 - elen;
        dangle -= f / df;
      }
    }
    aP.clear();
    for (unsigned int ip = 0; ip < np; ++ip) {
      const double angle = ip * dangle;
      aP.emplace_back(
        pitch * angle / (2*M_PI),
        rad0 * cos(angle),
        rad0 * sin(angle));
    }
    aS.resize(np);
    {
      const CVec3d v = (aP[1] - aP[0]).normalized();
      CVec3d s(0, 1, 0);
      s = (s - (s.dot(v)) * v).normalized();
      aS[0] = s;
    }
    for(unsigned int ip=1;ip<np-1;++ip) {
      const CMat3d R01 = Mat3_MinimumRotation(
        (aP[ip] - aP[ip-1]).normalized(),
        (aP[ip+1] - aP[ip]).normalized());
      aS[ip] = R01 * aS[ip-1];
    }
    aS[np-1] = CVec3d(1,0,0);
  }

 public:
  unsigned int num_points;

  //! axial position increment for the helix
  double pitch;

  double radius;

  double edge_length;
};

/*
void MakeProblemSetting_Spiral(
    std::vector<CVec3d> &aP0,
    std::vector<CVec3d> &aS0,
    std::vector<unsigned int> &aIP_HairRoot,
    const std::vector<HairDarbouxShape> &aHairShape);
*/

} // namespace delfem2

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/hair_darboux_util.cpp"
#endif

#endif  /* DFM2_HAIR_DARBOUX_H */
