/*
 * Copyright (c) 2021 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_RGD_V2M3_H
#define DFM2_RGD_V2M3_H

#include "delfem2/vec2.h"
#include "delfem2/mat3.h"

namespace delfem2 {


namespace rgd_v2m3{

CMat3d Mat3_Affine(
    const CVec2d& posl,
    const double theta,
    const CVec2d& posg);

}

class CRigidState2 {
public:
  bool is_fix;
  std::vector<CVec2d> shape;
  // below: derived from shape
  CVec2d posl;
  double mass;
  double I;
  //
  // below: change over time
  CVec2d posg;
  double theta;
  CVec2d velo;
  double omega;
  // below: temp data for PBD
  CVec2d posg_tmp;
  double theta_tmp;
};

class CContact {
public:
  unsigned int ir, jr;
  CVec2d Pi, Pjn, Njn;
  double lambda;
};

void Steptime_Rgd2(
    std::vector<CRigidState2>& aRS,
    double dt,
    const CVec2d& gravity);

}

#ifdef DFM2_HEADER_ONLY
  #include "delfem2/rgd_v2m3.cpp"
#endif

#endif