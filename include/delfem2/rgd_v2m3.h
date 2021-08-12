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

  //! velocity of shape for deforming rigid body (1-way deformation given by the user)
  std::vector<CVec2d> shape_velo;

  // below: set at the beginning. derived from shape
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

class CContactInfo2 {
public:
  //! index of rigid body A
  unsigned int irbA;

  //! index of rigid body B
  unsigned int irbB;

  //! contacting positioin of rigid body A (corner point)
  unsigned int ipA;

  //! contacting positioin of rigid body A (corner point)
  unsigned int ieB;

  //! ratio on the edge
  double reB;

  //! normal of rigid body B at the contacting point
  CVec2d Njn;

  //! magnitude of impulse
  double lambda;
};

void Steptime_Rgd2(
    std::vector<CRigidState2>& aRS,
    double dt,
    const CVec2d& gravity);

}

#ifndef DFM2_STATIC_LIBRARY
  #include "delfem2/rgd_v2m3.cpp"
#endif

#endif
