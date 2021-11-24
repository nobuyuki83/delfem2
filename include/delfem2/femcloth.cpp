/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/femcloth.h"

#include <cmath>


// --------------------------------------------------------


// compute energy and its 1st and 2nd derivative for contact against object
DFM2_INLINE void delfem2::WdWddW_Contact(
    double& W,  // (out) energy
    double dW[3], // (out) 1st derivative of energy
    double ddW[3][3], // (out) 2nd derivative of energy
    //
    const double c[3], // (in) deformed triangle vertex positions
    double stiff_contact,
    double contact_clearance,
    const CInput_Contact& input )
{
  double n[3];
  double pd = input.penetrationNormal(n[0],n[1],n[2], c[0],c[1],c[2]);
  pd += contact_clearance;
  if( pd  < 0 ){
    W = 0;
    dW[0] = 0;  dW[1] = 0;  dW[2] = 0;
    ddW[0][0] = 0;  ddW[0][1] = 0;  ddW[0][2] = 0;
    ddW[1][0] = 0;  ddW[1][1] = 0;  ddW[1][2] = 0;
    ddW[2][0] = 0;  ddW[2][1] = 0;  ddW[2][2] = 0;
    return;
  }
  W = 0.5*stiff_contact*pd*pd;
  
  dW[0] = -stiff_contact*pd*n[0];
  dW[1] = -stiff_contact*pd*n[1];
  dW[2] = -stiff_contact*pd*n[2];
  
  ddW[0][0] = stiff_contact*n[0]*n[0];
  ddW[0][1] = stiff_contact*n[0]*n[1];
  ddW[0][2] = stiff_contact*n[0]*n[2];
  ddW[1][0] = stiff_contact*n[1]*n[0];
  ddW[1][1] = stiff_contact*n[1]*n[1];
  ddW[1][2] = stiff_contact*n[1]*n[2];
  ddW[2][0] = stiff_contact*n[2]*n[0];
  ddW[2][1] = stiff_contact*n[2]*n[1];
  ddW[2][2] = stiff_contact*n[2]*n[2];
}


