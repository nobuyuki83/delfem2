//
//  inputs_garments.h
//  examples_glfwold_static
//
//  Created by Nobuyuki Umetani on 2020-04-06.
//

#ifndef inputs_garments_h
#define inputs_garments_h

#include "delfem2/garment.h"

namespace delfem2{

void Inputs_SmplTshirt
(std::string& name_cad_in_test_input,
 double& scale_adjust,
 std::vector<unsigned int>& aIESeam,
 double& mesher_edge_length,
 std::vector<CRigidTrans_2DTo3D>& aRT23)
{
  name_cad_in_test_input = "tshirt.svg";
  mesher_edge_length = 0.015;
  scale_adjust = 1.7;
  aIESeam = {
    15, 6,
    13, 0,
    4, 9,
    11, 2,
    20, 17,
    22, 25,
    1, 18,
    12, 19,
    5, 24,
    8, 23
  };
  { // initial position
    aRT23.resize(4);
    { // back body
      CRigidTrans_2DTo3D& rt23 = aRT23[0];
      rt23.org2 = CVec2d(0.189,-0.5)*scale_adjust;
      rt23.org3 = CVec3d(0.0, 0.1, -0.2);
      rt23.R.SetRotMatrix_Cartesian(0.0, 3.1415, 0.0);
    }
    { // front body
      CRigidTrans_2DTo3D& rt23 = aRT23[1];
      rt23.org2 = CVec2d(0.506,-0.5)*scale_adjust;
      rt23.org3 = CVec3d(0.0, 0.1, +0.2);
      rt23.R.SetIdentity();
    }
    { // front body
      CRigidTrans_2DTo3D& rt23 = aRT23[2];
      rt23.org2 = CVec2d(0.833,-0.45)*scale_adjust;
      rt23.org3 = CVec3d(+0.3, 0.3, +0.0);
      rt23.R.SetRotMatrix_BryantAngle(-M_PI*0.5, +M_PI*0.5, 0.0);
      rt23.radinv_x = 13;
    }
    { // front body
      CRigidTrans_2DTo3D& rt23 = aRT23[3];
      rt23.org2 = CVec2d(1.148,-0.45)*scale_adjust;
      rt23.org3 = CVec3d(-0.3, 0.3, +0.0);
      rt23.R.SetRotMatrix_BryantAngle(-M_PI*0.5, -M_PI*0.5, 0.0);
      rt23.radinv_x = 13;
    }
  }
}

void Inputs_SmplLtshirt
(std::string& name_cad_in_test_input,
 double& scale_adjust,
 std::vector<unsigned int>& aIESeam,
 double& mesher_edge_length,
 std::vector<CRigidTrans_2DTo3D>& aRT23)
{
  name_cad_in_test_input = "ltshirt.svg";
  scale_adjust = 1.7;
  aIESeam = {
    15, 6,
    13, 0,
    4, 9,
    11, 2,
    20, 17,
    22, 25,
    1, 18,
    12, 19,
    5, 24,
    8, 23
  };
  mesher_edge_length = 0.018;
  { // initial position
    aRT23.resize(4);
    { // back body
      CRigidTrans_2DTo3D& rt23 = aRT23[0];
      rt23.org2 = CVec2d(0.189-0.04553,-0.5+0.34009)*scale_adjust;
      rt23.org3 = CVec3d(0.0, 0.1, -0.15);
      rt23.R.SetRotMatrix_Cartesian(0.0, 3.1415, 0.0);
    }
    { // front body
      CRigidTrans_2DTo3D& rt23 = aRT23[1];
      rt23.org2 = CVec2d(0.506-0.04553,-0.5+0.34009)*scale_adjust;
      rt23.org3 = CVec3d(0.0, 0.1, +0.2);
      rt23.R.SetIdentity();
    }
    { // front body
      CRigidTrans_2DTo3D& rt23 = aRT23[2];
      rt23.org2 = CVec2d(0.833-0.04553,-0.45+0.34009)*scale_adjust;
      rt23.org3 = CVec3d(+0.3, 0.3, +0.0);
      rt23.R.SetRotMatrix_BryantAngle(-M_PI*0.5, +M_PI*0.5, -0.1);
      rt23.radinv_x = 13;
    }
    { // front body
      CRigidTrans_2DTo3D& rt23 = aRT23[3];
      rt23.org2 = CVec2d(1.148-0.04553,-0.45+0.34009)*scale_adjust;
      rt23.org3 = CVec3d(-0.3, 0.3, +0.0);
      rt23.R.SetRotMatrix_BryantAngle(-M_PI*0.5, -M_PI*0.5, +0.1);
      rt23.radinv_x = 13;
    }
  }
}

}

#endif /* inputs_garments_h */
