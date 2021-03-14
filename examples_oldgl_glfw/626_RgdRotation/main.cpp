/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

#include "delfem2/opengl/glfw/viewer3.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/mshuni.h"
#include "delfem2/opengl/old/color.h"
#include "delfem2/geo3_v23m34q.h"
#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mshprimitive.h"
#include <GLFW/glfw3.h>

namespace dfm2 = delfem2;

// ------------------------------------------------------

namespace delfem2 {

class CRigidBodyState
{
public:
  CVec3d pos;
  CMat3d R;
  CVec3d velo;
  CVec3d Omega;
public:
  CRigidBodyState Step(double dt, const std::vector<CVec3d>& vOpA) const {
    CRigidBodyState rb_out;
    rb_out.velo    = velo  + dt*vOpA[0];
    rb_out.Omega   = Omega + dt*vOpA[1];
    rb_out.pos     = pos   + dt*vOpA[2];
    CMat3d dR = dfm2::Mat3_RotCartesian(dt*vOpA[3]);
    if( dR.isNaN() ) dR.setZero();
    rb_out.R  = R*dR;
    return rb_out;
  }
};

class CRigidBodyInertia
{
public:
  double mass = 1.0;
  CMat3d Irot;
  CMat3d invIrot;
};

class CRigidBodyForceModel
{
public:
  void GetForceTorque(CVec3d& F, CVec3d& T) const{
    F.SetZero();
    T.SetZero();
  }
};

std::vector<CVec3d> VelocityRigidBody(
    const CRigidBodyState& rbs,
    const CRigidBodyInertia& rbi,
    const CRigidBodyForceModel& rbfm)
{
  CVec3d F,T;
  rbfm.GetForceTorque(F,T);
  std::vector<CVec3d> V(4);
  V[0] = (rbs.R*F)*(1.0/rbi.mass);
  V[1] = rbi.invIrot*( (rbs.Omega^(rbi.Irot*rbs.Omega)) + T);
  V[2] = rbs.velo;
  V[3] = rbs.Omega;
  return V;
}

CRigidBodyState StepTime_ForwardEuler(
    double dt,
    const CRigidBodyState& rbIn, // current rigid body
    const CRigidBodyInertia& rbInertia,
    const CRigidBodyForceModel& rbForceModel)
{
  const std::vector<CVec3d>& velo_vOpA = VelocityRigidBody(rbIn, rbInertia, rbForceModel);
  return rbIn.Step(dt, velo_vOpA);
}

CRigidBodyState StepTime_RungeKutta4(
    double dt,
    const CRigidBodyState& rb0, // current rigid body
    const CRigidBodyInertia& rbInertia,
    const CRigidBodyForceModel& rbForceModel)
{
  const std::vector<CVec3d>& vrb1 = VelocityRigidBody(rb0, rbInertia, rbForceModel);
  const CRigidBodyState& rb1 = rb0.Step(dt*0.5, vrb1);
  const std::vector<CVec3d>& vrb2 = VelocityRigidBody(rb1, rbInertia, rbForceModel);
  const CRigidBodyState& rb2 = rb0.Step(dt*0.5, vrb2);
  const std::vector<CVec3d>& vrb3 = VelocityRigidBody(rb2, rbInertia, rbForceModel);
  const CRigidBodyState& rb3 = rb0.Step(dt*1.0, vrb3);
  const std::vector<CVec3d>& vrb4 = VelocityRigidBody(rb3, rbInertia, rbForceModel);
  std::vector<CVec3d> vrb1234(4);
  vrb1234[0] = vrb1[0]+2*vrb2[0]+2*vrb3[0]+vrb4[0];
  vrb1234[1] = vrb1[1]+2*vrb2[1]+2*vrb3[1]+vrb4[1];
  vrb1234[2] = vrb1[2]+2*vrb2[2]+2*vrb3[2]+vrb4[2];
  vrb1234[3] = vrb1[3]+2*vrb2[3]+2*vrb3[3]+vrb4[3];
  return rb0.Step(dt/6.0, vrb1234);
}

}

// ------------------------------------------------

bool is_animation = true;

double dt;
dfm2::CRigidBodyState rbs;
dfm2::CRigidBodyInertia rbi;
dfm2::CRigidBodyForceModel rbfm;

// ------------------------------------------------


void myGlutDisplay(
    const std::vector<double>& aXYZ,
    const std::vector<unsigned int>& aTri)
{
  delfem2::opengl::DrawBackground();
  
  ::glMatrixMode(GL_MODELVIEW);
  ::glPushMatrix();
  {
    double mMV[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    rbs.R.AffineMatrixTrans(mMV);
    ::glMultMatrixd(mMV);
  }
  ::glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
  ::glDisable(GL_LIGHTING);
  ::glColor3d(0,0,0);
  dfm2::opengl::DrawMeshTri3D_Edge(aXYZ.data(),aXYZ.size()/3, aTri.data(), aTri.size()/3);
  ::glEnable(GL_LIGHTING);
  ::glColor3d(1,1,1);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ.data(), aTri.data(), aTri.size()/3);
  ::glPopMatrix();
}

int main(int argc,char* argv[])
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  dfm2::MeshTri3_Torus(aXYZ,aTri, 1.0, 0.2, 32, 8);

  rbi.mass = 1.0;
  {
    rbi.Irot = dfm2::CMat3d::Zero();
    dfm2::CVec3d ex(1,0,0), ey(0,1,0), ez(0,0,1);
    rbi.Irot += 1.0*dfm2::Mat3_OuterProduct(ex,ex);
    rbi.Irot += 3.0*dfm2::Mat3_OuterProduct(ey,ey);
    rbi.Irot += 5.0*dfm2::Mat3_OuterProduct(ez,ez);
  }
  rbi.invIrot = rbi.Irot.Inverse();
  
  rbs.pos = dfm2::CVec3d(0,0,0);
  rbs.R = dfm2::CMat3d::Identity();
  rbs.velo = dfm2::CVec3d(0,0,0);
  rbs.Omega = dfm2::CVec3d(1,1,1);
  
  dt = 0.05;
  
  // ---------------
  dfm2::opengl::CViewer3 viewer;
  viewer.Init_oldGL();
  viewer.camera.view_height = 1.5;
  viewer.camera.camera_rot_mode = delfem2::CCam3_OnAxisZplusLookOrigin<double>::CAMERA_ROT_MODE::TBALL;
  delfem2::opengl::setSomeLighting();

  rbs = StepTime_RungeKutta4(dt, rbs, rbi, rbfm);
  
  while(true){
    rbs = StepTime_RungeKutta4(dt,rbs,rbi,rbfm);
    viewer.DrawBegin_oldGL();
    myGlutDisplay(aXYZ,aTri);
    viewer.SwapBuffers();
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto CLOSE; }
  }
CLOSE:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);

  return 0;
}


