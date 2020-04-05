/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief SMPL model
 * @details skinning
 */

#include <cstdlib>
#include <random>
#include <ctime>
#include <GLFW/glfw3.h>

#include "delfem2/vec3.h"
#include "delfem2/mat3.h"
#include "delfem2/mat4.h"
#include "delfem2/quat.h"
#include "delfem2/mats.h"
#include "delfem2/mshtopo.h"
#include "delfem2/vecxitrsol.h"
//
#include "delfem2/v23m34q.h"
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/v3q_glold.h"
//
#include "delfem2/cnpy/smpl_cnpy.h"

#include "delfem2/rig_v3q.h"
#include "delfem2/opengl/rig_v3m3q_glold.h"

#include "delfem2/opengl/glfw/viewer_glfw.h"

namespace dfm2 = delfem2;

// -------------------


void Solve_MinRigging
(std::vector<dfm2::CRigBone>& aBone,
 const std::vector<dfm2::CTarget>& aTarget)
{
  std::vector<double> Lx, Ly, Lz;
  for(int ibs=0;ibs<aBone.size();++ibs){
    for(int idims=0;idims<3;++idims){
      dfm2::Rig_SensitivityBoneTransform_Eigen(Lx,Ly,Lz,
                                               ibs,idims,true,
                                               aBone);
    }
  }
  for(int idims=0;idims<3;++idims){
    dfm2::Rig_SensitivityBoneTransform_Eigen(Lx,Ly,Lz,
                                             0,idims,false,
                                             aBone);
  }
  
  std::vector<double> aC0; // [nC]
  std::vector<double> adC0; // [nC, nb*3 ]
  for(int it=0;it<aTarget.size();++it){
    dfm2::Rig_WdW_Target_Eigen(aC0,adC0,
                               aBone,aTarget[it],Lx,Ly,Lz);
  }
  
  const unsigned int nsns = Lx.size()/(aBone.size()*4);
  const unsigned int nC = aC0.size();
  
  class CSystemMatrix{
  public:
    CSystemMatrix(const std::vector<double>& adC_,
                  unsigned int nC_,
                  unsigned int nsns_) :
    adC(adC_), nC(nC_), nsns(nsns_)
    {
//      std::cout << "constructor reduced system matrix " << std::endl;
      assert(adC.size()==nsns*nC);
      tmpC0.resize(nC_);
    }
  public:
    void MatVec(double* y,
                double alpha, const double* x, double beta) const {
      dfm2::MatVec(tmpC0.data(),
                    adC.data(), nC, nsns, x);
      dfm2::MatTVec(y,
                   adC.data(), nC, nsns, tmpC0.data());
      for(int i=0;i<nsns;++i){ y[i] += (beta+0.01)*x[i]; }
    }
  public:
    const std::vector<double>& adC;
    unsigned int nC;
    unsigned int nsns;
    mutable std::vector<double> tmpC0;
  } mat(adC0, nC, nsns);
  
  std::vector<double> r(nsns,0.0);
  dfm2::MatTVec(r.data(),
                adC0.data(), nC, nsns, aC0.data());
  
  std::vector<double> u(nsns,0.0);
  std::vector<double> reshist = dfm2::Solve_CG(r.data(), u.data(),
                                               nsns, 1.0e-3, 100, mat);
//  std::cout << "convergence" << reshist.size() << std::endl;
  for(int ib=0;ib<aBone.size();++ib){
    dfm2::CVec3d vec_rot(u.data()+ib*3);
    dfm2::CQuatd dq = dfm2::Quat_CartesianAngle(-vec_rot);
    dfm2::CQuatd q0 = dq*dfm2::CQuatd(aBone[ib].quatRelativeRot);
    q0.CopyTo(aBone[ib].quatRelativeRot);
  }
  {
    dfm2::CVec3d vec_trans(u.data()+aBone.size()*3);
    aBone[0].transRelative[0] -= vec_trans.x();
    aBone[0].transRelative[1] -= vec_trans.y();
    aBone[0].transRelative[2] -= vec_trans.z();
  }
  dfm2::UpdateBoneRotTrans(aBone);
}

void Draw
(const std::vector<double>& aXYZ1,
 const std::vector<unsigned int>& aTri,
 const std::vector<dfm2::CRigBone>& aBone,
 const std::vector<dfm2::CTarget>& aTarget)
{
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_DEPTH_TEST);
  dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), aTri.data(), aTri.size()/3);
  //    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
  { // draw bone
    ::glDisable(GL_DEPTH_TEST);
    ::glDisable(GL_LIGHTING);
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    for(int it=0;it<aTarget.size();++it){
      const unsigned int ib = aTarget[it].ib;
      ::glColor3d(1,0,0);
      dfm2::opengl::myGlVertex(aBone[ib].Pos());
      ::glColor3d(1,0,0);
      dfm2::opengl::myGlVertex(aTarget[it].pos);
    }
    ::glEnd();
  }
  /*
   ::glDisable(GL_DEPTH_TEST);
   delfem2::opengl::DrawBone(aBone,
   -1, -1,
   0.01, 1.0);
   */
  ::glEnable(GL_DEPTH_TEST);
  //    dfm2::opengl::Draw_QuaternionsCoordinateAxes(aXYZ1,aQuat1,0.02);
}

// --------------------

int main()
{
  std::vector<double> aXYZ0;
  std::vector<double> aW;
  std::vector<dfm2::CRigBone> aBone;
  std::vector<unsigned int> aTri;
  {
    std::vector<int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl(aXYZ0,
                         aW,
                         aTri,
                         aIndBoneParent,
                         aJntRgrs,
                         std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    dfm2::Smpl2Rig(aBone,
                   aIndBoneParent, aXYZ0, aJntRgrs);
    
  }
  
  std::vector<double> aXYZ1 = aXYZ0;
  { // initalize pose
    /*
    for(int ibone=0;ibone<aBone.size();++ibone){
      dfm2::CQuatd::Random(0.2).CopyTo(aBone[ibone].quatRelativeRot);
    }
     */
    dfm2::UpdateBoneRotTrans(aBone);
    dfm2::Skinning_LBS(aXYZ1,
                       aXYZ0, aBone, aW);
  }
  
  std::vector<dfm2::CTarget> aTarget;
  {
    {
      dfm2::CTarget t;
      t.ib = 20;
      t.pos = aBone[t.ib].Pos();
      aTarget.push_back(t);
    }
    {
      dfm2::CTarget t;
      t.ib = 10;
      t.pos = aBone[t.ib].Pos();
      aTarget.push_back(t);
    }
  }
  std::vector< dfm2::CVec3d > aTargetOriginPos;
  for(int it=0;it<aTarget.size();++it){
    aTargetOriginPos.push_back(aTarget[it].pos);
  }
      
  // -----------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.camera_rot_mode = dfm2::CAMERA_ROT_TBALL;
  dfm2::opengl::setSomeLighting();

  int iframe = 0;
  while (true)
  {
    iframe++;
    {
      aTarget[0].pos = aTargetOriginPos[0] + 0.4*dfm2::CVec3d(1-cos(0.1*iframe), sin(0.1*iframe), 0.0);
      aTarget[1].pos = aTargetOriginPos[1] + 0.1*dfm2::CVec3d(sin(0.1*iframe), 1-cos(0.1*iframe), 0.0);
      Solve_MinRigging(aBone, aTarget);
      Skinning_LBS(aXYZ1,
                   aXYZ0, aBone, aW);
    }
    
    // -------------------
    viewer.DrawBegin_oldGL();
    Draw(aXYZ1,aTri,aBone,aTarget);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
