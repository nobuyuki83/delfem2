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
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/v3q_glold.h"
//
#include "delfem2/cnpy/smpl_cnpy.h"

#include "delfem2/rig_v3q.h"
#include "delfem2/opengl/rig_v3m3q_glold.h"

#include "delfem2/opengl/glfw/viewer_glfw.h"

#include "delfem2/opengl/tex_gl.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace dfm2 = delfem2;

// -------------------


void Solve_MinRigging
(std::vector<dfm2::CRigBone>& aBone,
 const std::vector<dfm2::CTarget>& aTarget)
{
  std::vector<double> Lx, Ly, Lz; // [nsns, nb*4]
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
  
  // -----------------------
  
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
  ::glDisable(GL_DEPTH_TEST);
  ::glDisable(GL_LIGHTING);
  ::glPointSize(10);
  ::glBegin(GL_POINTS);
  for(int it=0;it<aTarget.size();++it){
    const unsigned int ib = aTarget[it].ib;
    ::glColor3d(1,0,0);
    dfm2::opengl::myGlVertex(aBone[ib].Pos());
  }
  ::glEnd();
  // ------
  ::glEnable(GL_DEPTH_TEST);
  ::glBegin(GL_LINES);
  ::glColor3d(1,0,0);
  for(int it=0;it<aTarget.size();++it){
    dfm2::CVec3d p = aTarget[it].pos;
    dfm2::opengl::myGlVertex(p+10*dfm2::CVec3d(0,0,1));
    dfm2::opengl::myGlVertex(p-10*dfm2::CVec3d(0,0,1));
  }
  ::glEnd();
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
  dfm2::opengl::CTexRGB_Rect2D tex;
  {
    int width, height, channels;
    unsigned char *img = stbi_load((std::string(PATH_INPUT_DIR)+"/uglysweater.jpg").c_str(),
                                   &width, &height, &channels, 0);
    std::cout << width << " " << height << " " << channels << std::endl;
    tex.Initialize(width, height, img, "rgb");
    delete[] img;
    double scale = 0.5;
    tex.max_x = -scale;
    tex.min_x = +scale;
    tex.max_y = -scale*height/width;
    tex.min_y = +scale*height/width;
    tex.z = -0.5;
  }
  
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
    { // hip right
      dfm2::CTarget t;
      t.ib = 2;
      t.pos = aBone[t.ib].Pos();
      aTarget.push_back(t);
    }
    { // shoulder left
      dfm2::CTarget t;
      t.ib = 16;
      t.pos = dfm2::CVec3d(0.15, +0.25 , 0);
      aTarget.push_back(t);
    }
    { // shoulder right
      dfm2::CTarget t;
      t.ib = 17;
      t.pos = aBone[t.ib].Pos();
      aTarget.push_back(t);
    }
    { // elbow left
      dfm2::CTarget t;
      t.ib = 18;
      t.pos = dfm2::CVec3d(0.2, -0.0, 0);
      aTarget.push_back(t);
    }
    { // wrist left
      dfm2::CTarget t;
      t.ib = 20;
      t.pos = dfm2::CVec3d(0.28, 0.15, 0);
      aTarget.push_back(t);
    }
    { // wrist right
      dfm2::CTarget t;
      t.ib = 21;
      t.pos =  dfm2::CVec3d(-0.18, -0.22, 0);
      aTarget.push_back(t);
    }
  }
  std::vector< std::pair<dfm2::CVec3d,dfm2::CVec3d> > aTargetOriginPos;
  for(int it=0;it<aTarget.size();++it){
    unsigned int ib = aTarget[it].ib;
    aTargetOriginPos.push_back( std::make_pair(aTarget[it].pos,
                                               aBone[ib].Pos()) );
  }
  
  // -----------
//  Check(aBone, aTarget);
     
  // -----------
  dfm2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  viewer.nav.camera.camera_rot_mode = dfm2::CAMERA_ROT_YTOP;
  dfm2::opengl::setSomeLighting();
  tex.InitGL();

  int iframe = 0;
  while (true)
  {
    {
      double r = iframe*0.01;
      if( r > 1 ){ r = 1; }
      for(int it=0;it<aTarget.size();++it){
        aTarget[it].pos = r*aTargetOriginPos[it].first + (1-r)*aTargetOriginPos[it].second;
      }
      Solve_MinRigging(aBone, aTarget);
      Skinning_LBS(aXYZ1,
                   aXYZ0, aBone, aW);
    }
    if( iframe > 200 ){
      for(int ib=0;ib<aBone.size();++ib){
        dfm2::Quat_Identity(aBone[ib].quatRelativeRot);
      }
      dfm2::UpdateBoneRotTrans(aBone);
      for(int it=0;it<aTarget.size();++it){
        aTarget[it].pos = aTargetOriginPos[it].second;
      }
      Solve_MinRigging(aBone, aTarget);
      Skinning_LBS(aXYZ1,
                   aXYZ0, aBone, aW);
      iframe = 0;
    }
    
    // -------------------
    viewer.DrawBegin_oldGL();
    Draw(aXYZ1,aTri,aBone,aTarget);
    tex.Draw_oldGL();
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
    iframe++;
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
