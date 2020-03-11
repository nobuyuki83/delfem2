/*
* Copyright (c) 2019 Nobuyuki Umetani
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

/**
 * @brief SMPL model
 * @details hogehoge
 */

#include <cstdlib>
#include <random>
#include "cnpy/cnpy.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/quat.h"
#include "delfem2/mat3.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"

namespace dfm2 = delfem2;


void DrawJoints(
    const std::vector<double>& aJntPos,
    const std::vector<int>& aIndBoneParent)
{
  for(int ib=0;ib<aJntPos.size()/3;++ib){
    const double* p = aJntPos.data()+ib*3;
    ::glColor3d(1,0,0);
    ::glDisable(GL_LIGHTING);
    ::glDisable(GL_DEPTH_TEST);
    dfm2::opengl::DrawSphereAt(8, 8, 0.01, p[0], p[1], p[2]);
    int ibp = aIndBoneParent[ib];
    if( ibp == -1 ){ continue; }
    const double* pp = aJntPos.data()+ibp*3;
    ::glColor3d(0,0,0);
    ::glBegin(GL_LINES);
    ::glVertex3dv(p);
    ::glVertex3dv(pp);
    ::glEnd();
  }
}

int main()
{
  std::vector<double> aXYZ0;
  std::vector<double> aW;
  std::vector<unsigned int> aTri;
  std::vector<int> aIndBoneParent;
  std::vector<double> aJntRgrs;
  {
    cnpy::npz_t my_npz = cnpy::npz_load(std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    const unsigned int nP = my_npz["vertices_template"].shape[0];
    assert( my_npz["vertices_template"].shape[1] == 3 );
    aXYZ0 = my_npz["vertices_template"].as_vec<double>();
    {
      cnpy::NpyArray& npT = my_npz["face_indices"];
      assert( my_npz["face_indices"].shape[1] == 3 );
      aTri = npT.as_vec<unsigned>();
      for(auto &i: aTri){ i -= 1; }
    }
    const unsigned int nBone = my_npz["weights"].shape[1];
    assert( my_npz["weights"].shape[0] == nP );
    aW = my_npz["weights"].as_vec<double>();
    {
      const cnpy::NpyArray& npT = my_npz["kinematic_tree"];
      const int* tree = npT.data<int>();
      aIndBoneParent.assign(tree, tree+nBone);
    }
    assert( my_npz["joint_regressor"].shape[0] == nBone );
    assert( my_npz["joint_regressor"].shape[1] == nP );
    std::cout << my_npz["joint_regressor"].fortran_order << std::endl;
    aJntRgrs = my_npz["joint_regressor"].as_vec<double>();
    assert( aJntRgrs.size() == nBone*nP );
  }
    
  std::vector<double> aJntPos0;
  {
    const unsigned int nP = aXYZ0.size()/3;
    const unsigned int nBone = aIndBoneParent.size();
    aJntPos0.assign(nBone*3, 0.0);
    for(int ib=0;ib<nBone;++ib){
      aJntPos0[ib*3+0] = 0;
      aJntPos0[ib*3+1] = 0;
      aJntPos0[ib*3+2] = 0;
      for(int ip=0;ip<nP;++ip){
        aJntPos0[ib*3+0] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+0];
        aJntPos0[ib*3+1] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+1];
        aJntPos0[ib*3+2] += aJntRgrs[ip*nBone+ib]*aXYZ0[ip*3+2];
      }
      std::cout <<  ib << " ";
      std::cout << aJntPos0[ib*3+0] << " ";
      std::cout << aJntPos0[ib*3+1] << " ";
      std::cout << aJntPos0[ib*3+2] << std::endl;
    }
  }
  
  // ----------------------------------
  
  std::vector<double> aMat4AffineBone; // global affine matrix of a bone
  {
    const unsigned int nbone = aIndBoneParent.size();
    aMat4AffineBone.resize(nbone*16);
    for(unsigned int ibone=0;ibone<nbone;++ibone){
      dfm2::Mat4_Identity(aMat4AffineBone.data()+ibone*16);
    }
  }
  std::vector<double> aQuatRelativeRot;
  {
    const unsigned int nbone = aIndBoneParent.size();
    aQuatRelativeRot.resize(nbone*4);
    for(unsigned int ibone=0;ibone<nbone;++ibone){
      dfm2::Quat_Identity(aQuatRelativeRot.data()+ibone*4);
    }
  }
  std::vector<double> aXYZ1 = aXYZ0;
  std::vector<double> aJntPos1 = aJntPos0;
    
  // -----------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  dfm2::opengl::setSomeLighting();

  while (true)
  {
    {
      std::random_device rd;
      std::mt19937 mt(rd());
      std::uniform_real_distribution<double> dist(-0.1,+0.1);
      const unsigned int nBone = aIndBoneParent.size();
      for(int ibone=0;ibone<nBone;++ibone){
        double* q = aQuatRelativeRot.data()+ibone*4;
        q[0] = 1.0;
        q[1] = dist(mt);
        q[2] = dist(mt);
        q[3] = dist(mt);
        dfm2::Normalize_Quat(q);
      }
      const double trans_root[3] = {dist(mt), dist(mt),dist(mt)};
      dfm2::Mat4_ScaleRotTrans(aMat4AffineBone.data(),
                               1.0, aQuatRelativeRot.data(), trans_root);
      for(int ibone=1;ibone<nBone;++ibone){
        int ibp = aIndBoneParent[ibone];
        assert( ibp >= 0 && ibp < nBone );
        double p1[3];
        dfm2::Vec3_AffMat3Vec3Projection(p1,
                                         aMat4AffineBone.data()+ibp*16, aJntPos0.data()+ibone*3);
        double m0[16];
        dfm2::Mat4_AffineTranslation(m0,
                                     -p1[0], -p1[1], -p1[2]);
        double m1[16];
        dfm2::Mat4_Quat(m1,
                        aQuatRelativeRot.data()+ibone*4);
        double m2[16];
        dfm2::Mat4_AffineTranslation(m2,
                                     +p1[0], +p1[1], +p1[2]);
        double m3[16];
        dfm2::MatMat4(m3,
                      m1,m0);
        double m4[16];
        dfm2::MatMat4(m4,
                      m2,m3);
        dfm2::MatMat4(aMat4AffineBone.data()+ibone*16,
                      aMat4AffineBone.data()+ibp*16, m4);
      }
    }
    for(int ibone=0;ibone<aIndBoneParent.size();++ibone){
      dfm2::Vec3_AffMat3Vec3Projection(aJntPos1.data()+ibone*3,
                                       aMat4AffineBone.data()+ibone*16, aJntPos0.data()+ibone*3);
    }
    {
      const unsigned int nbone = aIndBoneParent.size();
      for(int ip=0;ip<aXYZ0.size()/3;++ip){
        const double* p0 = aXYZ0.data()+ip*3;
        double* p1 = aXYZ1.data()+ip*3;
        p1[0] = 0.0;  p1[1] = 0.0;  p1[2] = 0.0;
        for(int ibone=0;ibone<nbone;++ibone){
          double p2[3];
          dfm2::Vec3_AffMat3Vec3Projection(p2,
                                           aMat4AffineBone.data()+ibone*16, p0);
          p1[0] += aW[ip*nbone+ibone]*p2[0];
          p1[1] += aW[ip*nbone+ibone]*p2[1];
          p1[2] += aW[ip*nbone+ibone]*p2[2];
        }
      }
    }
    
    // -------------------
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    ::glEnable(GL_DEPTH_TEST);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), aTri.data(), aTri.size()/3);
    DrawJoints(aJntPos1, aIndBoneParent);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
