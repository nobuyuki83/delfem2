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
#include "delfem2/vecxitrsol.h"
#include "delfem2/quat.h"
#include "delfem2/mat3.h"
//
#include "delfem2/rig_v3q.h"
//
#include "delfem2/cnpy/smpl_cnpy.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/glfw_viewer.h"
#include "delfem2/opengl/glold_funcs.h"
#include "delfem2/opengl/glold_rig_v23q.h"

namespace dfm2 = delfem2;


int main()
{
  std::vector<double> aXYZ0;
  std::vector<double> aW;
  std::vector<unsigned int> aTri;
  std::vector<dfm2::CRigBone> aBone;
  {
    std::vector<int> aIndBoneParent;
    std::vector<double> aJntRgrs;
    dfm2::cnpy::LoadSmpl(aXYZ0,
                         aW,
                         aTri,
                         aIndBoneParent,
                         aJntRgrs,
                         std::string(PATH_INPUT_DIR)+"/smpl_model_f.npz");
    Smpl2Rig(aBone,
             aIndBoneParent, aXYZ0, aJntRgrs);
    
  }
  std::vector<double> aXYZ1 = aXYZ0;
    
  // -----------
  delfem2::opengl::CViewer_GLFW viewer;
  viewer.Init_oldGL();
  dfm2::opengl::setSomeLighting();

  int iframe = 0;
  while (true)
  {
    iframe = (iframe+1)%50;
    if( iframe ==0 ){
      std::random_device rd;
      std::mt19937 mt(rd());
      std::uniform_real_distribution<double> dist(-0.2,+0.2);
      const unsigned int nBone = aBone.size();
      for(int ibone=0;ibone<nBone;++ibone){
        double* q = aBone[ibone].quatRelativeRot;
        q[0] = 1.0;
        q[1] = dist(mt);
        q[2] = dist(mt);
        q[3] = dist(mt);
        dfm2::Normalize_Quat(q);
      }
      aBone[0].transRelative[0] = dist(mt);
      aBone[0].transRelative[1] = dist(mt);
      aBone[0].transRelative[2] = dist(mt);
      dfm2::UpdateBoneRotTrans(aBone);
      {
        for(int ip=0;ip<aXYZ0.size()/3;++ip){
          const double* p0 = aXYZ0.data()+ip*3;
          double* p1 = aXYZ1.data()+ip*3;
          p1[0] = 0.0;  p1[1] = 0.0;  p1[2] = 0.0;
          for(int ibone=0;ibone<nBone;++ibone){
            double p2[3];
            aBone[ibone].DeformSkin(p2, p0);
            p1[0] += aW[ip*nBone+ibone]*p2[0];
            p1[1] += aW[ip*nBone+ibone]*p2[1];
            p1[2] += aW[ip*nBone+ibone]*p2[2];
          }
        }
      }
    }
    
    // -------------------
    viewer.DrawBegin_oldGL();
    ::glEnable(GL_LIGHTING);
    ::glEnable(GL_DEPTH_TEST);
    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ1.data(), aTri.data(), aTri.size()/3);
    ::glDisable(GL_DEPTH_TEST);
    delfem2::opengl::DrawBone(aBone,
                              -1, -1,
                              0.01, 1.0);
//    dfm2::opengl::DrawJoints(aJntPos1, aIndBoneParent);
//    dfm2::opengl::DrawMeshTri3D_FaceNorm(aXYZ0.data(), aTri.data(), aTri.size()/3);
//    dfm2::opengl::DrawJoints(aJntPos0, aIndBoneParent);
    glfwSwapBuffers(viewer.window);
    glfwPollEvents();
    if( glfwWindowShouldClose(viewer.window) ){ goto EXIT; }
  }
EXIT:
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
