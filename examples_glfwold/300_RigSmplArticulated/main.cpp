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
#include "delfem2/rig_geo3.h"
//
#include "delfem2/cnpy/smpl_cnpy.h"

#include <GLFW/glfw3.h>
#include "delfem2/opengl/funcs_glold.h"
#include "delfem2/opengl/rigv3_glold.h"
#include "delfem2/opengl/glfw/viewer_glfw.h"

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
      for(int ibone=0;ibone<aBone.size();++ibone){
        dfm2::CQuatd::Random(0.2).CopyTo(aBone[ibone].quatRelativeRot);
      }
      dfm2::CVec3d::Random().CopyToScale(aBone[0].transRelative, 0.2);
      dfm2::UpdateBoneRotTrans(aBone);
      dfm2::Skinning_LBS(aXYZ1,
                         aXYZ0, aBone, aW);
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
